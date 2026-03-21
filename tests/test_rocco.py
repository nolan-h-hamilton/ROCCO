import itertools
import importlib
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pybedtools as pbt
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import rocco as rocco_module
from rocco import *

TEST_DIR = Path(__file__).resolve().parent
ROCCO_IMPL = importlib.import_module("rocco.rocco")
ROCCO_INFERENCE = importlib.import_module("rocco.inference")
ROCCO_READTRACKS = importlib.import_module("rocco.readtracks")
ROCCO_SCORES = importlib.import_module("rocco.scores")


@pytest.fixture
def test_setup():
    chromosomes = ["chr19", "chr21", "chrX"]
    chrom_ref_results = {
        "chr19": str(TEST_DIR / "ref_chr19.bed"),
        "chr21": str(TEST_DIR / "ref_chr21.bed"),
        "chrX": str(TEST_DIR / "ref_chrX.bed"),
    }
    matrices = np.load(TEST_DIR / "test_data.npz")
    intervals = np.load(TEST_DIR / "test_intervals.npz")
    return {
        "chromosomes": chromosomes,
        "chrom_ref_results": chrom_ref_results,
        "combined_ref_file": str(TEST_DIR / "combined_ref.bed"),
        "matrices": matrices,
        "intervals": intervals,
    }


def _bruteforce_penalized(scores, switch_costs, selection_penalty):
    scores = np.asarray(scores, dtype=float)
    switch_costs = np.asarray(switch_costs, dtype=float)
    n = len(scores)
    best_sol = None
    best_val = -np.inf
    best_count = None
    for bits in itertools.product([0, 1], repeat=n):
        sol = np.asarray(bits, dtype=np.uint8)
        penalized = (
            scores @ sol
            - np.sum(switch_costs * np.abs(np.diff(sol)))
            - (selection_penalty * np.sum(sol))
        )
        if (
            penalized > best_val
            or (
                np.isclose(penalized, best_val)
                and np.sum(sol) < best_count
            )
        ):
            best_sol = sol
            best_val = float(penalized)
            best_count = int(np.sum(sol))
    return best_sol, best_val, best_count


def _write_toy_bam(bam_path: Path, chrom_size: int = 500):
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "chr1", "LN": int(chrom_size)}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam_file:
        for name, start in [("r1", 100), ("r2", 125), ("r3", 160)]:
            segment = pysam.AlignedSegment()
            segment.query_name = name
            segment.query_sequence = "A" * 50
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = int(start)
            segment.mapping_quality = 60
            segment.cigar = ((0, 50),)
            segment.query_qualities = pysam.qualitystring_to_array("I" * 50)
            bam_file.write(segment)
    pysam.index(str(bam_path))


def _write_single_end_fragment_bam(
    bam_path: Path,
    fragment_length: int = 80,
    read_length: int = 30,
    chrom_size: int = 12000,
):
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "chr1", "LN": int(chrom_size)}],
    }
    unsorted_path = bam_path.with_name(f"{bam_path.stem}.unsorted.bam")
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as bam_file:
        for block_start in [200, 4200, 8200]:
            for offset in range(0, 120, 8):
                fwd = pysam.AlignedSegment()
                fwd.query_name = f"f{block_start}_{offset}"
                fwd.query_sequence = "A" * int(read_length)
                fwd.flag = 0
                fwd.reference_id = 0
                fwd.reference_start = int(block_start + offset)
                fwd.mapping_quality = 60
                fwd.cigar = ((0, int(read_length)),)
                fwd.query_qualities = pysam.qualitystring_to_array("I" * int(read_length))
                bam_file.write(fwd)

                rev = pysam.AlignedSegment()
                rev.query_name = f"r{block_start}_{offset}"
                rev.query_sequence = "T" * int(read_length)
                rev.flag = 16
                rev.reference_id = 0
                rev.reference_start = int(
                    block_start + offset + fragment_length - read_length
                )
                rev.mapping_quality = 60
                rev.cigar = ((0, int(read_length)),)
                rev.query_qualities = pysam.qualitystring_to_array("I" * int(read_length))
                bam_file.write(rev)
    pysam.sort("-o", str(bam_path), str(unsorted_path))
    os.remove(unsorted_path)
    pysam.index(str(bam_path))


@pytest.mark.correctness
def test_combine_chrom_results_no_names(test_setup):
    combined_outfile = combine_chrom_results(
        [str(x) for x in test_setup["chrom_ref_results"].values()],
        output_file="test_combined.bed",
    )
    assert os.path.exists(combined_outfile)
    combined_pbt = pbt.BedTool(combined_outfile)
    combined_ref_pbt = pbt.BedTool(test_setup["combined_ref_file"])
    combined_jaccard = round(
        combined_pbt.jaccard(combined_ref_pbt)["jaccard"], 5
    )
    assert combined_jaccard > 0.99
    os.remove(combined_outfile)


@pytest.mark.correctness
def test_score_central_tendency_chrom():
    X = np.array([[1, 2, 3, 4, 5], [2, 3, 4, 5, 6], [3, 4, 5, 6, 10]])
    scores = score_central_tendency_chrom(X, method="quantile")
    assert np.allclose(scores, np.array([2, 3, 4, 5, 6]))

    scores = score_central_tendency_chrom(X, method="mean")
    assert np.allclose(scores, np.array([2, 3, 4, 5, 7]))

    X = np.array(
        [
            [0, 0, 0, 0, 0],
            [2, 3, 4, 5, 6],
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 6],
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 6],
            [1, 2, 3, 4, 5],
            [5, 3, 4, 5, 6],
            [1, 2, 3, 4, 5],
            [1000, 1000, 1000, 1000, 1000],
        ]
    )
    scores = score_central_tendency_chrom(X, method="tmean", tprop=0.11)
    assert np.allclose(scores, np.array([1.875, 2.5, 3.5, 4.5, 5.5]))


@pytest.mark.correctness
def test_score_dispersion_chrom():
    X = np.zeros(shape=(11, 3))
    X[:, 0] = [0, 1, 1, 1, 1, 5, 25, 25, 25, 25, 100]
    X[:, 1] = np.ones(11)
    X[:, 2] = [1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1]

    scores = score_dispersion_chrom(X, method="mad")
    assert np.allclose(scores, np.array([5, 0, 0]))

    scores = score_dispersion_chrom(X, method="std")
    assert np.allclose(scores, np.array([27.89265, 0, 0.481045]), rtol=1e-4, atol=1e-4)

    scores = score_dispersion_chrom(X, method="iqr")
    assert np.allclose(scores, np.array([24, 0, 1]))


@pytest.mark.correctness
def test_score_loci_wls_emphasizes_consistent_signal():
    X = np.array(
        [
            [1.0, 1.1, 8.0, 8.1],
            [0.9, 1.0, 7.9, 8.2],
            [1.1, 0.95, 8.2, 7.8],
        ]
    )
    scores, details = score_loci_wls(X, return_details=True)
    assert scores.shape == (4,)
    assert np.all(details["sample_weights"] > 0)
    assert np.mean(scores[2:]) > np.mean(scores[:2])
    assert np.mean(details["z_scores"][2:]) > np.mean(details["z_scores"][:2])


@pytest.mark.correctness
def test_score_loci_wls_log_scales_input():
    scores, details = score_loci_wls(
        np.array([[1.0, 15.0]]),
        lower_bound_z=0.0,
        return_details=True,
    )
    assert details["input_scale"] == "log2p1"
    assert np.allclose(details["sample_intercepts"], np.array([2.5]))
    assert np.allclose(details["sample_baselines"], np.array([2.5]))
    assert np.allclose(details["mean"], np.array([-1.5, 1.5]))
    assert np.allclose(scores, np.array([-1.5, 1.5]))


@pytest.mark.correctness
def test_score_loci_wls_low_memory_keeps_centered_matrix_in_float32():
    scores, details = score_loci_wls(
        np.array([[1.0, 3.0, 7.0], [1.2, 2.8, 6.5]]),
        low_memory=True,
        return_details=True,
    )
    assert scores.dtype == np.float64
    assert details["centered_matrix"].dtype == np.float32
    assert np.all(np.isfinite(details["centered_matrix"]))


@pytest.mark.correctness
def test_native_wls_matches_python_on_tied_large_matrix():
    centered = np.zeros((3, 250000), dtype=np.float64)
    native_scores, native_details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        lower_bound_z=1.0,
        prior_df=5.0,
    )
    python_scores, python_details = ROCCO_INFERENCE._score_centered_wls_matrix_python(
        centered,
        lower_bound_z=1.0,
        prior_df=5.0,
    )
    assert np.allclose(native_scores, python_scores)
    assert np.allclose(
        native_details["sample_weights"],
        python_details["sample_weights"],
    )
    assert np.allclose(native_details["mean"], python_details["mean"])
    assert np.allclose(
        native_details["standard_error"],
        python_details["standard_error"],
    )
    assert np.allclose(native_details["z_scores"], python_details["z_scores"])


@pytest.mark.correctness
def test_consenrich_crossfit_local_baseline_tracks_broad_background():
    x = np.arange(129, dtype=np.float64)
    broad = 2.5 * np.exp(-0.5 * ((x - 64.0) / 18.0) ** 2)
    spike = 5.0 * np.exp(-0.5 * ((x - 64.0) / 2.5) ** 2)
    y = broad + spike
    baseline = ROCCO_INFERENCE._consenrich_crossfit_whittaker_baseline(
        y,
        block_size=41,
    )
    residual = y - baseline

    shoulder_idx = 46
    peak_idx = 64
    assert baseline.shape == y.shape
    assert baseline[shoulder_idx] > 0.5 * broad[shoulder_idx]
    assert residual[peak_idx] > 3.0 * max(residual[shoulder_idx], 1.0e-6)


@pytest.mark.correctness
def test_exact_dp_matches_bruteforce():
    rng = np.random.default_rng(7)
    scores = rng.normal(size=9)
    switch_costs = rng.uniform(0.2, 1.3, size=8)
    for selection_penalty in (-0.5, 0.0, 0.6, 1.4):
        sol, penalized_objective, selected_count = solve_penalized_chain(
            scores,
            switch_costs,
            selection_penalty,
        )
        brute_sol, brute_val, brute_count = _bruteforce_penalized(
            scores,
            switch_costs,
            selection_penalty,
        )
        assert np.array_equal(sol, brute_sol)
        assert np.isclose(penalized_objective, brute_val)
        assert selected_count == brute_count


@pytest.mark.correctness
def test_solve_chrom_exact_respects_budget():
    scores = np.array([0.5, 1.5, 1.4, -0.2, 3.0, 2.8, -0.1, 0.1])
    solution, objective, details = solve_chrom_exact(
        scores,
        budget=0.375,
        gamma=1.0,
        return_details=True,
    )
    assert solution.dtype == np.uint8
    assert np.sum(solution) <= int(np.floor(len(scores) * 0.375))
    assert np.isclose(
        objective,
        objective_value(
            solution,
            scores,
            build_switch_costs(scores, gamma=1.0),
        ),
    )
    assert details["selected_fraction"] <= 0.375


@pytest.mark.correctness
def test_empirical_bayes_budget_shrinkage():
    candidate_counts = {"chr1": 5, "chr2": 50, "chr3": 15}
    total_counts = {"chr1": 1000, "chr2": 1000, "chr3": 1000}
    budgets, meta = estimate_empirical_bayes_budgets(
        candidate_counts,
        total_counts,
    )
    raw = {chrom: candidate_counts[chrom] / total_counts[chrom] for chrom in candidate_counts}
    assert meta["prior_strength"] > 0
    assert meta["prior_dispersion"] >= meta["min_prior_dispersion"]
    assert meta["posterior_summary"] == "mean_minus_sd"
    assert budgets["chr1"] < raw["chr1"]
    assert budgets["chr2"] < raw["chr2"]
    assert budgets["chr1"] < budgets["chr3"] < budgets["chr2"]


@pytest.mark.correctness
def test_empirical_bayes_budget_dispersion_hits_binomial_floor_when_rates_match():
    candidate_counts = {"chr1": 10, "chr2": 20, "chr3": 30}
    total_counts = {"chr1": 100, "chr2": 200, "chr3": 300}
    budgets, meta = estimate_empirical_bayes_budgets(
        candidate_counts,
        total_counts,
    )
    assert meta["prior_dispersion"] >= meta["min_prior_dispersion"]
    assert meta["prior_dispersion_at_floor"]
    assert meta["posterior_summary"] == "mean_minus_sd"
    assert meta["prior_fit_method"] == "weak_pooled_prior"
    assert meta["observed_raw_budget_var"] <= meta["theoretical_min_raw_budget_var"]
    assert budgets["chr1"] < 0.1
    assert budgets["chr2"] < 0.1
    assert budgets["chr3"] < 0.1
    assert budgets["chr1"] < budgets["chr2"] < budgets["chr3"]


@pytest.mark.correctness
def test_budget_nonnull_fraction_from_empirical_null_reports_soft_count_metadata():
    x = np.arange(512, dtype=np.float64)
    peak1 = 6.0 * np.exp(-0.5 * ((x - 120.0) / 15.0) ** 2)
    peak2 = 5.5 * np.exp(-0.5 * ((x - 320.0) / 15.0) ** 2)
    chrom_matrix = np.vstack(
        [
            0.25 + peak1 + peak2 + 0.05 * np.sin(x / 13.0),
            0.20 + 0.95 * peak1 + 1.05 * peak2 + 0.04 * np.cos(x / 15.0),
            0.22 + 1.1 * peak1 + 0.9 * peak2 + 0.05 * np.sin(x / 17.0),
        ]
    )
    scores, details = score_loci_wls(
        chrom_matrix,
        return_details=True,
    )
    fraction, meta = estimate_budget_nonnull_fraction_from_empirical_null(
        details["centered_matrix"],
        observed_scores=scores,
        dependence_lag_hint=16,
        num_null_draws=6,
        return_details=True,
    )
    assert 0.0 < fraction <= 1.0
    assert np.isclose(fraction, meta["nonnull_fraction"])
    assert 0.0 <= meta["observed_positive_fraction"] <= 1.0
    assert 0.0 <= meta["null_positive_fraction"] <= 1.0
    assert meta["observed_excess_mass"] > meta["null_excess_mass"] > 0.0
    assert meta["observed_excess_units"] > meta["null_excess_units"] > 0.0
    assert meta["effective_count"] > 0.0
    assert 1.0 <= meta["effective_total_count"] <= meta["num_loci"]
    assert meta["autocorrelation_time"] >= 1.0
    assert meta["ess_max_lag"] == 64.0
    assert meta["null_method"] == "dependent_wild_residual_bootstrap"
    assert meta["num_null_draws"] == 6.0
    assert meta["max_null_draws"] == 6.0
    assert not meta["adaptive_stop"]
    assert meta["wild_bandwidth"] >= 8.0
    assert meta["null_excess_units_sd"] > 0.0
    assert meta["null_reference_mean_positive_consensus"] >= 0.0
    assert meta["negative_support_size"] > 0.0
    assert 0.0 < meta["negative_fraction"] <= 1.0


@pytest.mark.correctness
def test_budget_nonnull_fraction_parallel_matches_serial():
    if "fork" not in mp.get_all_start_methods():
        pytest.skip("parallel budget null smoke test requires fork support")

    x = np.arange(512, dtype=np.float64)
    peak = 6.0 * np.exp(-0.5 * ((x - 200.0) / 18.0) ** 2)
    chrom_matrix = np.vstack(
        [
            0.2 + peak + 0.03 * np.sin(x / 19.0),
            0.22 + 0.95 * peak + 0.04 * np.cos(x / 23.0),
            0.18 + 1.05 * peak + 0.03 * np.sin(x / 29.0),
        ]
    )
    scores, details = score_loci_wls(chrom_matrix, return_details=True)
    serial_fraction, serial_meta = estimate_budget_nonnull_fraction_from_empirical_null(
        details["centered_matrix"],
        observed_scores=scores,
        dependence_lag_hint=16,
        num_null_draws=6,
        num_processes=1,
        return_details=True,
    )
    parallel_fraction, parallel_meta = estimate_budget_nonnull_fraction_from_empirical_null(
        details["centered_matrix"],
        observed_scores=scores,
        dependence_lag_hint=16,
        num_null_draws=6,
        num_processes=2,
        return_details=True,
    )
    assert np.isclose(serial_fraction, parallel_fraction)
    assert np.isclose(serial_meta["null_excess_units"], parallel_meta["null_excess_units"])
    assert serial_meta["num_null_draws"] == parallel_meta["num_null_draws"] == 6.0


@pytest.mark.correctness
def test_empirical_bayes_budget_single_chrom_uses_default_center():
    budgets, meta = estimate_empirical_bayes_budgets(
        {"chr1": 0},
        {"chr1": 0},
    )
    assert np.isclose(meta["genome_wide_budget"], 0.05)
    assert meta["posterior_summary"] == "mean_minus_sd"
    assert np.isclose(budgets["chr1"], 0.05)


@pytest.mark.correctness
def test_context_size_and_gamma_track_peak_width():
    x = np.arange(512, dtype=float)
    narrow = (
        4.0 * np.exp(-0.5 * ((x - 140.0) / 4.0) ** 2)
        + 3.5 * np.exp(-0.5 * ((x - 340.0) / 5.0) ** 2)
    )
    wide = (
        4.0 * np.exp(-0.5 * ((x - 140.0) / 11.0) ** 2)
        + 3.5 * np.exp(-0.5 * ((x - 340.0) / 13.0) ** 2)
    )
    narrow_ctx = estimate_context_size(narrow)
    wide_ctx = estimate_context_size(wide)
    narrow_gamma, narrow_meta = estimate_gamma_from_scores(narrow)
    wide_gamma, wide_meta = estimate_gamma_from_scores(wide)

    assert wide_ctx[0] > narrow_ctx[0]
    assert narrow_gamma > 0
    assert wide_gamma > narrow_gamma
    assert wide_meta["context_size_point"] > narrow_meta["context_size_point"]


@pytest.mark.correctness
def test_length_bins_do_not_get_finer_than_100bp():
    lengths = np.arange(50, 275, 25, dtype=np.int64)
    binned, representatives = ROCCO_SCORES._assign_length_bins(lengths, max_bins=24)

    assert binned.shape == lengths.shape
    assert representatives.size <= 2


@pytest.mark.correctness
def test_build_chrom_cache_uses_largest_chrom_for_shared_context(monkeypatch):
    chrom_lengths = {"chr_small": 120, "chr_big": 240}
    context_calls = []
    gamma_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        n = chrom_lengths[chrom]
        return np.arange(n, dtype=float), np.zeros((n, 2), dtype=float)

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        n = chrom_matrix.shape[0]
        scores = np.linspace(0.0, 3.0, n, dtype=float)
        return scores, {
            "centered_matrix": np.zeros((n, 2), dtype=float),
            "local_baseline_window": 101,
        }

    def fake_budget_estimator(centered_matrix, observed_scores, **kwargs):
        return 0.05, {"effective_total_count": float(observed_scores.shape[0])}

    def fake_estimate_context_size(vals, **kwargs):
        context_calls.append(int(np.asarray(vals).shape[0]))
        return 12, 10, 14, {
            "feature_indices": np.array([10, 20], dtype=np.int64),
            "feature_scores_log": np.array([1.0, 0.8], dtype=float),
            "num_features": 2,
            "tau_sq_hat": 0.1,
            "best_order": 1,
        }

    def fake_estimate_gamma_from_scores(scores, **kwargs):
        gamma_calls.append(kwargs.get("context_size_summary"))
        n = int(np.asarray(scores).shape[0])
        return float(n), {
            "gamma_scale": 0.5,
            "signal_scale": 1.0,
            "signal_log_var": 0.1,
            "context_size_point": 12,
            "context_size_lower": 10,
            "context_size_upper": 14,
            "width_log_var": 0.1,
            "log_gamma_var": 0.1,
            "num_features": 2,
            "context_shared": 1,
            "feature_detection_order": 1,
        }

    def fake_shrink_gamma_estimates(gamma_hats, gamma_variances):
        return dict(gamma_hats), {"prior_mean": 0.0, "prior_tau_sq": 0.0}

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fake_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_budget_nonnull_fraction_from_wild_bootstrap_null",
        fake_budget_estimator,
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_context_size",
        fake_estimate_context_size,
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_gamma_from_scores",
        fake_estimate_gamma_from_scores,
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "shrink_gamma_estimates",
        fake_shrink_gamma_estimates,
    )

    args = {
        "chrom_sizes_file": None,
        "step": 50,
        "round_digits": 5,
        "effective_genome_size": None,
        "norm_method": "rpkm",
        "min_mapping_score": 0,
        "flag_include": None,
        "flag_exclude": None,
        "extend_reads": 0,
        "center_reads": False,
        "ignore_for_norm": [],
        "scale_factor": 1.0,
        "threads": 1,
        "score_lower_bound_z": 1.0,
        "score_prior_df": 5.0,
        "budget_null_draws": 4,
        "gamma": None,
        "gamma_strategy": "peak_width",
        "gamma_scale": 0.5,
        "gamma_context_min_span": 3,
        "gamma_context_max_span": 64,
        "gamma_context_band_z": 1.0,
        "gamma_context_max_order": 5,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr_small", "chr_big"],
        [],
        args,
        None,
    )

    assert context_calls == [240]
    assert gamma_calls == [(12, 10, 14), (12, 10, 14)]
    assert chrom_cache["chr_small"]["gamma_meta"]["context_source_chrom"] == "chr_big"
    assert chrom_cache["chr_big"]["gamma_meta"]["context_source_chrom"] == "chr_big"


@pytest.mark.correctness
def test_gamma_shrinkage_pulls_extremes_toward_genome_wide_center():
    shrunk, meta = shrink_gamma_estimates(
        {"chr1": 1.0, "chr2": 100.0, "chr3": 10.0},
        {"chr1": 0.05, "chr2": 0.05, "chr3": 0.05},
    )
    assert meta["num_chromosomes"] == 3
    assert shrunk["chr1"] > 1.0
    assert shrunk["chr2"] < 100.0
    assert shrunk["chr1"] < shrunk["chr3"] < shrunk["chr2"]


@pytest.mark.correctness
def test_empirical_null_survival_has_finite_sample_correction():
    empirical_null = EmpiricalNull(np.array([0.1, 0.2, 0.3]))
    assert np.isclose(empirical_null.survival(-1.0), 1.0)
    assert np.isclose(empirical_null.survival(0.3), 0.5)
    assert np.isclose(empirical_null.survival(10.0), 0.25)


@pytest.mark.correctness
def test_length_binning_reduces_ecdf_count():
    binned_lengths, representatives = ROCCO_SCORES._assign_length_bins(
        np.array([50, 52, 55, 90, 95, 400, 420]),
        max_bins=3,
    )
    assert np.unique(binned_lengths).size <= 3
    assert representatives.size <= 3
    assert binned_lengths[0] == binned_lengths[1]
    assert binned_lengths[-1] == binned_lengths[-2]


@pytest.mark.correctness
def test_generate_chrom_matrix_counts_bam_with_native_backend(tmp_path):
    bam_path = tmp_path / "toy.bam"
    chrom_sizes_path = tmp_path / "toy.sizes"
    chrom_sizes_path.write_text("chr1\t500\n", encoding="utf-8")
    _write_toy_bam(bam_path)

    intervals, count_matrix = generate_chrom_matrix(
        "chr1",
        [str(bam_path)],
        str(chrom_sizes_path),
        step=50,
        effective_genome_size=150,
        norm_method="RPGC",
        flag_exclude=0,
        ignore_for_norm=[],
        round_digits=6,
    )

    assert intervals.tolist() == [100, 150, 200]
    assert count_matrix.shape == (1, 3)
    assert np.allclose(count_matrix[0], np.array([2.0, 2.0, 1.0]))


@pytest.mark.correctness
def test_generate_chrom_matrix_low_memory_uses_float32(tmp_path):
    bam_path = tmp_path / "toy_lowmem.bam"
    chrom_sizes_path = tmp_path / "toy_lowmem.sizes"
    chrom_sizes_path.write_text("chr1\t500\n", encoding="utf-8")
    _write_toy_bam(bam_path)

    _, count_matrix = generate_chrom_matrix(
        "chr1",
        [str(bam_path)],
        str(chrom_sizes_path),
        step=50,
        effective_genome_size=150,
        norm_method="RPGC",
        flag_exclude=0,
        ignore_for_norm=[],
        round_digits=6,
        low_memory=True,
    )

    assert count_matrix.dtype == np.float32
    assert np.all(np.isfinite(count_matrix))


@pytest.mark.correctness
def test_native_fragment_length_estimation_for_single_end_bam(tmp_path):
    bam_path = tmp_path / "single_end_frag.bam"
    _write_single_end_fragment_bam(bam_path, fragment_length=80, read_length=30)

    fragment_length = ROCCO_READTRACKS._hts_counts.get_alignment_fragment_length(
        str(bam_path),
        max_insert_size=200,
        block_size=256,
        rolling_chunk_size=8,
        lag_step=1,
        early_exit=32,
        fallback=0,
    )

    assert 70 <= int(fragment_length) <= 90


@pytest.mark.correctness
def test_single_end_fragment_inference_is_used_for_counting(tmp_path):
    bam_path = tmp_path / "single_end_count.bam"
    _write_single_end_fragment_bam(bam_path, fragment_length=80, read_length=30)

    metadata = ROCCO_READTRACKS._get_bam_count_metadata(
        str(bam_path),
        step=25,
        norm_method="RPGC",
        effective_genome_size=2000,
        ignore_for_norm=[],
        extend_reads=0,
        num_processors=1,
    )

    assert metadata["paired_end"] is False
    assert metadata["paired_end_mode"] is False
    assert int(metadata["read_length"]) == 30
    assert int(metadata["norm_read_length"]) >= 70
    assert int(metadata["resolved_extend_bp"]) >= 70


@pytest.mark.correctness
def test_raw_count_matrix_uses_native_interval_counter(tmp_path):
    bam_path = tmp_path / "toy.bam"
    peak_path = tmp_path / "toy_peaks.bed"
    output_path = tmp_path / "toy_counts.tsv"
    _write_toy_bam(bam_path)
    peak_path.write_text(
        "chr1\t90\t110\n"
        "chr1\t100\t150\n"
        "chr1\t150\t210\n",
        encoding="utf-8",
    )

    count_matrix_path = raw_count_matrix(
        [str(bam_path)],
        str(peak_path),
        str(output_path),
    )
    rows = Path(count_matrix_path).read_text(encoding="utf-8").strip().splitlines()
    assert rows[0] == "peak_name\ttoy"
    assert rows[1] == "chr1_90_110\t1"
    assert rows[2] == "chr1_100_150\t2"
    assert rows[3] == "chr1_150_210\t2"


@pytest.mark.correctness
def test_prepare_args_low_memory_uses_conservative_defaults(monkeypatch):
    parser = ROCCO_IMPL._build_parser()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rocco",
            "-i",
            "fake.bam",
            "-s",
            "fake.sizes",
            "--norm_method",
            "CPM",
            "--low_memory",
        ],
    )
    args = ROCCO_IMPL._prepare_args(parser)
    assert args["low_memory"] is True
    assert 1 <= int(args["threads"]) <= 4
    assert int(args["budget_null_draws"]) == 16


@pytest.mark.correctness
def test_removed_scoring_api_deleted():
    assert not hasattr(rocco_module, "parsig")
    assert not hasattr(rocco_module, "score_chrom_linear")
    assert not hasattr(rocco_module, "apply_filter")
    assert not hasattr(rocco_module, "apply_transformation")
    assert not hasattr(rocco_module, "resolve_transformations")


@pytest.mark.correctness
def test_small_end_to_end_subset(test_setup):
    chrom = "chr21"
    scores, details = score_loci_wls(
        test_setup["matrices"][chrom][:, :5000],
        return_details=True,
    )
    assert scores.shape == (5000,)
    budget_fraction, budget_meta = estimate_budget_nonnull_fraction_from_empirical_null(
        details["centered_matrix"],
        observed_scores=scores,
        dependence_lag_hint=101,
        num_null_draws=6,
        return_details=True,
    )
    assert 0.0 <= budget_fraction <= 1.0
    assert budget_meta["effective_count"] >= 0.0
    assert 1.0 <= budget_meta["effective_total_count"] <= budget_meta["num_loci"]
    candidate_mask = candidate_mask_from_wls(details["z_scores"], tail_z=2.0)
    assert candidate_mask.shape == (5000,)

    solution, objective, _ = solve_chrom_exact(
        scores,
        budget=0.02,
        gamma=1.0,
        return_details=True,
    )
    assert np.sum(solution) <= int(np.floor(5000 * 0.02))
    chrom_outfile = chrom_solution_to_bed(
        chrom,
        test_setup["intervals"][chrom][:5000],
        solution,
        ID="subset",
    )
    assert os.path.exists(chrom_outfile)
    os.remove(chrom_outfile)
    assert np.isfinite(objective)


@pytest.mark.correctness
def test_no_input_no_args():
    result = subprocess.run(["rocco"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode == 0
    assert "usage:" in result.stdout.decode()


@pytest.mark.correctness
def test_no_input_listed():
    result = subprocess.run(["rocco", "--input_files"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode != 0
    assert "usage:" in result.stderr.decode()


@pytest.mark.correctness
def test_unrecognized_arg():
    result = subprocess.run(["rocco", "--unrecognized_arg"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode != 0
    assert "unrecognized" in result.stderr.decode()


@pytest.mark.correctness
def test_bedgraph_input(tmp_path):
    bedgraph_path = tmp_path / "test.bedgraph"
    bedgraph_path.write_text("chr1\t0\t50\t1.0\n", encoding="utf-8")
    chrom_sizes_path = Path(__file__).resolve().parents[1] / "rocco" / "hg38.sizes"
    result = subprocess.run(
        [
            "rocco",
            "--input_files",
            str(bedgraph_path),
            "--chrom_sizes_file",
            str(chrom_sizes_path),
            "--effective_genome_size",
            "2913022398",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode != 0
    assert "bam" in result.stderr.decode().lower()


@pytest.mark.correctness
def test_removed_solver_flag_is_rejected():
    parser = ROCCO_IMPL._build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--solver", "pdlp"])
    assert exc_info.value.code != 0


@pytest.mark.correctness
@pytest.mark.parametrize(
    "flag",
    [
        "--transform_locratio",
        "--transform_logpc",
        "--transform_savgol",
        "--transform_medfilt",
    ],
)
def test_heuristic_transform_flags_are_rejected(flag):
    parser = ROCCO_IMPL._build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([flag])
    assert exc_info.value.code != 0
