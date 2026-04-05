import itertools
import importlib
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyBigWig
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
ROCCO_VERSION = importlib.import_module("rocco._version").__version__


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
        if penalized > best_val or (
            np.isclose(penalized, best_val) and np.sum(sol) < best_count
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
                fwd.query_qualities = pysam.qualitystring_to_array(
                    "I" * int(read_length)
                )
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
                rev.query_qualities = pysam.qualitystring_to_array(
                    "I" * int(read_length)
                )
                bam_file.write(rev)
    pysam.sort("-o", str(bam_path), str(unsorted_path))
    os.remove(unsorted_path)
    pysam.index(str(bam_path))


def _write_toy_bigwig(
    bigwig_path: Path,
    chrom_size: int = 500,
    entries: list[tuple[int, int, float]] | None = None,
):
    if entries is None:
        entries = [
            (0, 50, 0.0),
            (50, 100, 1.5),
            (100, 150, -0.5),
            (150, 200, 2.0),
        ]
    bw = pyBigWig.open(str(bigwig_path), "w")
    try:
        bw.addHeader([("chr1", int(chrom_size))])
        bw.addEntries(
            ["chr1"] * len(entries),
            [int(start) for start, _, _ in entries],
            ends=[int(end) for _, end, _ in entries],
            values=[float(value) for _, _, value in entries],
        )
    finally:
        bw.close()


def _load_bed_records(bed_file: str) -> list[tuple[str, int, int]]:
    records = []
    with open(bed_file, "r", encoding="utf-8") as handle:
        for line in handle:
            line_ = line.strip()
            if line_ == "":
                continue
            chrom, start, end = line_.split("\t")[0:3]
            records.append((str(chrom), int(start), int(end)))
    return records


def _interval_jaccard(
    left_records: list[tuple[str, int, int]],
    right_records: list[tuple[str, int, int]],
) -> float:
    chroms = sorted(
        set(chrom for chrom, _, _ in left_records)
        | set(chrom for chrom, _, _ in right_records)
    )
    overlap_total = 0
    union_total = 0
    for chrom in chroms:
        left = sorted(
            [(start, end) for chrom_, start, end in left_records if chrom_ == chrom]
        )
        right = sorted(
            [(start, end) for chrom_, start, end in right_records if chrom_ == chrom]
        )
        i = 0
        j = 0
        while i < len(left) and j < len(right):
            start = max(left[i][0], right[j][0])
            end = min(left[i][1], right[j][1])
            if end > start:
                overlap_total += end - start
            if left[i][1] <= right[j][1]:
                i += 1
            else:
                j += 1
        chrom_union = sum(end - start for start, end in left) + sum(
            end - start for start, end in right
        )
        union_total += chrom_union
    union_total -= overlap_total
    if union_total <= 0:
        return 0.0
    return float(overlap_total) / float(union_total)


@pytest.mark.correctness
def test_combine_chrom_results_no_names(test_setup):
    combined_outfile = combine_chrom_results(
        [str(x) for x in test_setup["chrom_ref_results"].values()],
        output_file="test_combined.bed",
    )
    assert os.path.exists(combined_outfile)
    combined_jaccard = round(
        _interval_jaccard(
            _load_bed_records(combined_outfile),
            _load_bed_records(test_setup["combined_ref_file"]),
        ),
        5,
    )
    assert combined_jaccard > 0.99
    os.remove(combined_outfile)


@pytest.mark.correctness
def test_score_loci_wls_log_scales_input():
    scores, details = score_loci_wls(
        np.array([[1.0, 15.0]]),
        lower_bound_z=0.0,
        return_details=True,
    )
    assert details["input_scale"] == "log2p1"
    assert "sample_intercepts" not in details
    assert "sample_baselines" not in details
    assert np.allclose(details["mean"], np.array([-1.5, 1.5]))
    assert np.allclose(details["z_scores"], np.array([-0.67449076, 0.67449076]))
    assert np.allclose(scores, np.array([-0.67449076, 0.67449076]))


@pytest.mark.correctness
def test_score_loci_wls_explicit_min_effect_shrinks_standardized_score():
    scores, details = score_loci_wls(
        np.array([[1.0, 15.0]]),
        min_effect=0.5,
        return_details=True,
    )
    assert np.isclose(details["min_effect"], 0.5)
    assert scores[1] < details["z_scores"][1]
    assert scores[0] < details["z_scores"][0]


@pytest.mark.correctness
def test_score_loci_wls_precision_floor_raises_standard_error():
    centered = np.array(
        [
            [0.05, 1.0, 1.0, 0.05],
            [0.04, 1.0, 1.0, 0.04],
            [0.06, 1.0, 1.0, 0.06],
        ],
        dtype=np.float64,
    )
    low_floor_scores, low_floor_details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        prior_df=6.0,
        precision_floor_ratio=0.0,
    )
    high_floor_scores, high_floor_details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        prior_df=6.0,
        precision_floor_ratio=0.25,
    )
    assert np.isclose(high_floor_details["precision_floor_ratio"], 0.25)
    assert np.all(
        high_floor_details["standard_error"] >= low_floor_details["standard_error"]
    )
    assert np.all(high_floor_scores <= low_floor_scores)


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
def test_native_wls_scores_tied_large_matrix():
    if ROCCO_INFERENCE._wls_native is None:
        pytest.skip("native _wls backend is not built in this environment")
    centered = np.zeros((3, 250000), dtype=np.float64)
    scores, details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        lower_bound_z=1.0,
        prior_df=5.0,
    )
    assert scores.shape == (250000,)
    assert np.allclose(details["mean"], 0.0)
    assert np.allclose(details["z_scores"], 0.0)
    assert np.allclose(scores, -1.0)
    assert np.all(details["standard_error"] > 0.0)


@pytest.mark.correctness
def test_native_wls_downweights_noisy_track_locally():
    if ROCCO_INFERENCE._wls_native is None:
        pytest.skip("native _wls backend is not built in this environment")
    x = np.linspace(-4.0, 4.0, 513, dtype=np.float64)
    smooth = 0.9 * np.sin(x) + 0.15 * np.cos(2.0 * x)
    noisy = smooth.copy()
    noisy_region = slice(180, 333)
    noisy[noisy_region] += 0.75 * np.where(
        (np.arange(noisy_region.stop - noisy_region.start) % 2) == 0,
        1.0,
        -1.0,
    )
    centered = np.vstack([smooth, noisy])
    _, details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        lower_bound_z=0.0,
        prior_df=6.0,
        spatial_window=31,
    )
    simple_mean = centered.mean(axis=0)
    quiet_region = slice(40, 140)
    assert np.mean(
        np.abs(details["mean"][noisy_region] - smooth[noisy_region])
    ) < np.mean(np.abs(simple_mean[noisy_region] - smooth[noisy_region]))
    assert np.mean(details["standard_error"][noisy_region]) > np.mean(
        details["standard_error"][quiet_region]
    )


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
    raw = {
        chrom: candidate_counts[chrom] / total_counts[chrom]
        for chrom in candidate_counts
    }
    assert meta["prior_strength"] > 0
    assert meta["prior_dispersion"] >= meta["min_prior_dispersion"]
    assert meta["posterior_summary"] == "beta_quantile"
    assert np.isclose(meta["posterior_quantile"], 0.01)
    assert budgets["chr1"] < raw["chr1"]
    assert budgets["chr2"] < raw["chr2"]
    assert budgets["chr1"] < budgets["chr3"] < budgets["chr2"]


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
def test_empirical_bayes_budget_single_chrom_uses_default_center():
    budgets, meta = estimate_empirical_bayes_budgets(
        {"chr1": 0},
        {"chr1": 0},
    )
    assert np.isclose(meta["genome_wide_budget"], 0.05)
    assert meta["posterior_summary"] == "beta_quantile"
    assert np.isclose(meta["posterior_quantile"], 0.01)
    assert 0.0 < budgets["chr1"] < meta["genome_wide_budget"]


@pytest.mark.correctness
def test_empirical_bayes_budget_lower_quantile_is_more_conservative():
    candidate_counts = {"chr1": 5, "chr2": 50, "chr3": 15, "chr4": 30}
    total_counts = {"chr1": 1000, "chr2": 1000, "chr3": 1000, "chr4": 1000}
    conservative, _ = estimate_empirical_bayes_budgets(
        candidate_counts,
        total_counts,
        posterior_quantile=0.20,
    )
    less_conservative, _ = estimate_empirical_bayes_budgets(
        candidate_counts,
        total_counts,
        posterior_quantile=0.40,
    )
    for chrom in candidate_counts:
        assert conservative[chrom] <= less_conservative[chrom]


@pytest.mark.correctness
def test_context_size_and_gamma_track_peak_width():
    x = np.arange(512, dtype=float)
    narrow = 4.0 * np.exp(-0.5 * ((x - 140.0) / 4.0) ** 2) + 3.5 * np.exp(
        -0.5 * ((x - 340.0) / 5.0) ** 2
    )
    wide = 4.0 * np.exp(-0.5 * ((x - 140.0) / 11.0) ** 2) + 3.5 * np.exp(
        -0.5 * ((x - 340.0) / 13.0) ** 2
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
def test_context_size_is_stable_under_local_ripple():
    x = np.arange(512, dtype=float)
    wide = 4.0 * np.exp(-0.5 * ((x - 140.0) / 11.0) ** 2) + 3.5 * np.exp(
        -0.5 * ((x - 340.0) / 13.0) ** 2
    )
    ripple = 0.35 * np.sin(x / 1.7) * np.exp(-0.5 * ((x - 340.0) / 18.0) ** 2)
    wide_noisy = np.clip(wide + ripple, 0.0, None)

    clean_ctx = estimate_context_size(wide)
    noisy_ctx = estimate_context_size(wide_noisy)

    assert noisy_ctx[0] > 0
    assert abs(noisy_ctx[0] - clean_ctx[0]) <= 0.35 * clean_ctx[0]


@pytest.mark.correctness
def test_length_bins_do_not_get_finer_than_100bp():
    lengths = np.arange(50, 275, 25, dtype=np.int64)
    binned, representatives = ROCCO_SCORES._assign_length_bins(lengths, max_bins=24)

    assert binned.shape == lengths.shape
    assert representatives.size <= 2


@pytest.mark.correctness
def test_build_chrom_cache_partially_pools_context_sizes(monkeypatch):
    chrom_lengths = {"chr_small": 120, "chr_big": 240}
    context_calls = []
    gamma_calls = []
    context_edges = []
    gamma_edges = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        n = chrom_lengths[chrom]
        return np.arange(n, dtype=float), np.zeros((n, 2), dtype=float)

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        n = chrom_matrix.shape[0]
        scores = np.linspace(0.0, 3.0, n, dtype=float)
        return scores, {
            "centered_matrix": np.zeros((n, 2), dtype=float),
            "local_baseline_window": 101,
            "mean": np.linspace(10.0, 13.0, n, dtype=float),
        }

    def fake_budget_estimator(centered_matrix, observed_scores, **kwargs):
        return 0.05, {"effective_total_count": float(observed_scores.shape[0])}

    def fake_estimate_context_size(vals, **kwargs):
        n = int(np.asarray(vals).shape[0])
        context_calls.append(n)
        context_edges.append((float(np.asarray(vals)[0]), float(np.asarray(vals)[-1])))
        if n <= 120:
            point, lower, upper = 6, 4, 8
        else:
            point, lower, upper = 24, 20, 28
        return (
            point,
            lower,
            upper,
            {
                "feature_indices": np.array([10, 20], dtype=np.int64),
                "feature_scores_log": np.array([1.0, 0.8], dtype=float),
                "num_features": 2,
                "tau_sq_hat": 0.1,
                "context_log_mean": float(np.log(point)),
                "context_log_mean_var": 0.5,
                "context_max_span": 64,
                "best_order": 1,
            },
        )

    def fake_estimate_gamma_from_scores(scores, **kwargs):
        gamma_calls.append(kwargs.get("context_size_summary"))
        gamma_edges.append(
            (float(np.asarray(scores)[0]), float(np.asarray(scores)[-1]))
        )
        n = int(np.asarray(scores).shape[0])
        return float(n), {
            "gamma_scale": 0.5,
            "signal_scale": 1.0,
            "context_size_point": 12,
            "context_size_lower": 10,
            "context_size_upper": 14,
            "num_features": 2,
            "feature_detection_order": 1,
        }

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
        "input_track_type": "bam",
        "score_lower_bound_z": 1.0,
        "score_prior_df": 5.0,
        "score_precision_floor_ratio": 0.01,
        "budget_null_draws": 4,
        "gamma": None,
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

    assert context_calls == [120, 240]
    assert all(np.allclose(edge, (10.0, 13.0)) for edge in context_edges)
    assert gamma_calls[0] != gamma_calls[1]
    assert all(np.allclose(edge, (10.0, 13.0)) for edge in gamma_edges)
    assert (
        chrom_cache["chr_small"]["gamma_meta"]["context_pooling"] == "log_width_partial"
    )
    assert (
        chrom_cache["chr_big"]["gamma_meta"]["context_pooling"] == "log_width_partial"
    )
    assert chrom_cache["chr_small"]["gamma_meta"]["raw_context_size_point"] == 6
    assert chrom_cache["chr_big"]["gamma_meta"]["raw_context_size_point"] == 24
    assert 6 < chrom_cache["chr_small"]["gamma_meta"]["context_size_point"] < 24
    assert 6 < chrom_cache["chr_big"]["gamma_meta"]["context_size_point"] < 24


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
def test_generate_chrom_matrix_reads_bigwig_scores_directly(tmp_path):
    bw_path = tmp_path / "toy.bw"
    chrom_sizes_path = tmp_path / "toy_bw.sizes"
    chrom_sizes_path.write_text("chr1\t500\n", encoding="utf-8")
    _write_toy_bigwig(bw_path)

    intervals, score_matrix = generate_chrom_matrix(
        "chr1",
        [str(bw_path)],
        str(chrom_sizes_path),
        step=999,
        round_digits=6,
    )

    assert intervals.tolist() == [0, 50, 100, 150]
    assert score_matrix.shape == (1, 4)
    assert np.allclose(score_matrix[0], np.array([0.0, 1.5, -0.5, 2.0]))


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
        "chr1\t90\t110\n" "chr1\t100\t150\n" "chr1\t150\t210\n",
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
def test_build_chrom_cache_uses_bigwig_scores_directly(monkeypatch):
    direct_budget_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        return np.array([0, 50, 100, 150], dtype=int), np.array(
            [
                [0.0, 2.0, 1.0, 0.0],
                [0.0, 3.0, 2.0, 0.0],
            ],
            dtype=float,
        )

    def fail_score_loci_wls(*args, **kwargs):
        raise AssertionError("bigWig inputs should bypass WLS scoring")

    def fake_budget_estimator(scores, **kwargs):
        direct_budget_calls.append(np.asarray(scores, dtype=float))
        return 0.05, {"effective_total_count": float(len(scores))}

    def fake_estimate_context_size(vals, **kwargs):
        return 6, 4, 8, {
            "feature_indices": np.array([1, 2], dtype=np.int64),
            "feature_scores_log": np.array([1.0, 0.8], dtype=float),
            "num_features": 2,
            "tau_sq_hat": 0.1,
            "context_log_mean": float(np.log(6.0)),
            "context_log_mean_var": 0.5,
            "context_max_span": 64,
            "best_order": 1,
        }

    def fake_estimate_gamma_from_scores(scores, **kwargs):
        return 3.0, {
            "gamma_scale": float(kwargs["gamma_scale"]),
            "signal_scale": 1.0,
            "context_size_point": int(kwargs["context_size_summary"][0]),
            "context_size_lower": int(kwargs["context_size_summary"][1]),
            "context_size_upper": int(kwargs["context_size_summary"][2]),
            "num_features": 2,
            "feature_detection_order": 1,
        }

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fail_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_budget_nonnull_fraction_from_score_track",
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

    args = {
        "chrom_sizes_file": None,
        "step": 50,
        "round_digits": 5,
        "effective_genome_size": None,
        "norm_method": "RPGC",
        "min_mapping_score": 0,
        "flag_include": None,
        "flag_exclude": None,
        "extend_reads": 0,
        "center_reads": False,
        "ignore_for_norm": [],
        "scale_factor": 1.0,
        "threads": 1,
        "input_track_type": "bigwig",
        "score_lower_bound_z": 1.0,
        "score_prior_df": 5.0,
        "score_min_effect": None,
        "score_precision_floor_ratio": 0.01,
        "budget_null_draws": 4,
        "gamma": None,
        "gamma_scale": 0.5,
        "gamma_context_min_span": 3,
        "gamma_context_max_span": 64,
        "gamma_context_band_z": 1.0,
        "gamma_context_max_order": 5,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr1"],
        ["track1.bw", "track2.bw"],
        args,
        None,
    )

    assert len(direct_budget_calls) == 1
    assert np.allclose(direct_budget_calls[0], np.array([0.0, 2.5, 1.5, 0.0]))
    assert np.allclose(chrom_cache["chr1"]["scores"], np.array([0.0, 2.5, 1.5, 0.0]))
    assert chrom_cache["chr1"]["gamma"] == 3.0


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
def test_version_flag():
    result = subprocess.run(
        [sys.executable, "-m", "rocco.rocco", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode == 0
    assert result.stdout.decode().strip() == f"rocco {ROCCO_VERSION}"


@pytest.mark.correctness
def test_no_input_listed():
    result = subprocess.run(
        ["rocco", "--input_files"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    assert result.returncode != 0
    assert "usage:" in result.stderr.decode()


@pytest.mark.correctness
def test_unrecognized_arg():
    result = subprocess.run(
        ["rocco", "--unrecognized_arg"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    assert result.returncode != 0
    assert "unrecognized" in result.stderr.decode()
