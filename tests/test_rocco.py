import itertools
import importlib
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest
try:
    import pyBigWig
except ImportError:  # pragma: no cover - depends on optional extra
    pyBigWig = None

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import rocco as rocco_module
from rocco import *

TEST_DIR = Path(__file__).resolve().parent
ROCCO_IMPL = importlib.import_module("rocco.rocco")
ROCCO_INFERENCE = importlib.import_module("rocco.inference")
ROCCO_READTRACKS = importlib.import_module("rocco.readtracks")
ROCCO_SCORES = importlib.import_module("rocco.scores")
ROCCO_VERSION = importlib.import_module("rocco._version").__version__
BIGWIG_AVAILABLE = pyBigWig is not None


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
        "matrices": matrices,
        "intervals": intervals,
    }


def _bruteforce_penalized(scores, switch_costs, selection_penalty):
    scores = np.asarray(scores, dtype=float)
    switch_costs = np.asarray(switch_costs, dtype=float)
    selection_penalty = np.asarray(selection_penalty, dtype=float)
    n = len(scores)
    best_sol = None
    best_val = -np.inf
    best_count = None
    for bits in itertools.product([0, 1], repeat=n):
        sol = np.asarray(bits, dtype=np.uint8)
        if selection_penalty.ndim == 0:
            selected_penalty = float(selection_penalty) * np.sum(sol)
        else:
            selected_penalty = selection_penalty @ sol
        penalized = (
            scores @ sol
            - np.sum(switch_costs * np.abs(np.diff(sol)))
            - selected_penalty
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
    if pyBigWig is None:
        raise ImportError("pyBigWig is required for bigWig tests")
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
    expected_records = []
    for ref_file in test_setup["chrom_ref_results"].values():
        expected_records.extend(_load_bed_records(ref_file))
    combined_jaccard = round(
        _interval_jaccard(
            _load_bed_records(combined_outfile),
            expected_records,
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
def test_write_narrowpeak_summit_offsets_uses_wls_mean_centers(tmp_path):
    peak_file = tmp_path / "peaks.bed"
    peak_file.write_text(
        "chr1\t100\t250\nchr1\t150\t250\nchr2\t0\t50\n",
        encoding="utf-8",
    )
    summit_track_file = tmp_path / "chr1_summit_track.npz"
    np.savez(
        summit_track_file,
        starts=np.array([100, 150, 200], dtype=np.int64),
        centers=np.array([125, 175, 225], dtype=np.int64),
        mean=np.array([1.0, 5.0, 2.0], dtype=np.float32),
    )
    chrom_cache = {
        "chr1": {"summit_track_file": str(summit_track_file)},
        "chr2": {"summit_track_file": None},
    }
    output_file = tmp_path / "pointsource.tsv"
    ROCCO_IMPL._write_narrowpeak_summit_offsets(
        str(peak_file),
        chrom_cache,
        str(output_file),
    )
    assert output_file.read_text(encoding="utf-8").splitlines() == [
        "chr1_100_250\t75",
        "chr1_150_250\t25",
        "chr2_0_50\t-1",
    ]


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
def test_solve_chrom_exact_soft_budget_modes():
    scores = np.array([0.5, 1.5, 1.4, -0.2, 3.0, 2.8, -0.1, 0.1])
    cases = [
        (None, None, "unpenalized", 0.0, None),
        (0.0, None, "soft_selection_penalty", np.quantile(scores, 1.0), 0),
        (0.25, None, "soft_selection_penalty", np.quantile(scores, 0.75), 2),
        (1.0, None, "soft_selection_penalty", 0.0, 8),
        (np.array([0.0, 1.0]), 0.6, "manual_selection_penalty", 0.6, None),
    ]
    for (
        budget,
        selection_penalty,
        budget_mode,
        expected_penalty,
        expected_target_count,
    ) in cases:
        solution, objective, details = solve_chrom_exact(
            scores,
            budget=budget,
            gamma=1.0,
            selection_penalty=selection_penalty,
            return_details=True,
        )
        assert solution.dtype == np.uint8
        assert details["budget_mode"] == budget_mode
        assert details["soft_budget_penalty"] == pytest.approx(expected_penalty)
        assert details["selection_penalty"] == pytest.approx(expected_penalty)
        assert np.isclose(
            objective,
            objective_value(
                solution,
                scores,
                build_switch_costs(scores, gamma=1.0),
            ),
        )
        if budget is not None and selection_penalty is None:
            assert details["budget"] == pytest.approx(
                float(np.mean(np.atleast_1d(budget)))
            )
            assert details["budget_target_count"] == expected_target_count
        else:
            assert "budget" not in details


@pytest.mark.correctness
def test_solve_chrom_exact_block_budget_preserves_boundary_switch_cost():
    scores = np.array([2.0, 2.0, 3.0, 3.0])
    budget_blocks = np.array([[0.0, 2.0, 0.0], [2.0, 4.0, 1.0]])
    solution, objective, details = solve_chrom_exact(
        scores,
        budget_blocks=budget_blocks,
        gamma=5.0,
        return_details=True,
    )
    penalty_track = np.array([2.0, 2.0, 0.0, 0.0])
    brute_sol, brute_val, brute_count = _bruteforce_penalized(
        scores,
        build_switch_costs(scores, gamma=5.0),
        penalty_track,
    )
    assert np.array_equal(solution, brute_sol)
    assert np.array_equal(solution, np.ones(4, dtype=np.uint8))
    assert np.isclose(details["penalized_objective"], brute_val)
    assert details["selected_count"] == brute_count
    assert np.isclose(
        objective,
        objective_value(solution, scores, build_switch_costs(scores, gamma=5.0)),
    )
    assert details["budget_mode"] == "block_soft_selection_penalty"
    assert details["budget_block_count"] == 2
    assert np.array_equal(details["budget_block_starts"], np.array([0, 2]))
    assert np.array_equal(details["budget_block_ends"], np.array([2, 4]))
    assert np.allclose(details["budget_block_penalties"], np.array([2.0, 0.0]))
    with pytest.raises(ValueError, match="finite integers"):
        solve_chrom_exact(
            scores,
            budget_blocks=np.array([[0.0, 2.5, 0.0], [2.5, 4.0, 1.0]]),
            gamma=5.0,
            return_details=True,
        )

    _, _, manual_details = solve_chrom_exact(
        scores,
        budget_blocks=np.array([[0.0, 2.0, 2.0], [2.0, 4.0, -1.0]]),
        gamma=5.0,
        selection_penalty=2.5,
        return_details=True,
    )
    assert manual_details["budget_mode"] == "manual_selection_penalty"
    assert manual_details["selection_penalty"] == pytest.approx(2.5)
    assert "budget" not in manual_details


@pytest.mark.correctness
@pytest.mark.parametrize(("n_loci", "num_null_blocks"), [(10, 4), (3, 4), (1, 4)])
def test_budget_block_bounds_cover_chromosome(n_loci, num_null_blocks):
    bounds = ROCCO_IMPL._make_null_blocks(n_loci, num_null_blocks)
    starts = np.array([start for start, _ in bounds], dtype=int)
    ends = np.array([end for _, end in bounds], dtype=int)
    lengths = ends - starts

    assert len(bounds) == min(n_loci, num_null_blocks)
    assert starts[0] == 0
    assert ends[-1] == n_loci
    assert np.array_equal(ends[:-1], starts[1:])
    assert np.all(lengths > 0)
    assert int(np.max(lengths) - np.min(lengths)) <= 1


@pytest.mark.correctness
def test_estimate_correlation_length_cases():
    rng = np.random.default_rng(43)
    iid = rng.normal(size=1024)
    ar1 = np.empty(1024, dtype=np.float64)
    noise = rng.normal(size=1024)
    ar1[0] = noise[0]
    for idx in range(1, ar1.size):
        ar1[idx] = (0.9 * ar1[idx - 1]) + noise[idx]

    iid_span, iid_meta = ROCCO_INFERENCE.estimate_correlation_length(iid, step_bp=50)
    ar1_span, ar1_meta = ROCCO_INFERENCE.estimate_correlation_length(ar1, step_bp=50)

    assert iid_span == 25
    assert ar1_span > iid_span
    assert iid_meta["correlation_length_method"] == "acf_crossing"
    assert ar1_meta["correlation_length_intervals"] == ar1_span
    with pytest.raises(ValueError, match="more loci"):
        ROCCO_INFERENCE.estimate_correlation_length(np.arange(16.0), step_bp=50)


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
def test_budget_nonnull_fraction_from_wild_bootstrap_reports_soft_count_metadata():
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
    fraction, meta = estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
        details["centered_matrix"],
        observed_scores=scores,
        correlation_length=16,
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
    assert meta["wild_bandwidth"] == 16.0
    assert meta["dwb_bandwidth"] == 16.0
    assert meta["correlation_length_intervals"] == 16.0
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
def test_resolve_budgets_clips_final_budgets_to_fixed_range():
    chrom_cache = {
        "chr1": {"budget_count_hat": 0.0, "total_count": 1000.0},
        "chr2": {"budget_count_hat": 2.0, "total_count": 1000.0},
        "chr3": {"budget_count_hat": 120.0, "total_count": 1000.0},
        "chr4": {"budget_count_hat": 400.0, "total_count": 1000.0},
    }
    budgets, meta = ROCCO_IMPL._resolve_budgets(
        chrom_cache,
        {
            "budget_posterior_quantile": 0.01,
            "budget": None,
            "scale_chrom_budgets": 1.0,
        },
    )
    assert meta["posterior_summary"] == "beta_quantile"
    for chrom in chrom_cache:
        assert 0.001 <= budgets[chrom] <= 0.25


@pytest.mark.correctness
def test_resolve_budgets_pools_chrom_block_units(monkeypatch):
    chrom_cache = {
        "chr1": {
            "budget_count_hat": 3.0,
            "total_count": 20.0,
            "budget_rate_meta": {
                "budget_blocks": [
                    {
                        "block_id": 0,
                        "start_idx": 0,
                        "stop_idx": 5,
                        "start_bp": 0,
                        "stop_bp": 250,
                        "num_loci": 5,
                        "budget_fraction_hat": 0.1,
                        "budget_count_hat": 1.0,
                        "total_count": 10.0,
                    },
                    {
                        "block_id": 1,
                        "start_idx": 5,
                        "stop_idx": 10,
                        "start_bp": 250,
                        "stop_bp": 500,
                        "num_loci": 5,
                        "budget_fraction_hat": 0.2,
                        "budget_count_hat": 2.0,
                        "total_count": 10.0,
                    },
                ]
            },
        },
        "chr2": {
            "budget_count_hat": 7.0,
            "total_count": 20.0,
            "budget_rate_meta": {
                "budget_blocks": [
                    {
                        "block_id": 0,
                        "start_idx": 0,
                        "stop_idx": 4,
                        "start_bp": 0,
                        "stop_bp": 200,
                        "num_loci": 4,
                        "budget_fraction_hat": 0.3,
                        "budget_count_hat": 3.0,
                        "total_count": 10.0,
                    },
                    {
                        "block_id": 1,
                        "start_idx": 4,
                        "stop_idx": 8,
                        "start_bp": 200,
                        "stop_bp": 400,
                        "num_loci": 4,
                        "budget_fraction_hat": 0.4,
                        "budget_count_hat": 4.0,
                        "total_count": 10.0,
                    },
                ]
            },
        },
    }
    captured = {}

    def fake_estimate_empirical_bayes_budgets(
        candidate_counts,
        total_counts,
        posterior_quantile,
    ):
        captured["candidate_counts"] = dict(candidate_counts)
        captured["total_counts"] = dict(total_counts)
        return (
            {key: 0.01 * (idx + 1) for idx, key in enumerate(candidate_counts)},
            {"genome_wide_budget": 0.1, "posterior_summary": "fake"},
        )

    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_empirical_bayes_budgets",
        fake_estimate_empirical_bayes_budgets,
    )

    budgets, meta = ROCCO_IMPL._resolve_budgets(
        chrom_cache,
        {
            "budget_posterior_quantile": 0.01,
            "budget": None,
            "scale_chrom_budgets": 1.0,
        },
    )

    assert list(captured["candidate_counts"]) == [
        "chr1:0",
        "chr1:1",
        "chr2:0",
        "chr2:1",
    ]
    assert list(captured["candidate_counts"].values()) == [1.0, 2.0, 3.0, 4.0]
    assert list(captured["total_counts"].values()) == [10.0, 10.0, 10.0, 10.0]
    assert budgets["chr1"]["mode"] == "block"
    assert budgets["chr1"]["budget"] == pytest.approx(0.015)
    assert np.allclose(
        budgets["chr1"]["block_bounds"],
        np.array([[0.0, 5.0], [5.0, 10.0]]),
    )
    assert np.allclose(budgets["chr2"]["block_budgets"], np.array([0.03, 0.04]))
    assert np.allclose(
        ROCCO_IMPL._solver_budget_blocks(budgets["chr1"]),
        np.array([[0.0, 5.0, 0.01], [5.0, 10.0, 0.02]]),
    )
    assert budgets["chr2"]["blocks"][1]["budget"] == pytest.approx(0.04)
    assert meta["budget_unit_count"] == 4


@pytest.mark.correctness
def test_solve_cached_chromosome_passes_block_budgets(monkeypatch, tmp_path):
    captured = {}

    def fake_solve_chrom_exact(
        scores,
        budget=None,
        gamma=0.25,
        selection_penalty=None,
        return_details=False,
        budget_blocks=None,
    ):
        captured["scores"] = np.asarray(scores, dtype=float).copy()
        captured["budget"] = budget
        captured["budget_blocks"] = np.asarray(budget_blocks, dtype=float).copy()
        captured["gamma"] = gamma
        captured["selection_penalty"] = selection_penalty
        captured["return_details"] = return_details
        return np.array([1, 0, 1, 0], dtype=np.uint8), 12.0, {
            "budget_mode": "block_soft_selection_penalty",
            "soft_budget_penalty": 1.5,
            "selected_count": 2,
        }

    def fake_chrom_solution_to_bed(*args, **kwargs):
        captured["bed_args"] = args
        return str(tmp_path / "chr1.bed")

    monkeypatch.setattr(ROCCO_IMPL, "solve_chrom_exact", fake_solve_chrom_exact)
    monkeypatch.setattr(
        ROCCO_IMPL,
        "chrom_solution_to_bed",
        fake_chrom_solution_to_bed,
    )

    ROCCO_IMPL._CHROM_SOLVE_PROCESS_STATE = {
        "chrom_cache": {
            "chr1": {
                "scores": np.array([0.0, 1.0, 2.0, 3.0]),
                "intervals": np.array([0, 50, 100, 150]),
                "gamma": 2.0,
            }
        },
        "chrom_budgets": {
            "chr1": {
                "mode": "block",
                "budget": 0.15,
                "block_bounds": np.array([[0.0, 2.0], [2.0, 4.0]]),
                "block_budgets": np.array([0.1, 0.2]),
            }
        },
        "selection_penalty": None,
        "min_length_bp": None,
        "run_id": "test",
    }
    try:
        chrom, objective, meta, solution, outfile = ROCCO_IMPL._solve_cached_chromosome(
            "chr1"
        )
    finally:
        ROCCO_IMPL._CHROM_SOLVE_PROCESS_STATE = None

    assert chrom == "chr1"
    assert objective == pytest.approx(12.0)
    assert meta["budget_mode"] == "block_soft_selection_penalty"
    assert np.array_equal(solution, np.array([1, 0, 1, 0], dtype=np.uint8))
    assert outfile == str(tmp_path / "chr1.bed")
    assert captured["budget"] is None
    assert captured["gamma"] == pytest.approx(2.0)
    assert captured["selection_penalty"] is None
    assert captured["return_details"] is True
    assert np.allclose(
        captured["budget_blocks"],
        np.array([[0.0, 2.0, 0.1], [2.0, 4.0, 0.2]]),
    )


@pytest.mark.correctness
def test_length_bins_do_not_get_finer_than_100bp():
    lengths = np.arange(50, 275, 25, dtype=np.int64)
    binned, representatives = ROCCO_SCORES._assign_length_bins(lengths, max_bins=24)

    assert binned.shape == lengths.shape
    assert representatives.size <= 2


@pytest.mark.correctness
def test_build_chrom_cache_uses_fixed_gamma_and_correlation_length(monkeypatch):
    chrom_lengths = {"chr_small": 120, "chr_big": 240}
    budget_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        n = chrom_lengths[chrom]
        return 50 * np.arange(n, dtype=float), np.zeros((2, n), dtype=float)

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        n = chrom_matrix.shape[1]
        scores = np.linspace(0.0, 3.0, n, dtype=float)
        return scores, {
            "centered_matrix": np.zeros((2, n), dtype=float),
            "mean": np.linspace(10.0, 13.0, n, dtype=float),
            "z_scores": scores + 1.0,
        }

    def fake_estimate_correlation_length(scores, **kwargs):
        return 7, {"correlation_length_intervals": 7}

    def fake_budget_estimator(centered_matrix, observed_scores, **kwargs):
        observed_scores_ = np.asarray(observed_scores, dtype=float)
        observed_len = int(observed_scores_.shape[0])
        budget_calls.append(
            {
                "centered_shape": np.asarray(centered_matrix).shape,
                "observed_scores": observed_scores_.copy(),
                "observed_len": observed_len,
                "kwargs": kwargs,
            }
        )
        return 0.05, {
            "budget_count_hat": 0.05 * float(observed_len),
            "effective_total_count": float(observed_len),
            "dwb_bandwidth": float(kwargs["correlation_length"]),
            "num_null_draws": float(kwargs["num_null_draws"]),
        }

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fake_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_correlation_length",
        fake_estimate_correlation_length,
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_budget_nonnull_fraction_from_wild_bootstrap_null",
        fake_budget_estimator,
    )

    args = {
        "chrom_sizes_file": None,
        "step": 50,
        "round_digits": 5,
        "effective_genome_size": None,
        "norm_method": "RPKM",
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
        "num_null_blocks": 4,
        "gamma": 2.5,
        "peak_mode": None,
        "dependence_span": None,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr_small", "chr_big"],
        [],
        args,
    )

    assert chrom_cache["chr_small"]["gamma"] == 2.5
    assert chrom_cache["chr_big"]["gamma"] == 2.5
    assert chrom_cache["chr_small"]["correlation_length_intervals"] == 7
    assert chrom_cache["chr_big"]["dwb_bandwidth"] == 7
    assert [call["observed_len"] for call in budget_calls] == [
        30,
        30,
        30,
        30,
        60,
        60,
        60,
        60,
    ]
    assert [call["centered_shape"] for call in budget_calls] == [
        (2, 30),
        (2, 30),
        (2, 30),
        (2, 30),
        (2, 60),
        (2, 60),
        (2, 60),
        (2, 60),
    ]
    assert np.allclose(
        np.concatenate([call["observed_scores"] for call in budget_calls[:4]]),
        np.linspace(0.0, 3.0, 120, dtype=float),
    )
    assert np.allclose(
        np.concatenate([call["observed_scores"] for call in budget_calls[4:]]),
        np.linspace(0.0, 3.0, 240, dtype=float),
    )
    assert [call["kwargs"]["correlation_length"] for call in budget_calls] == [7] * 8
    assert [call["kwargs"]["step_bp"] for call in budget_calls] == [50] * 8
    assert [call["kwargs"]["random_seed"] for call in budget_calls] == [
        0,
        1009,
        2018,
        3027,
        0,
        1009,
        2018,
        3027,
    ]
    assert chrom_cache["chr_small"]["budget_rate_meta"]["num_null_blocks"] == 4
    assert chrom_cache["chr_small"]["budget_rate_meta"]["budget_blocks"][0][
        "stop_idx"
    ] == 30
    assert chrom_cache["chr_big"]["total_count"] == pytest.approx(240.0)


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("peak_mode", "expected_summit_calls"),
    [(None, 0), ("narrow", 1), ("broad", 0), ("both", 1)],
)
def test_build_chrom_cache_peak_modes(monkeypatch, peak_mode, expected_summit_calls):
    summit_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        return 50 * np.arange(64, dtype=float), np.zeros((2, 64), dtype=float)

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        scores = np.linspace(-1.0, 2.5, chrom_matrix.shape[1], dtype=float)
        return scores, {
            "centered_matrix": np.zeros_like(chrom_matrix, dtype=float),
            "mean": scores.copy(),
            "z_scores": scores + 2.0,
        }

    def fake_estimate_correlation_length(scores, **kwargs):
        return 5, {"correlation_length_intervals": 5}

    def fake_budget_estimator(centered_matrix, observed_scores, **kwargs):
        observed_len = int(np.asarray(observed_scores).shape[0])
        return 0.05, {
            "budget_count_hat": 0.05 * float(observed_len),
            "effective_total_count": float(observed_len),
            "dwb_bandwidth": 5.0,
            "num_null_draws": float(kwargs["num_null_draws"]),
        }

    def fake_summit_track(chrom, intervals, effect_mean):
        summit_calls.append((chrom, intervals.copy(), effect_mean.copy()))
        return None

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fake_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_correlation_length",
        fake_estimate_correlation_length,
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimate_budget_nonnull_fraction_from_wild_bootstrap_null",
        fake_budget_estimator,
    )
    monkeypatch.setattr(ROCCO_IMPL, "_cpy_narrowpeak_summit_track", fake_summit_track)

    args = {
        "chrom_sizes_file": None,
        "step": 50,
        "round_digits": 5,
        "effective_genome_size": None,
        "norm_method": "RPKM",
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
        "gamma": 0.25,
        "peak_mode": peak_mode,
        "dependence_span": None,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr1"],
        [],
        args,
    )

    assert len(summit_calls) == expected_summit_calls
    assert "effect_mean" not in chrom_cache["chr1"]
    assert chrom_cache["chr1"]["gamma"] == pytest.approx(0.25)
    assert chrom_cache["chr1"]["correlation_length_intervals"] == 5


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
@pytest.mark.skipif(not BIGWIG_AVAILABLE, reason="pyBigWig not installed")
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

    def fake_estimate_correlation_length(scores, **kwargs):
        return 3, {"correlation_length_intervals": 3}

    def fake_budget_estimator(scores, **kwargs):
        scores_ = np.asarray(scores, dtype=float)
        direct_budget_calls.append((scores_, kwargs))
        return 0.05, {
            "budget_count_hat": 0.05 * float(scores_.shape[0]),
            "effective_total_count": float(scores_.shape[0]),
            "dwb_bandwidth": float(kwargs["correlation_length"]),
            "num_null_draws": float(kwargs["num_null_draws"]),
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
        "estimate_correlation_length",
        fake_estimate_correlation_length,
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
        "num_null_blocks": 4,
        "gamma": 3.0,
        "peak_mode": None,
        "dependence_span": None,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr1"],
        ["track1.bw", "track2.bw"],
        args,
    )

    assert len(direct_budget_calls) == 4
    assert [call[0].shape[0] for call in direct_budget_calls] == [1, 1, 1, 1]
    assert np.allclose(
        np.concatenate([call[0] for call in direct_budget_calls]),
        np.array([0.0, 2.5, 1.5, 0.0]),
    )
    assert [call[1]["correlation_length"] for call in direct_budget_calls] == [3] * 4
    assert [call[1]["random_seed"] for call in direct_budget_calls] == [
        0,
        1009,
        2018,
        3027,
    ]
    assert np.allclose(chrom_cache["chr1"]["scores"], np.array([0.0, 2.5, 1.5, 0.0]))
    assert chrom_cache["chr1"]["gamma"] == 3.0
    assert chrom_cache["chr1"]["budget_rate_meta"]["num_null_blocks"] == 4


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
    assert int(args["num_null_blocks"]) == 4


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("extra_args", "message"),
    [
        (["--gamma", "-0.1"], "gamma"),
        (
            [
                "--peak_mode",
                "broad",
                "--score_lower_bound_z",
                "1.0",
                "--broad_score_lower_bound_z",
                "2.0",
            ],
            "cannot exceed",
        ),
        (["--num_null_blocks", "0"], "num_null_blocks"),
        (["-i", "fake.bw", "--peak_mode", "broad"], "BAM inputs"),
    ],
)
def test_prepare_args_validation_failures(monkeypatch, extra_args, message):
    parser = ROCCO_IMPL._build_parser()
    argv = [
        "rocco",
        "-i",
        "fake.bam",
        "-s",
        "fake.sizes",
        "--norm_method",
        "CPM",
    ]
    if extra_args and extra_args[0] == "-i":
        argv = ["rocco", "-s", "fake.sizes"] + extra_args
    else:
        argv.extend(extra_args)
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(ValueError, match=message):
        ROCCO_IMPL._prepare_args(parser)


@pytest.mark.correctness
def test_small_end_to_end_subset(test_setup):
    chrom = "chr21"
    scores, details = score_loci_wls(
        test_setup["matrices"][chrom][:, :5000],
        return_details=True,
    )
    assert scores.shape == (5000,)
    budget_fraction, budget_meta = estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
        details["centered_matrix"],
        observed_scores=scores,
        correlation_length=101,
        num_null_draws=6,
        return_details=True,
    )
    assert 0.0 <= budget_fraction <= 1.0
    assert budget_meta["effective_count"] >= 0.0
    assert 1.0 <= budget_meta["effective_total_count"] <= budget_meta["num_loci"]
    candidate_mask = candidate_mask_from_wls(details["z_scores"], tail_z=2.0)
    assert candidate_mask.shape == (5000,)

    solution, objective, solve_details = solve_chrom_exact(
        scores,
        budget=0.02,
        gamma=1.0,
        return_details=True,
    )
    assert solve_details["budget_mode"] == "soft_selection_penalty"
    assert solve_details["budget"] == pytest.approx(0.02)
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
    result = subprocess.run(
        [sys.executable, "-m", "rocco.rocco"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
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
        [sys.executable, "-m", "rocco.rocco", "--input_files"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode != 0
    assert "usage:" in result.stderr.decode()


@pytest.mark.correctness
def test_unrecognized_arg():
    result = subprocess.run(
        [sys.executable, "-m", "rocco.rocco", "--unrecognized_arg"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert result.returncode != 0
    assert "unrecognized" in result.stderr.decode()
