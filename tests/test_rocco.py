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
from scipy import signal

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


def _stationaryBlockMeta(numLoci, blockLength, kwargs):
    return dict(budgetFraction=0.05, effectiveCount=0.05 * numLoci,
        effectiveTotalCount=float(numLoci), numLoci=numLoci,
        bootstrapMethod="stationary_bootstrap", bootstrapBlockLength=blockLength,
        numBootstrap=kwargs["numBootstrap"], randomSeed=kwargs["randomSeed"],
        thresholdZ=kwargs["thresholdZ"], useLocalBootstrapRadius=kwargs["useLocalBootstrapRadius"],
        localRadiusIntervals=min(numLoci - 1, blockLength))


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
        np.array([[1.0, 3.0, 15.0]]),
        interval_bp=50,
        lower_bound_z=0.0,
        return_details=True,
    )
    assert details["input_scale"] == "log2p1"
    assert "sample_intercepts" not in details
    assert "sample_baselines" not in details
    assert np.allclose(details["mean"], np.array([-4.0 / 3.0, -1.0 / 3.0, 5.0 / 3.0]))
    assert details["centeringWindowBP"] == 1_250_000
    assert details["centeringWindowBins"] == 3
    expected = (details["mean"] + 1.0) / (details["standard_error"] + 1.0)
    assert np.allclose(details["z_scores"], expected)
    assert np.allclose(scores, expected + 1.0)


@pytest.mark.correctness
def test_score_loci_wls_explicit_min_effect_shrinks_standardized_score():
    scores, details = score_loci_wls(
        np.array([[1.0, 3.0, 15.0]]),
        interval_bp=50,
        min_effect=0.5,
        return_details=True,
    )
    assert np.isclose(details["min_effect"], 0.5)
    expected = (details["mean"] - 0.5 + 1.0) / (details["standard_error"] + 1.0)
    assert np.allclose(scores, expected + 1.0)


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
def test_score_loci_wls_returns_centered_matrix_for_dependence_estimation():
    scores, details = score_loci_wls(
        np.array([[1.0, 3.0, 7.0], [1.2, 2.8, 6.5]]),
        interval_bp=50,
        return_details=True,
    )
    assert scores.dtype == np.float64
    assert details["centered_matrix"].dtype == np.float64
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
        "chr1": {"summitTrackFile": str(summit_track_file)},
        "chr2": {"summitTrackFile": None},
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
    centered = np.zeros((3, 250000), dtype=np.float64)
    scores, details = ROCCO_INFERENCE._score_centered_wls_matrix(
        centered,
        lower_bound_z=1.0,
        prior_df=5.0,
    )
    assert scores.shape == (250000,)
    assert np.allclose(details["mean"], 0.0)
    expected = (details["mean"] + 1.0) / (details["standard_error"] + 1.0)
    assert np.allclose(details["z_scores"], expected)
    assert np.allclose(scores, expected)
    assert np.all(details["standard_error"] > 0.0)


@pytest.mark.correctness
def test_native_wls_downweights_noisy_track_locally():
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
@pytest.mark.parametrize(("n_loci", "interval_bp"), [(31, 50), (32, 100), (17, 250000)])
def test_order_zero_savgol_prefix_sum_parity(n_loci, interval_bp):
    from scipy.signal import savgol_filter

    x = np.arange(n_loci, dtype=np.float64)
    values = np.vstack((np.sin(x / 3.0), (x - 7.0) ** 2))
    window_bins = ROCCO_INFERENCE._resolve_centering_window_bins(n_loci, interval_bp)
    baseline = ROCCO_INFERENCE._savgol_order_zero_baseline(values, window_bins)
    expected = savgol_filter(
        values,
        window_length=window_bins,
        polyorder=0,
        axis=1,
        mode="interp",
    )

    assert np.allclose(baseline, expected, rtol=1.0e-12, atol=1.0e-12)
    assert window_bins % 2 == 1
    assert window_bins <= n_loci

    assert ROCCO_INFERENCE._resolve_centering_window_bins(30000, 50) == 25001
    assert ROCCO_INFERENCE._resolve_centering_window_bins(20000, 100) == 12501


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
    assert meta["priorStrength"] > 0
    assert meta["priorDispersion"] >= meta["minimumPriorDispersion"]
    assert meta["posteriorSummary"] == "beta_quantile"
    assert np.isclose(meta["posteriorQuantile"], 0.01)
    assert budgets["chr1"] < raw["chr1"]
    assert budgets["chr2"] < raw["chr2"]
    assert budgets["chr1"] < budgets["chr3"] < budgets["chr2"]


@pytest.mark.correctness
def test_stationary_bootstrap_native_contract():
    native_draw = ROCCO_INFERENCE._wls.stationaryNullBootstrapDraw
    template = np.arange(20, dtype=np.float64)
    expected_sources = {
        7: [7, 5, 1, 2, 3, 4, 5, 6, 7, 9, 6, 13, 14, 15, 10, 11, 12, 13, 19, 0],
        2: [2, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 14, 13, 16, 17, 18, 19, 18, 19],
        -1: [11, 4, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 17, 4, 5, 6, 7, 18, 19],
    }
    restart_indices = np.array([0, 1, 2, 9, 10, 11, 13, 14, 18])
    for radius, source in expected_sources.items():
        source = np.asarray(source, dtype=np.float64)
        rng = np.random.Generator(np.random.PCG64(0))
        draw = native_draw(template, 2, rng, radius)
        np.testing.assert_array_equal(draw, source - np.mean(source))
        assert draw.dtype == np.float64 and draw.flags.c_contiguous
        assert np.mean(draw) == pytest.approx(0.0, abs=1.0e-15)
        assert rng.bit_generator.random_raw() == 12646017539498340653
        if radius >= 0:
            assert np.all(np.abs(source[restart_indices] - restart_indices) <= radius)
    source = expected_sources[7]
    np.testing.assert_array_equal(np.diff(source[2:9]), np.ones(6))
    np.testing.assert_array_equal(source[18:20], np.array([19.0, 0.0]))

    invalid_calls = ((np.empty(0), 2, np.random.default_rng(0), -1),
        (np.zeros((2, 2)), 2, np.random.default_rng(0), -1),
        (np.array([0.0, np.nan]), 2, np.random.default_rng(0), -1),
        (template, 0, np.random.default_rng(0), -1),
        (template, 2, np.random.default_rng(0), -2), (template, 2, object(), -1))
    for call in invalid_calls:
        with pytest.raises((TypeError, ValueError)):
            native_draw(*call)


@pytest.mark.correctness
def test_stationary_budget_estimator_contract(monkeypatch):
    negative = -np.linspace(0.1, 10.0, 1_000)
    positive = np.linspace(0.05, 25.0, 1_000)
    positive[::100] = 5.0
    scores = np.empty(2_000)
    scores[0::2], scores[1::2] = negative, positive
    template, meta = ROCCO_INFERENCE._prepareStationaryNullTemplate(scores, 0.0)
    positive_mask = scores > 0.0
    ranks = ROCCO_INFERENCE.stats.rankdata(scores[positive_mask], method="average")
    reflected = scores.copy()
    reflected[positive_mask] = np.quantile(
        -scores[scores < 0.0], (ranks - 0.5) / ranks.size,
        method="interpolated_inverted_cdf")
    caps = np.quantile(reflected, (0.001, 0.999), method="interpolated_inverted_cdf")
    expected = np.clip(reflected, *caps)
    expected -= np.mean(expected)
    lower_scale = np.median(-scores[scores < 0.0]) / ROCCO_INFERENCE.stats.norm.ppf(0.75)
    expected *= lower_scale / np.sqrt(np.mean(expected**2))
    np.testing.assert_allclose(template, expected)
    assert set(meta) == {"templateLowerSize", "templateClipLower", "templateClipUpper"}
    assert (meta["templateClipLower"], meta["templateClipUpper"]) == pytest.approx(caps)
    counts = np.arange(1, 9)
    draw_index = [0]
    radius_calls = []
    def fake_draw(template, block_length, rng, radius):
        radius_calls.append(radius)
        draw = np.full(8, -2.0)
        draw[-counts[draw_index[0]] :] = 2.0
        draw_index[0] += 1
        return draw
    with monkeypatch.context() as patch:
        patch.setattr(ROCCO_INFERENCE, "_estimateStationaryNull", lambda _: (0.0, 1.0))
        patch.setattr(
            ROCCO_INFERENCE,
            "_prepareStationaryNullTemplate",
            lambda *_: (np.array([-2, -1, 0, 0, 0, 0, 1, 2], dtype=float), {}),
        )
        patch.setattr(ROCCO_INFERENCE._wls, "stationaryNullBootstrapDraw", fake_draw)
        budget, strict_meta = estimateStationaryBootstrapBudget(
            np.array([-5, -4, -3, -2, 0, 0, 0, 1], dtype=float),
            3,
            thresholdZ=0.0,
            numBootstrap=8,
            returnDetails=True,
        )
    null_occupancies = counts / 8.0
    null_sd = np.std(null_occupancies, ddof=1)
    assert strict_meta["trackTailOccupancy"] == 1.0 / 8.0
    assert strict_meta["nullTailOccupancy"] == pytest.approx(np.mean(null_occupancies))
    assert strict_meta["nullTailOccupancySD"] == pytest.approx(null_sd)
    assert strict_meta["nullTailMCSE"] == pytest.approx(null_sd / np.sqrt(8.0))
    assert budget == strict_meta["budgetFraction"] == max(strict_meta["signedTailExcess"], 0.0)
    assert strict_meta["localRadiusIntervals"] == 5 and radius_calls == [20_000] * 8
    rng = np.random.default_rng(20260831)
    size, rho = 8_192, 0.75
    null_tracks = (
        rng.normal(size=size),
        signal.lfilter(
            [np.sqrt(1.0 - rho**2)], [1.0, -rho], rng.normal(size=size)
        ),
    )
    for null_track in null_tracks:
        medians = []
        for dose in (0.0, 0.01, 0.04):
            planted = null_track.copy()
            planted[: int(dose * size)] += 7.0
            values = [
                estimateStationaryBootstrapBudget(planted, 32, numBootstrap=64,
                    randomSeed=seed, returnDetails=True)[1]["signedTailExcess"]
                for seed in (17, 31, 59)
            ]
            replay = estimateStationaryBootstrapBudget(planted, 32, numBootstrap=64,
                randomSeed=17, returnDetails=True)[1]["signedTailExcess"]
            assert replay == values[0]
            medians.append(float(np.median(values)))
        assert abs(medians[0]) <= 0.005
        assert medians[0] <= medians[1] <= medians[2]
        assert medians[2] - medians[0] >= 0.025

    invalid_calls = (([], 2), (np.zeros((2, 2)), 2), ([0.0, np.inf], 2), (scores, 0))
    for call in invalid_calls:
        with pytest.raises(ValueError):
            estimateStationaryBootstrapBudget(*call)


@pytest.mark.correctness
def test_empirical_bayes_budget_single_chrom_uses_default_center():
    budgets, meta = estimate_empirical_bayes_budgets(
        {"chr1": 0},
        {"chr1": 0},
    )
    assert np.isclose(meta["genomeWideBudgetFraction"], 0.05)
    assert meta["posteriorSummary"] == "beta_quantile"
    assert np.isclose(meta["posteriorQuantile"], 0.01)
    assert 0.0 < budgets["chr1"] < meta["genomeWideBudgetFraction"]


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
        chrom: {
            "budgetMetadata": {
                "budgetBlocks": [{
                    "blockID": 0,
                    "numLoci": 1000,
                    "budgetFraction": count / 1000.0,
                    "effectiveCount": count,
                    "effectiveTotalCount": 1000.0,
                }]
            }
        }
        for chrom, count in {"chr1": 0.0, "chr2": 2.0, "chr3": 120.0,
                             "chr4": 400.0}.items()
    }
    budgets, meta = ROCCO_IMPL._resolve_budgets(
        chrom_cache,
        {
            "budget_posterior_quantile": 0.01,
            "budget": None,
            "scale_chrom_budgets": 1.0,
        },
    )
    assert meta["posteriorSummary"] == "beta_quantile"
    for chrom in chrom_cache:
        assert 0.001 <= budgets[chrom]["posteriorBudgetFraction"] <= 0.25


@pytest.mark.correctness
def test_resolve_budgets_pools_chrom_block_units(monkeypatch):
    def blockMetadata(blockID, startIndex, stopIndex, effectiveCount):
        return {
            "blockID": blockID,
            "startIndex": startIndex,
            "stopIndex": stopIndex,
            "startBP": 50 * startIndex,
            "stopBP": 50 * stopIndex,
            "numLoci": stopIndex - startIndex,
            "budgetFraction": effectiveCount / 10.0,
            "effectiveCount": effectiveCount,
            "effectiveTotalCount": 10.0,
        }

    chrom_cache = {
        "chr1": {
            "budgetMetadata": {
                "budgetBlocks": [
                    blockMetadata(0, 0, 5, 1.0),
                    blockMetadata(1, 5, 10, 2.0),
                ]
            },
        },
        "chr2": {
            "budgetMetadata": {
                "budgetBlocks": [
                    blockMetadata(0, 0, 4, 3.0),
                    blockMetadata(1, 4, 8, 4.0),
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
            {"genomeWideBudgetFraction": 0.1, "posteriorSummary": "fake"},
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
    assert budgets["chr1"]["posteriorBudgetFraction"] == pytest.approx(0.015)
    assert [block["posteriorBudgetFraction"]
            for block in budgets["chr2"]["budgetBlocks"]] == [0.03, 0.04]
    assert np.allclose(
        ROCCO_IMPL._solver_budget_blocks(budgets["chr1"]),
        np.array([[0.0, 5.0, 0.01], [5.0, 10.0, 0.02]]),
    )
    assert budgets["chr2"]["budgetBlocks"][1]["budgetFraction"] == pytest.approx(0.4)
    assert meta["pooledBudgetBlockCount"] == 4


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
        return (
            np.array([1, 0, 1, 0], dtype=np.uint8),
            12.0,
            {
                "budget_mode": "block_soft_selection_penalty",
                "soft_budget_penalty": 1.5,
                "selected_count": 2,
            },
        )

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
                "posteriorBudgetFraction": 0.15,
                "budgetBlocks": [
                    {"startIndex": 0, "stopIndex": 2,
                     "posteriorBudgetFraction": 0.1},
                    {"startIndex": 2, "stopIndex": 4,
                     "posteriorBudgetFraction": 0.2},
                ],
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
def test_build_chrom_cache_bam_uses_one_score_estimator_path(monkeypatch):
    chrom_lengths = {"chr1": 8, "chr2": 12}
    score_calls = []
    dependence_calls = []
    budget_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        n = chrom_lengths[chrom]
        matrix = np.vstack((np.arange(n), np.arange(n) + 10.0))
        return 50 * np.arange(n, dtype=int), matrix

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        score_calls.append(kwargs.copy())
        n = chrom_matrix.shape[1]
        scores = np.linspace(0.0, 3.0, n, dtype=float)
        centered = np.asarray(chrom_matrix, dtype=float) + 100.0
        return scores, {
            "centered_matrix": centered,
            "mean": np.linspace(10.0, 13.0, n, dtype=float),
            "z_scores": scores + 1.0,
            "centeringMethod": "savgolOrder0",
            "centeringWindowBP": 1_250_000,
            "centeringWindowBins": n if n % 2 else n - 1,
        }

    def fake_choose_dependence_span(matrices, coordinates, step_bp, **kwargs):
        dependence_calls.append(
            {
                "matrices": {
                    chromosome: np.asarray(matrix).copy()
                    for chromosome, matrix in matrices.items()
                },
                "coordinates": {
                    chromosome: np.asarray(values).copy()
                    for chromosome, values in coordinates.items()
                },
                "step_bp": step_bp,
                "kwargs": kwargs,
            }
        )
        return (
            3,
            2,
            4,
            {
                "estimateBP": 150.0,
                "lowerBP": 100.0,
                "upperBP": 200.0,
                "workingSpanBP": 550.0,
                "workingSpanIntervals": 11,
            },
        )

    def fake_budget_estimator(score_track, block_length, **kwargs):
        scores = np.asarray(score_track, dtype=float)
        observed_len = int(scores.size)
        budget_calls.append((scores.copy(), block_length, kwargs.copy()))
        return 0.05, _stationaryBlockMeta(observed_len, block_length, kwargs)

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fake_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL, "choose_dependence_span", fake_choose_dependence_span
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimateStationaryBootstrapBudget",
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
        "score_min_effect": None,
        "score_precision_floor_ratio": 0.01,
        "budget_bootstrap_draws": 8,
        "budget_threshold_z": 2.0,
        "budget_bootstrap_seed": 42,
        "no_local_bootstrap_radius": False,
        "num_null_blocks": 4,
        "gamma": 2.5,
        "peak_mode": None,
        "window_bp": 50000,
        "window_count": 256,
        "working_quantile": 0.9,
        "bootstrap_draws": 500,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr1", "chr2"],
        ["track1.bam", "track2.bam"],
        args,
    )

    assert len(dependence_calls) == 1
    assert dependence_calls[0]["step_bp"] == 50
    assert [call["interval_bp"] for call in score_calls] == [50, 50]
    for chromosome, n_loci in chrom_lengths.items():
        expected_centered = (
            np.vstack((np.arange(n_loci), np.arange(n_loci) + 10.0)) + 100.0
        )
        assert np.array_equal(
            dependence_calls[0]["matrices"][chromosome], expected_centered
        )
        chromosome_budget_calls = budget_calls[:4]
        budget_calls = budget_calls[4:]
        assert np.allclose(
            np.concatenate([call[0] for call in chromosome_budget_calls]),
            np.linspace(0.0, 3.0, n_loci),
        )
        assert [call[1] for call in chromosome_budget_calls] == [11] * 4
        assert [call[2]["numBootstrap"] for call in chromosome_budget_calls] == [8] * 4
        diagnostics = chrom_cache[chromosome]["dependenceDiagnostics"]
        assert diagnostics["estimateBP"] == 150.0
        assert diagnostics["workingSpanIntervals"] == 11
        assert diagnostics["workingSpanBP"] == 550.0
        assert chrom_cache[chromosome]["gamma"] == 2.5
    all_budget_blocks = [
        block
        for chromosome in chrom_cache.values()
        for block in chromosome["budgetMetadata"]["budgetBlocks"]
    ]
    assert [block["bootstrapBlockLength"] for block in all_budget_blocks] == [11] * 8
    assert all({"blockID", "startIndex", "stopIndex", "startBP", "stopBP",
                "numLoci", "budgetFraction", "effectiveCount",
                "effectiveTotalCount"} <= set(block) for block in all_budget_blocks)
    assert all(not any("_" in key for key in block) for block in all_budget_blocks)


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("peak_mode", "expected_summit_calls"),
    [(None, 0), ("narrow", 1), ("broad", 0), ("both", 1)],
)
def test_build_chrom_cache_peak_modes(monkeypatch, peak_mode, expected_summit_calls):
    summit_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        return 50 * np.arange(64, dtype=int), np.zeros((1, 64), dtype=float)

    def fake_score_loci_wls(chrom_matrix, **kwargs):
        scores = np.linspace(-1.0, 2.5, chrom_matrix.shape[1], dtype=float)
        return scores, {
            "centered_matrix": np.zeros_like(chrom_matrix, dtype=float),
            "mean": scores.copy(),
            "z_scores": scores + 2.0,
            "centeringMethod": "savgolOrder0",
            "centeringWindowBP": 1_250_000,
            "centeringWindowBins": 63,
        }

    def fake_choose_dependence_span(*args, **kwargs):
        return (
            5,
            4,
            6,
            {
                "estimateBP": 250.0,
                "lowerBP": 200.0,
                "upperBP": 300.0,
                "workingSpanBP": 250.0,
                "workingSpanIntervals": 5,
            },
        )

    def fake_budget_estimator(score_track, block_length, **kwargs):
        observed_len = int(np.asarray(score_track).size)
        return 0.05, _stationaryBlockMeta(observed_len, block_length, kwargs)

    def fake_summit_track(chrom, intervals, effect_mean):
        summit_calls.append((chrom, intervals.copy(), effect_mean.copy()))
        return None

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fake_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_IMPL, "choose_dependence_span", fake_choose_dependence_span
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimateStationaryBootstrapBudget",
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
        "score_min_effect": None,
        "score_precision_floor_ratio": 0.01,
        "budget_bootstrap_draws": 8,
        "budget_threshold_z": 2.0,
        "budget_bootstrap_seed": 42,
        "no_local_bootstrap_radius": False,
        "num_null_blocks": 4,
        "gamma": 0.25,
        "peak_mode": peak_mode,
        "window_bp": 50000,
        "window_count": 256,
        "working_quantile": 0.9,
        "bootstrap_draws": 500,
    }

    chrom_cache = ROCCO_IMPL._build_chrom_cache(
        ["chr1"],
        ["track.bam"],
        args,
    )

    assert len(summit_calls) == expected_summit_calls
    assert "effectMean" not in chrom_cache["chr1"]
    assert chrom_cache["chr1"]["gamma"] == pytest.approx(0.25)
    assert chrom_cache["chr1"]["dependenceDiagnostics"]["estimateBP"] == 250.0


@pytest.mark.correctness
def test_broad_parent_merging_uses_working_span_not_median_radius(monkeypatch):
    intervals = 50 * np.arange(8, dtype=np.int64)
    weak_solution = np.array([1, 1, 0, 0, 0, 0, 1, 1], dtype=np.uint8)

    def fake_solve_chrom_exact(*args, **kwargs):
        return (
            weak_solution,
            0.0,
            {
                "budget_mode": "unpenalized",
                "soft_budget_penalty": 0.0,
                "selected_count": 4,
            },
        )

    monkeypatch.setattr(ROCCO_IMPL, "solve_chrom_exact", fake_solve_chrom_exact)
    chrom_cache = {
        "chr1": {
            "solution": np.array([1, 0, 0, 0, 0, 0, 1, 0], dtype=np.uint8),
            "intervals": intervals,
            "zScores": np.ones(8),
            "gamma": 1.0,
            "intervalBP": 50,
            "dependenceDiagnostics": {"workingSpanIntervals": 3},
        }
    }
    records, block_map = ROCCO_IMPL._build_broad_parent_records(
        chrom_cache,
        {"chr1": {"posteriorBudgetFraction": 0.05}},
        {
            "broad_score_lower_bound_z": 0.0,
            "broad_max_gap_bp": None,
            "broad_min_peak_bp": 1,
        },
    )

    assert records == [("chr1", 0, 400)]
    assert block_map[("chr1", 0, 400)] == [(0, 50), (300, 350)]


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
def test_generate_chrom_matrix_smooths_each_bam_by_fragment_length(
    monkeypatch, tmp_path
):
    bam_a = tmp_path / "frag75.bam"
    bam_b = tmp_path / "frag10.bam"
    chrom_sizes_path = tmp_path / "toy_frag_smooth.sizes"
    bam_a.write_bytes(b"")
    bam_b.write_bytes(b"")
    chrom_sizes_path.write_text("chr1\t125\n", encoding="utf-8")
    monkeypatch.setattr(ROCCO_READTRACKS, "_BAM_COUNT_METADATA_CACHE", {})

    class FakeNativeCounts:
        def is_alignment_paired_end(self, *args, **kwargs):
            return False

        def get_alignment_read_length(self, *args, **kwargs):
            return 25

        def get_alignment_mapped_read_count(self, *args, **kwargs):
            return 1000000, 1000000

        def get_alignment_fragment_length(self, bam_file, **kwargs):
            return 75 if bam_file == str(bam_a) else 10

        def get_alignment_chrom_range(self, *args, **kwargs):
            return 0, 125

        def count_alignment_region(self, bam_file, *args, **kwargs):
            if bam_file == str(bam_a):
                return np.array([0.0, 0.0, 4.0, 8.0, 0.0], dtype=np.float32)
            return np.array([0.0, 5.0, 0.0, 0.0, 0.0], dtype=np.float32)

    monkeypatch.setattr(ROCCO_READTRACKS, "_hts_counts", FakeNativeCounts())

    intervals, count_matrix = generate_chrom_matrix(
        "chr1",
        [str(bam_a), str(bam_b)],
        str(chrom_sizes_path),
        step=25,
        norm_method="CPM",
        flag_exclude=0,
        ignore_for_norm=[],
        round_digits=6,
        num_processors=1,
    )

    assert intervals.tolist() == [25, 50, 75, 100]
    assert np.allclose(
        count_matrix,
        np.array(
            [
                [1.333333, 4.0, 4.0, 4.0],
                [5.0, 0.0, 0.0, 0.0],
            ]
        ),
    )


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
    _write_toy_bigwig(
        bw_path,
        entries=[
            (0, 50, -7.1256789),
            (100, 150, 12.75),
            (150, 200, 0.25),
        ],
    )

    intervals, score_matrix = generate_chrom_matrix(
        "chr1",
        [str(bw_path)],
        str(chrom_sizes_path),
        step=999,
        round_digits=2,
    )

    assert intervals.tolist() == [0, 50, 100, 150]
    assert score_matrix.shape == (1, 4)
    assert score_matrix[0, 0] == pytest.approx(-7.1256789, abs=5.0e-7)
    assert np.isnan(score_matrix[0, 1])
    assert np.allclose(score_matrix[0, 2:], np.array([12.75, 0.25]))


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
def test_get_ecdf_counts_sampled_intervals_in_native_batches(monkeypatch, tmp_path):
    bam_a = tmp_path / "a.bam"
    bam_b = tmp_path / "b.bam"
    bam_a.write_bytes(b"")
    bam_b.write_bytes(b"")
    native_calls = []

    class FakeNativeCounts:
        def count_alignment_intervals(
            self,
            alignment_path,
            chromosomes,
            starts,
            ends,
            **kwargs,
        ):
            native_calls.append(
                (alignment_path, list(chromosomes), list(starts), list(ends), kwargs)
            )
            if alignment_path == str(bam_a):
                return np.array([2.0, 4.0], dtype=np.float32)
            return np.array([6.0, 8.0], dtype=np.float32)

    def fail_alignment_file(*args, **kwargs):
        raise AssertionError("ECDF sampling should use native interval batches")

    monkeypatch.setattr(
        ROCCO_SCORES,
        "_random_intervals",
        lambda *args, **kwargs: [("chr1", 10, 110), ("chr2", 20, 120)],
    )
    monkeypatch.setattr(ROCCO_SCORES, "_hts_counts", FakeNativeCounts())
    monkeypatch.setattr(ROCCO_SCORES.pysam, "AlignmentFile", fail_alignment_file)

    empirical_null = ROCCO_SCORES.get_ecdf(
        [str(bam_a), str(bam_b)],
        length=100,
        chrom_sizes_file=str(tmp_path / "chrom.sizes"),
        nsamples=2,
        sample_scaling_constants=[0.5, 2.0],
        seed=7,
        row_scale=100,
        pc=1,
        thread_count=3,
    )

    transformed = np.log2(np.array([[2.0, 13.0], [3.0, 17.0]]))
    expected = np.sort(np.percentile(transformed, 75, axis=1))
    assert np.allclose(empirical_null.values, expected)
    assert len(native_calls) == 2
    for _, chromosomes, starts, ends, kwargs in native_calls:
        assert chromosomes == ["chr1", "chr2"]
        assert starts == [10, 20]
        assert ends == [110, 120]
        assert kwargs["one_read_per_bin"] == 1
        assert kwargs["thread_count"] == 3
        assert (kwargs["min_mapping_quality"], kwargs["flag_exclude"]) == (20, 3844)
        assert kwargs["count_mode"] == "coverage"


@pytest.mark.correctness
def test_multi_ecdf_serial_path_avoids_pool(monkeypatch, tmp_path):
    bam_path = tmp_path / "toy.bam"
    bam_path.write_bytes(b"")
    calls = []

    def fake_get_ecdf(
        bam_files,
        length,
        chrom_sizes_file,
        nsamples,
        sample_scaling_constants,
        seed,
        null_stat,
        trim_proportion,
        row_scale,
        pc,
        thread_count,
    ):
        calls.append(
            (
                bam_files,
                int(length),
                chrom_sizes_file,
                int(nsamples),
                sample_scaling_constants,
                seed,
                row_scale,
                pc,
                thread_count,
            )
        )
        return ROCCO_SCORES.EmpiricalNull(np.array([float(length)]))

    def fail_get_context(*args, **kwargs):
        raise AssertionError("proc=1 should not create a multiprocessing pool")

    monkeypatch.setattr(ROCCO_SCORES, "get_ecdf", fake_get_ecdf)
    monkeypatch.setattr(ROCCO_SCORES.multiprocessing, "get_context", fail_get_context)

    result = ROCCO_SCORES.multi_ecdf(
        [str(bam_path)],
        np.array([200, 100, 200]),
        str(tmp_path / "chrom.sizes"),
        nsamples_per_length=9,
        sample_scaling_constants=[1.5],
        seed=11,
        proc=1,
        row_scale=50,
        pc=2,
        thread_count=4,
    )

    assert list(result.keys()) == [100, 200]
    assert [call[1] for call in calls] == [100, 200]
    assert all(call[3:] == (9, [1.5], 11, 50, 2, 4) for call in calls)


@pytest.mark.correctness
def test_score_peaks_regenerates_stale_count_matrix(monkeypatch, tmp_path):
    bam_path = tmp_path / "toy.bam"
    peak_path = tmp_path / "peaks.bed"
    count_path = tmp_path / "counts.tsv"
    output_path = tmp_path / "scored.bed"
    bam_path.write_bytes(b"")
    peak_path.write_text(
        "chr1\t0\t100\n" "chr1\t100\t200\n",
        encoding="utf-8",
    )
    count_path.write_text(
        "peak_name\ttoy\n" "chr1_0_100\t1\n",
        encoding="utf-8",
    )
    regenerate_calls = []

    def fake_raw_count_matrix(bam_files, peak_file, output_file, bed_columns=3):
        regenerate_calls.append((bam_files, peak_file, output_file, bed_columns))
        Path(output_file).write_text(
            "peak_name\ttoy\n" "chr1_0_100\t4\n" "chr1_100_200\t6\n",
            encoding="utf-8",
        )
        return output_file

    class FakeAlignmentFile:
        mapped = 100

        def __init__(self, *args, **kwargs):
            pass

        def count(self, *args, **kwargs):
            return 0

        def close(self):
            pass

    def fake_multi_ecdf(bam_files, lengths, chrom_sizes_file, **kwargs):
        return {
            int(length): ROCCO_SCORES.EmpiricalNull(np.array([0.0, 100.0]))
            for length in np.asarray(lengths)
        }

    monkeypatch.setattr(ROCCO_SCORES, "raw_count_matrix", fake_raw_count_matrix)
    monkeypatch.setattr(ROCCO_SCORES, "get_read_length", lambda *args, **kwargs: 100)
    monkeypatch.setattr(ROCCO_SCORES.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(ROCCO_SCORES, "multi_ecdf", fake_multi_ecdf)

    with ROCCO_SCORES.pd.option_context("mode.copy_on_write", True):
        scores, _, _ = ROCCO_SCORES.score_peaks(
            [str(bam_path)],
            chrom_sizes_file=str(tmp_path / "chrom.sizes"),
            peak_file=str(peak_path),
            count_matrix_file=str(count_path),
            effective_genome_size=10000,
            output_file=str(output_path),
            ecdf_nsamples=2,
            proc=1,
        )

    assert len(regenerate_calls) == 1
    assert scores.shape == (2,)
    output_rows = output_path.read_text(encoding="utf-8").splitlines()
    assert len(output_rows) == 2
    assert all(row.split("\t")[8] == "-1" for row in output_rows)


@pytest.mark.correctness
def test_build_chrom_cache_preserves_bigwig_rows_and_bypasses_wls(monkeypatch):
    raw_matrix = None
    dependence_calls = []
    direct_budget_calls = []

    def fake_generate_chrom_matrix(chrom, *args, **kwargs):
        return np.array([0, 50, 100, 150], dtype=int), raw_matrix.copy()

    def fail_score_loci_wls(*args, **kwargs):
        raise AssertionError("bigWig inputs should bypass WLS scoring")

    def fail_bam_centering(*args, **kwargs):
        raise AssertionError("bigWig inputs should bypass BAM centering")

    def fake_choose_dependence_span(matrices, coordinates, step_bp, **kwargs):
        dependence_calls.append(
            {
                "matrix": np.asarray(matrices["chr1"]).copy(),
                "step_bp": step_bp,
            }
        )
        return (
            3,
            2,
            4,
            {
                "estimateBP": 150.0,
                "lowerBP": 100.0,
                "upperBP": 200.0,
                "workingSpanBP": 550.0,
                "workingSpanIntervals": 11,
            },
        )

    def fake_budget_estimator(scores, block_length, **kwargs):
        scores_ = np.asarray(scores, dtype=float)
        direct_budget_calls.append((scores_, block_length, kwargs))
        return 0.05, _stationaryBlockMeta(scores_.size, block_length, kwargs)

    monkeypatch.setattr(ROCCO_IMPL, "generate_chrom_matrix", fake_generate_chrom_matrix)
    monkeypatch.setattr(ROCCO_IMPL, "score_loci_wls", fail_score_loci_wls)
    monkeypatch.setattr(
        ROCCO_INFERENCE, "_savgol_order_zero_baseline", fail_bam_centering
    )
    monkeypatch.setattr(
        ROCCO_IMPL,
        "estimateStationaryBootstrapBudget",
        fake_budget_estimator,
    )
    monkeypatch.setattr(
        ROCCO_IMPL, "choose_dependence_span", fake_choose_dependence_span
    )

    cases = [
        (
            np.array([[-2.125, 5.5, 12.75, 0.25]]),
            np.array([-2.125, 5.5, 12.75, 0.25]),
            "direct",
        ),
        (
            np.array([[-2.125, np.nan, 12.75, 0.25], [-1.25, 5.5, np.nan, 2.25]]),
            np.array([-1.6875, 5.5, 12.75, 1.25]),
            "directMedian",
        ),
    ]
    for case_index, (raw_values, expected_scores, scoring_method) in enumerate(cases):
        raw_matrix = raw_values
        input_files = [f"track{index}.bw" for index in range(raw_values.shape[0])]
        args = {
            "chrom_sizes_file": None,
            "step": 50,
            "round_digits": 2,
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
            "budget_bootstrap_draws": 8,
            "budget_threshold_z": 2.0,
            "budget_bootstrap_seed": 42,
            "no_local_bootstrap_radius": False,
            "num_null_blocks": 4,
            "gamma": 3.0,
            "peak_mode": None,
            "window_bp": 50000,
            "window_count": 256,
            "working_quantile": 0.9,
            "bootstrap_draws": 500,
        }
        budget_start = len(direct_budget_calls)
        chrom_cache = ROCCO_IMPL._build_chrom_cache(["chr1"], input_files, args)
        case_budget_calls = direct_budget_calls[budget_start:]

        assert np.array_equal(
            dependence_calls[case_index]["matrix"], raw_values, equal_nan=True
        )
        assert np.allclose(chrom_cache["chr1"]["scores"], expected_scores)
        assert chrom_cache["chr1"]["signalPreparation"] == "asProvided"
        assert chrom_cache["chr1"]["scoringMethod"] == scoring_method
        assert "centeringMethod" not in chrom_cache["chr1"]
        assert len(case_budget_calls) == 4
        assert np.allclose(
            np.concatenate([call[0] for call in case_budget_calls]), expected_scores
        )
        assert [call[1] for call in case_budget_calls] == [11] * 4
        assert [call[2]["numBootstrap"] for call in case_budget_calls] == [8] * 4
        diagnostics = chrom_cache["chr1"]["dependenceDiagnostics"]
        assert diagnostics["estimateBP"] == 150.0
        assert diagnostics["workingSpanIntervals"] == 11


@pytest.mark.correctness
def test_prepare_args_low_memory_preserves_bootstrap_draws(monkeypatch):
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
            "--no_local_bootstrap_radius",
        ],
    )
    args = ROCCO_IMPL._prepare_args(parser)
    assert args["low_memory"] is True
    assert 1 <= int(args["threads"]) <= 4
    assert int(args["budget_bootstrap_draws"]) == 64
    assert args["budget_threshold_z"] == pytest.approx(2.0)
    assert int(args["budget_bootstrap_seed"]) == 42
    assert args["no_local_bootstrap_radius"] is True
    assert int(args["num_null_blocks"]) == 4
    assert args["working_quantile"] == pytest.approx(0.95)


@pytest.mark.correctness
def test_ignore_for_norm_distinguishes_default_empty_and_values(tmp_path):
    bamFile = tmp_path / "fake.bam"
    bamFile.touch()
    parser = ROCCO_IMPL._build_parser()
    for extraArgs, expected in (
        ([], ["chrX", "chrY", "chrM"]),
        (["--ignore_for_norm"], []),
        (["--ignore_for_norm", "chrM"], ["chrM"]),
    ):
        args = vars(parser.parse_args(["-i", str(bamFile), *extraArgs]))
        args["input_track_type"] = "bam"
        ROCCO_IMPL._prepare_inputs(args)
        assert args["ignore_for_norm"] == expected


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
        (
            [
                "-i",
                "fake.bw",
                "--peak_mode",
                "broad",
            ],
            "BAM inputs",
        ),
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
        interval_bp=50,
        return_details=True,
    )
    assert scores.shape == (5000,)
    budget_fraction, budget_meta = estimateStationaryBootstrapBudget(
        scores,
        101,
        numBootstrap=8,
        returnDetails=True,
    )
    assert 0.0 <= budget_fraction <= 1.0
    assert budget_meta["effectiveCount"] >= 0.0
    assert 1.0 <= budget_meta["effectiveTotalCount"] <= budget_meta["numLoci"]
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
