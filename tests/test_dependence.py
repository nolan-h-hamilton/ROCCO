import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rocco.dependence import (
    choose_dependence_span,
    estimate_dependence_radius_for_window,
)

STEP_BP = 50
WINDOW_BINS = 1000


def _gaussian_window(seed, sigma_bins=4.0):
    noise = np.random.default_rng(seed).normal(size=WINDOW_BINS)
    return gaussian_filter1d(noise, sigma=sigma_bins, mode="wrap")


@pytest.mark.correctness
def test_finite_pair_estimator_synthetic_accuracy_and_invariance():
    sigma_bins = 4.0
    complete_radii = []
    missing_radii = []
    for seed in range(100, 131):
        values = _gaussian_window(seed, sigma_bins)
        complete = estimate_dependence_radius_for_window(values[None, :], STEP_BP)
        missing_values = values.copy()
        missing_values[np.random.default_rng(seed + 1000).random(WINDOW_BINS) < 0.2] = (
            np.nan
        )
        missing = estimate_dependence_radius_for_window(
            missing_values[None, :], STEP_BP
        )
        duplicated = estimate_dependence_radius_for_window(
            np.vstack((values, values)), STEP_BP
        )
        complete_radii.append(complete["gaussianEquivalentRadiusBP"])
        missing_radii.append(missing["gaussianEquivalentRadiusBP"])
        duplication_drift = (
            abs(
                duplicated["gaussianEquivalentRadiusBP"]
                - complete["gaussianEquivalentRadiusBP"]
            )
            / complete["gaussianEquivalentRadiusBP"]
        )
        assert duplication_drift < 0.01

    expected_radius_bp = 3.0 * sigma_bins * STEP_BP
    complete_median = float(np.median(complete_radii))
    missing_median = float(np.median(missing_radii))
    assert abs(complete_median - expected_radius_bp) / expected_radius_bp <= 0.10
    assert abs(missing_median - complete_median) / complete_median <= 0.15

    AR1_medians = []
    for phi in (0.9, -0.9):
        crossings = []
        for seed in range(100, 131):
            rng = np.random.default_rng(seed)
            innovations = rng.normal(size=WINDOW_BINS)
            values = np.empty(WINDOW_BINS, dtype=np.float64)
            values[0] = innovations[0]
            for bin_index in range(1, WINDOW_BINS):
                values[bin_index] = phi * values[bin_index - 1] + innovations[bin_index]
            result = estimate_dependence_radius_for_window(values[None, :], STEP_BP)
            crossings.append(result["rawCrossingLagBP"])
        AR1_medians.append(float(np.median(crossings)))
    AR1_difference = abs(AR1_medians[0] - AR1_medians[1])
    assert AR1_difference <= max(STEP_BP, 0.10 * max(AR1_medians))


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("signal_kind", "right_censored"),
    [("damped", False), ("persistent", True)],
)
def test_finite_pair_estimator_reports_oscillation_and_censoring(
    signal_kind, right_censored
):
    bins = np.arange(WINDOW_BINS, dtype=np.float64)
    if signal_kind == "persistent":
        values = np.cos(2.0 * np.pi * bins * STEP_BP / 1000.0)
    else:
        kernel_bp = np.arange(400, dtype=np.float64) * STEP_BP
        kernel = np.exp(-kernel_bp / 500.0) * np.cos(2.0 * np.pi * kernel_bp / 1000.0)
        values = np.convolve(
            np.random.default_rng(52).normal(size=WINDOW_BINS),
            kernel,
            mode="full",
        )[:WINDOW_BINS]

    result = estimate_dependence_radius_for_window(values[None, :], STEP_BP)

    assert result["rightCensored"] is right_censored
    assert result["dominantACFPeriodBP"] is not None
    assert 0.0 < result["oscillationStrength"] <= 1.0
    assert result["postCrossingACFRevival"] > 0.0
    assert result["finitePairCountMin"] > 0.0
    if right_censored:
        assert result["rawCrossingLagBP"] is None
        assert result["lagsBP"][-1] == 49_950


@pytest.mark.correctness
def test_choose_dependence_span_gaussian_radius_and_genome_contract():
    matrices = {}
    coordinates = {}
    seed = 400
    for chromosome_index in range(1, 5):
        track_windows = [[] for _ in range(8)]
        for _ in range(5):
            for windows in track_windows:
                windows.append(_gaussian_window(seed))
                seed += 1
        chromosome = f"chr{chromosome_index}"
        matrices[chromosome] = np.vstack(
            [np.concatenate(windows) for windows in track_windows]
        )
        coordinates[chromosome] = STEP_BP * np.arange(5 * WINDOW_BINS)

    estimate, lower, upper, diagnostics = choose_dependence_span(
        matrices,
        coordinates,
        STEP_BP,
        window_count=20,
        bootstrap_draws=40,
        random_seed=34,
    )

    expected_unconstrained_radius_bp = 3.0 * 4.0 * STEP_BP
    assert (
        abs(
            diagnostics["unconstrainedEstimateBP"]
            - expected_unconstrained_radius_bp
        )
        / expected_unconstrained_radius_bp
        <= 0.10
    )
    assert diagnostics["estimateBP"] == pytest.approx(2500.0)
    assert diagnostics["lowerBP"] >= 2500.0
    assert diagnostics["upperBP"] >= 2500.0
    assert diagnostics["workingSpanBP"] >= 2500.0
    assert diagnostics["minimumCorrelationRadiusApplied"] is True
    assert estimate == math.ceil(diagnostics["estimateBP"] / STEP_BP)
    assert lower == math.ceil(diagnostics["lowerBP"] / STEP_BP)
    assert upper == math.ceil(diagnostics["upperBP"] / STEP_BP)
    assert diagnostics["workingSpanIntervals"] == math.ceil(
        diagnostics["workingSpanBP"] / STEP_BP
    )
    assert diagnostics["workingQuantile"] == pytest.approx(0.95)
    assert diagnostics["windowCountSelected"] == 20
    assert diagnostics["selectedAutosomeCount"] == 4
    assert diagnostics["censorFraction"] == pytest.approx(0.0)


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("values", "kwargs", "message"),
    [
        (np.arange(10.0), {}, "two-dimensional"),
        (np.empty((1, 0)), {}, "must not be empty"),
        (np.ones((1, 20)), {"step_bp": 0}, "positive"),
        (np.ones((1, 20)), {"acf_threshold": 1.0}, "lie in"),
        (np.ones((1, 20)), {}, "nonconstant finite"),
    ],
)
def test_window_estimator_rejects_invalid_inputs(values, kwargs, message):
    step_bp = kwargs.pop("step_bp", STEP_BP)
    with pytest.raises((TypeError, ValueError), match=message):
        estimate_dependence_radius_for_window(values, step_bp, **kwargs)


@pytest.mark.correctness
def test_choose_dependence_span_uses_prior_for_insufficient_support():
    matrices = {"chr1": _gaussian_window(9)[None, :]}
    coordinates = {"chr1": STEP_BP * np.arange(WINDOW_BINS)}

    estimate, lower, upper, diagnostics = choose_dependence_span(
        matrices,
        coordinates,
        STEP_BP,
        window_count=20,
        bootstrap_draws=20,
    )
    assert (estimate, lower, upper) == (100, 100, 100)
    assert diagnostics["usedPrior"] is True
    assert diagnostics["workingSpanIntervals"] == 100


@pytest.mark.correctness
@pytest.mark.parametrize(
    ("matrices", "coordinates", "kwargs", "message"),
    [
        ({}, {}, {}, "must not be empty"),
        (
            {"chr1": np.arange(WINDOW_BINS)},
            {"chr1": STEP_BP * np.arange(WINDOW_BINS)},
            {},
            "two-dimensional",
        ),
        (
            {"chr1": np.ones((1, WINDOW_BINS))},
            {},
            {},
            "keys must match",
        ),
        (
            {"chr1": np.ones((1, WINDOW_BINS))},
            {"chr1": STEP_BP * np.arange(WINDOW_BINS)},
            {"window_count": 19},
            "at least 20",
        ),
        (
            {"chr1": np.ones((1, WINDOW_BINS))},
            {"chr1": STEP_BP * np.arange(WINDOW_BINS)},
            {"bootstrap_draws": 19},
            "at least 20",
        ),
    ],
)
def test_choose_dependence_span_rejects_malformed_inputs(
    matrices, coordinates, kwargs, message
):
    with pytest.raises((TypeError, ValueError), match=message):
        choose_dependence_span(
            matrices,
            coordinates,
            STEP_BP,
            **kwargs,
        )


@pytest.mark.correctness
def test_persistent_oscillation_terminates_in_subprocess():
    script = "\n".join(
        [
            "import numpy as np",
            "from rocco.dependence import estimate_dependence_radius_for_window",
            "x = np.cos(2 * np.pi * np.arange(1000) / 20)",
            "result = estimate_dependence_radius_for_window(x[None, :], 50)",
            "assert result['rightCensored']",
        ]
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=5,
    )
    assert completed.returncode == 0, completed.stderr
