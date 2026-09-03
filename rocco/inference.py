r"""
ROCCO: Inference
==================================================================================

Functions for scoring loci and setting chromosome-specific budgets.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Tuple

import numpy as np
from scipy import optimize, special, stats

from . import _wls

logger = logging.getLogger(__name__)
_CENTERING_WINDOW_BP = 1_250_000
_TINY = float(np.finfo(np.float64).tiny)


def _robust_scale(values: np.ndarray, floor: float = 1.0e-6) -> float:
    values_ = np.asarray(values, dtype=np.float64)
    if values_.size == 0:
        return float(floor)
    mad = np.median(np.abs(values_ - np.median(values_)))
    return float(max(mad * 1.4826, floor))


def _resolve_centering_window_bins(n_loci: int, interval_bp: int) -> int:
    n_loci_ = int(n_loci)
    interval_bp_ = int(interval_bp)
    if n_loci_ <= 0:
        raise ValueError("`chrom_matrix` must contain at least one locus")
    if interval_bp_ <= 0:
        raise ValueError("`interval_bp` must be positive")
    target_bins = int(np.floor((_CENTERING_WINDOW_BP / float(interval_bp_)) + 0.5))
    target_bins = max(1, target_bins)
    if (target_bins % 2) == 0:
        target_bins += 1
    if target_bins > n_loci_:
        target_bins = n_loci_ if (n_loci_ % 2) == 1 else n_loci_ - 1
    return int(max(1, target_bins))


def _savgol_order_zero_baseline(
    values: np.ndarray,
    window_bins: int,
) -> np.ndarray:
    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("`values` must be two-dimensional")
    n_loci = int(matrix.shape[1])
    window_bins_ = int(window_bins)
    if window_bins_ <= 0 or window_bins_ > n_loci or (window_bins_ % 2) == 0:
        raise ValueError("`window_bins` must be an odd value within the track length")

    prefix_sums = np.empty((matrix.shape[0], n_loci + 1), dtype=np.float64)
    prefix_sums[:, 0] = 0.0
    np.cumsum(matrix, axis=1, out=prefix_sums[:, 1:])
    rolling_means = (
        prefix_sums[:, window_bins_:] - prefix_sums[:, : n_loci - window_bins_ + 1]
    ) / float(window_bins_)
    half_window = window_bins_ // 2
    baseline = np.empty_like(matrix, dtype=np.float64)
    baseline[:, half_window : n_loci - half_window] = rolling_means
    if half_window > 0:
        baseline[:, :half_window] = rolling_means[:, :1]
        baseline[:, n_loci - half_window :] = rolling_means[:, -1:]
    return baseline


def _score_centered_wls_matrix(
    centered_matrix: np.ndarray,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    spatial_window: int | None = None,
    precision_floor_ratio: float = 0.01,
) -> tuple[np.ndarray, Dict[str, np.ndarray | float]]:
    centered = np.asarray(centered_matrix, dtype=np.float64)
    if centered.ndim != 2:
        raise ValueError("`centered_matrix` must be two-dimensional")
    if centered.shape[0] == 0 or centered.shape[1] == 0:
        raise ValueError("`centered_matrix` must be non-empty")
    precision_floor_ratio_ = float(max(precision_floor_ratio, 0.0))
    # run EB munc + WLS to get per-locus scores and details for the constrained optimization (DP)
    result_native = _wls.score_centered_wls(
        centered,
        lower_bound_z=float(lower_bound_z),
        prior_df=float(prior_df),
        min_effect=min_effect,
        spatial_window=31 if spatial_window is None else int(spatial_window),
        precision_floor_ratio=precision_floor_ratio_,
    )
    (
        scores_arr,
        mean_arr,
        raw_var_arr,
        prior_var_arr,
        moderated_var_arr,
        se_arr,
        total_df,
        resolved_window,
    ) = result_native
    se = np.asarray(se_arr, dtype=np.float64)
    mean = np.asarray(mean_arr, dtype=np.float64)
    shifted_z_scores = (mean + 1.0) / (se + 1.0)
    result = (
        np.asarray(scores_arr, dtype=np.float64),
        {
            "mean": mean,
            "raw_variance": np.asarray(raw_var_arr, dtype=np.float64),
            "prior_variance": np.asarray(prior_var_arr, dtype=np.float64),
            "moderated_variance": np.asarray(moderated_var_arr, dtype=np.float64),
            "standard_error": se,
            "z_scores": shifted_z_scores,
            "min_effect": float(0.0 if min_effect is None else max(min_effect, 0.0)),
            "precision_floor_ratio": float(precision_floor_ratio_),
            "degrees_of_freedom": np.full(
                centered.shape[1],
                float(total_df),
                dtype=np.float64,
            ),
            "prior_spatial_window": float(resolved_window),
        },
    )

    scores, details = result
    if (
        not np.all(np.isfinite(scores))
        or not np.all(np.isfinite(details["mean"]))
        or not np.all(np.isfinite(details["raw_variance"]))
        or not np.all(np.isfinite(details["prior_variance"]))
        or not np.all(np.isfinite(details["moderated_variance"]))
        or not np.all(np.isfinite(details["standard_error"]))
        or not np.all(np.isfinite(details["z_scores"]))
    ):
        raise ValueError("EB scoring produced non-finite values")
    return scores, details


def score_loci_wls(
    chrom_matrix: np.ndarray,
    interval_bp: int,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    return_details: bool = False,
) -> np.ndarray | Tuple[np.ndarray, Dict[str, Any]]:
    r"""Score loci with an EB-moderated summary on baseline-corrected log signal."""
    matrix = np.asarray(chrom_matrix, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("`chrom_matrix` must be two-dimensional")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("`chrom_matrix` must be non-empty")
    if np.any(~np.isfinite(matrix)):
        raise ValueError("`chrom_matrix` contains non-finite values")

    window_bins = _resolve_centering_window_bins(matrix.shape[1], interval_bp)
    matrix = np.log2(np.clip(matrix, 0.0, None) + 1.0)
    baselines = _savgol_order_zero_baseline(
        matrix,
        window_bins=window_bins,
    )
    centered = matrix - baselines
    del matrix
    del baselines
    scores, core_details = _score_centered_wls_matrix(
        centered,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
    )
    if not return_details:
        return scores.astype(np.float64)

    details = {
        "input_scale": "log2p1",
        "centeringMethod": "savgolOrder0",
        "centeringWindowBP": int(_CENTERING_WINDOW_BP),
        "centeringWindowBins": int(window_bins),
        "mean": np.asarray(core_details["mean"], dtype=np.float64),
        "raw_variance": np.asarray(core_details["raw_variance"], dtype=np.float64),
        "prior_variance": np.asarray(core_details["prior_variance"], dtype=np.float64),
        "moderated_variance": np.asarray(
            core_details["moderated_variance"], dtype=np.float64
        ),
        "standard_error": np.asarray(core_details["standard_error"], dtype=np.float64),
        "z_scores": np.asarray(core_details["z_scores"], dtype=np.float64),
        "min_effect": float(core_details["min_effect"]),
        "precision_floor_ratio": float(core_details["precision_floor_ratio"]),
        "prior_spatial_window": int(core_details["prior_spatial_window"]),
        "degrees_of_freedom": np.asarray(
            core_details["degrees_of_freedom"], dtype=np.float64
        ),
        "centered_matrix": centered,
    }
    return scores.astype(np.float64), details


def _standardize_wls_z_scores(
    z_scores: np.ndarray,
) -> tuple[np.ndarray, float]:
    z_scores_ = np.asarray(z_scores, dtype=np.float64)
    if z_scores_.ndim != 1:
        raise ValueError("`z_scores` must be one-dimensional")
    finite = np.isfinite(z_scores_)
    if not np.any(finite):
        raise ValueError("`z_scores` must contain at least one finite value")
    z_scores_finite = z_scores_[finite]
    negative = z_scores_finite[z_scores_finite <= 0.0]
    null_scale = (
        _robust_scale(np.concatenate((negative, -negative)))
        if negative.size > 0
        else _robust_scale(z_scores_finite)
    )
    standardized = np.zeros_like(z_scores_, dtype=np.float64)
    # we standardize by empirical null scale, where null scale is
    # from the negative residuals so the null is always centered
    # at zero and the standardized scores are one-sided pos. exceedances
    standardized[finite] = z_scores_[finite] / max(null_scale, 1.0e-6)
    return standardized, float(null_scale)


def candidate_mask_from_wls(
    z_scores: np.ndarray,
    tail_z: float = 2.0,
    min_signal: float = 0.0,
) -> np.ndarray:
    r"""Turn WLS z-scores into a one-sided exceedance mask.

    ROCCO first estimates a robust null width from the non-positive side,
    rescales the z-scores by that width, and then marks loci with
    :math:`\tilde z_j > z_0`.
    """
    z_scores_ = np.asarray(z_scores, dtype=np.float64)
    standardized, _ = _standardize_wls_z_scores(z_scores_)
    mask = standardized > float(tail_z)
    if min_signal > 0:
        mask &= z_scores_ > float(min_signal)
    return mask


def _estimate_effective_sample_size(
    values: np.ndarray,
    max_lag: int,
) -> tuple[float, float, int]:
    r"""Estimate ESS from the integrated autocorrelation time of a 1D series.

    ROCCO uses

    .. math::

       n_{\mathrm{eff}} = n / \tau_{\mathrm{int}},
       \qquad
       \tau_{\mathrm{int}} = 1 + 2 \sum_{k \ge 1} \rho_k

    with Geyer's initial-positive-sequence truncation so the variance inflation
    factor stays nonnegative.
    """
    values_ = np.asarray(values, dtype=np.float64)
    if values_.ndim != 1:
        raise ValueError("`values` must be one-dimensional")
    n_loci = int(values_.size)
    if n_loci < 4:
        return float(max(1, n_loci)), 1.0, 0

    centered = values_ - float(np.mean(values_))
    var0 = float(np.mean(centered * centered))
    if not np.isfinite(var0) or var0 <= 1.0e-12:
        return float(n_loci), 1.0, 0

    max_lag_ = int(min(max(2, max_lag), n_loci - 1))
    n_fft = 1 << int(np.ceil(np.log2((2 * n_loci) - 1)))
    spectrum = np.fft.rfft(centered, n=n_fft)
    acov = np.fft.irfft(
        spectrum * np.conjugate(spectrum),
        n=n_fft,
    )[: max_lag_ + 1]
    acov /= np.arange(n_loci, n_loci - max_lag_ - 1, -1, dtype=np.float64)
    if not np.isfinite(acov[0]) or acov[0] <= 1.0e-12:
        return float(n_loci), 1.0, 0

    acf = np.clip(acov[1:] / acov[0], -1.0, 1.0)
    tau_int = 1.0
    lags_used = 0
    # Geyer keeps :math:`\tau_{\mathrm{int}}` on the valid side by
    # summing adjacent autocorrelation pairs until the first nonpositive pair
    for lag_idx in range(0, int(acf.size), 2):
        rho_k = float(acf[lag_idx])
        rho_k1 = float(acf[lag_idx + 1]) if (lag_idx + 1) < acf.size else 0.0
        pair_sum = rho_k + rho_k1
        if not np.isfinite(pair_sum) or pair_sum <= 0.0:
            break
        tau_int += 2.0 * pair_sum
        lags_used = int(min(max_lag_, lag_idx + 2))

    effective_n = float(np.clip(n_loci / max(tau_int, 1.0), 1.0, n_loci))
    return effective_n, float(tau_int), int(lags_used)


def _halfSampleMode(sortedValues: np.ndarray) -> float:
    values = np.asarray(sortedValues, dtype=np.float64)
    size = int(values.size)
    if size == 0:
        raise ValueError("`sortedValues` must be non-empty")
    if size == 1:
        return float(values[0])
    if size == 2:
        return float(np.mean(values))
    if size == 3:
        leftWidth = float(values[1] - values[0])
        rightWidth = float(values[2] - values[1])
        if leftWidth <= rightWidth:
            return float(np.mean(values[:2]))
        return float(np.mean(values[1:]))

    window = int(np.ceil(size / 2.0))
    bestStart = 0
    bestWidth = float(values[window - 1] - values[0])
    for start in range(1, size - window + 1):
        width = float(values[start + window - 1] - values[start])
        if width < bestWidth:
            bestWidth = width
            bestStart = start
    return _halfSampleMode(values[bestStart : bestStart + window])


def _estimateStationaryNull(scoreTrack: np.ndarray) -> tuple[float, float]:
    scores = np.asarray(scoreTrack, dtype=np.float64)
    bulkCutoff = float(
        np.quantile(
            scores,
            0.60,
            method="interpolated_inverted_cdf",
        )
    )
    lowerBulk = scores[scores <= bulkCutoff]
    lowerBulk = np.sort(lowerBulk)
    nullCenter = _halfSampleMode(lowerBulk)
    bulkResiduals = lowerBulk - nullCenter
    bulkMAD = 1.4826 * float(
        np.median(np.abs(bulkResiduals - np.median(bulkResiduals)))
    )
    bulkIQR = float(stats.iqr(bulkResiduals, rng=(25, 75))) / 1.349
    bulkSD = float(np.std(bulkResiduals, ddof=1))
    centerScale = float(max(bulkMAD, bulkIQR, bulkSD, 1.0e-6))

    for _ in range(20):
        standardized = (scores - nullCenter) / (3.0 * centerScale)
        centerMask = np.abs(standardized) < 1.0
        if np.count_nonzero(centerMask) < 8:
            raise ValueError("Stationary null center requires eight central residuals")
        centerWeights = (1.0 - standardized[centerMask] ** 2) ** 2
        nextCenter = float(
            np.dot(centerWeights, scores[centerMask]) / np.sum(centerWeights)
        )
        if abs(nextCenter - nullCenter) <= 1.0e-10 * max(centerScale, 1.0):
            nullCenter = nextCenter
            break
        nullCenter = nextCenter

    lowerMagnitudes = nullCenter - scores[scores < nullCenter]
    if lowerMagnitudes.size < 8:
        raise ValueError("Stationary null estimation requires eight lower residuals")
    nullScale = float(np.median(lowerMagnitudes) / stats.norm.ppf(0.75))
    if nullScale <= _TINY:
        raise ValueError("Stationary lower-residual scale must be positive")
    return float(nullCenter), float(nullScale)


def _prepareStationaryNullTemplate(
    scoreTrack: np.ndarray,
    nullCenter: float,
) -> tuple[np.ndarray, dict[str, float | int]]:
    centered = np.asarray(scoreTrack, dtype=np.float64) - float(nullCenter)
    lowerMagnitudes = -centered[centered < 0.0]
    if lowerMagnitudes.size < 8:
        raise ValueError("Stationary null template requires eight lower residuals")
    lowerScale = float(np.median(lowerMagnitudes) / stats.norm.ppf(0.75))
    if lowerScale <= _TINY:
        raise ValueError("Stationary lower-residual scale must be positive")

    template = centered.copy()
    positiveMask = centered > 0.0
    positiveCount = int(np.count_nonzero(positiveMask))
    if positiveCount:
        positiveRanks = stats.rankdata(centered[positiveMask], method="average")
        positiveQuantiles = (positiveRanks - 0.5) / float(positiveCount)
        sortedLowerMagnitudes = np.sort(lowerMagnitudes)
        virtualIndices = positiveQuantiles * float(lowerMagnitudes.size) - 1.0
        lowerIndices = np.floor(virtualIndices).astype(np.int64)
        fractions = virtualIndices - lowerIndices
        upperIndices = np.clip(
            lowerIndices + 1,
            0,
            lowerMagnitudes.size - 1,
        )
        lowerIndices = np.clip(
            lowerIndices,
            0,
            lowerMagnitudes.size - 1,
        )
        template[positiveMask] = (
            sortedLowerMagnitudes[lowerIndices] * (1.0 - fractions)
            + sortedLowerMagnitudes[upperIndices] * fractions
        )

    lowerCap, upperCap = np.quantile(
        template,
        (0.001, 0.999),
        method="interpolated_inverted_cdf",
    )
    template = np.clip(template, lowerCap, upperCap)
    template -= float(np.mean(template))
    templateRMS = float(np.sqrt(np.mean(template * template)))
    if templateRMS <= _TINY:
        raise ValueError("Stationary null template scale must be positive")
    template *= lowerScale / templateRMS

    metadata: dict[str, float | int] = {
        "templateLowerSize": int(lowerMagnitudes.size),
        "templateClipLower": float(lowerCap),
        "templateClipUpper": float(upperCap),
    }
    return np.ascontiguousarray(template, dtype=np.float64), metadata


def estimateStationaryBootstrapBudget(
    scoreTrack: np.ndarray,
    bootstrapBlockLength: int,
    *,
    stepBP: int = 50,
    thresholdZ: float = 2.0,
    numBootstrap: int = 64,
    randomSeed: int = 42,
    useLocalBootstrapRadius: bool = True,
    returnDetails: bool = False,
) -> float | tuple[float, dict[str, Any]]:
    scores = np.asarray(scoreTrack, dtype=np.float64)
    if scores.ndim != 1:
        raise ValueError("`scoreTrack` must be one-dimensional")
    if scores.size == 0:
        raise ValueError("`scoreTrack` must be non-empty")
    if np.any(~np.isfinite(scores)):
        raise ValueError("`scoreTrack` values must be finite")
    if isinstance(bootstrapBlockLength, (bool, np.bool_)) or not isinstance(
        bootstrapBlockLength,
        (int, np.integer),
    ):
        raise ValueError("`bootstrapBlockLength` must be a positive integer")
    blockLength = int(bootstrapBlockLength)
    if blockLength <= 0:
        raise ValueError("`bootstrapBlockLength` must be a positive integer")
    if isinstance(stepBP, (bool, np.bool_)) or not isinstance(
        stepBP,
        (int, np.integer),
    ):
        raise ValueError("`stepBP` must be a positive integer")
    stepBPValue = int(stepBP)
    if stepBPValue <= 0:
        raise ValueError("`stepBP` must be a positive integer")
    if isinstance(thresholdZ, (bool, np.bool_)):
        raise ValueError("`thresholdZ` must be finite and non-negative")
    thresholdZValue = float(thresholdZ)
    if not np.isfinite(thresholdZValue) or thresholdZValue < 0.0:
        raise ValueError("`thresholdZ` must be finite and non-negative")
    if isinstance(numBootstrap, (bool, np.bool_)) or not isinstance(
        numBootstrap,
        (int, np.integer),
    ):
        raise ValueError("`numBootstrap` must be an integer of at least 8")
    numBootstrapValue = int(numBootstrap)
    if numBootstrapValue < 8:
        raise ValueError("`numBootstrap` must be an integer of at least 8")
    if isinstance(randomSeed, (bool, np.bool_)) or not isinstance(
        randomSeed,
        (int, np.integer),
    ):
        raise ValueError("`randomSeed` must be a non-negative integer")
    randomSeedValue = int(randomSeed)
    if randomSeedValue < 0:
        raise ValueError("`randomSeed` must be a non-negative integer")
    if not isinstance(useLocalBootstrapRadius, (bool, np.bool_)):
        raise ValueError("`useLocalBootstrapRadius` must be boolean")
    if not isinstance(returnDetails, (bool, np.bool_)):
        raise ValueError("`returnDetails` must be boolean")

    nullCenter, nullScale = _estimateStationaryNull(scores)
    template, templateMetadata = _prepareStationaryNullTemplate(
        scores,
        nullCenter,
    )
    tailAlpha = float(stats.norm.sf(thresholdZValue))
    thresholdQuantile = float(stats.norm.cdf(thresholdZValue))
    thresholdOffset = float(
        max(
            np.quantile(
                template,
                thresholdQuantile,
                method="interpolated_inverted_cdf",
            ),
            0.0,
        )
    )
    threshold = float(nullCenter + thresholdOffset)
    tailTrack = scores > threshold
    trackTailOccupancy = float(np.mean(tailTrack))

    useLocalRadius = bool(useLocalBootstrapRadius)
    maxLocalRadiusIntervals = int(1_000_000 // stepBPValue) if useLocalRadius else -1
    if useLocalRadius:
        if blockLength >= scores.size:
            localRadiusIntervals = int(scores.size - 1)
        else:
            localRadiusIntervals = int(
                min(
                    np.ceil(np.sqrt(float(blockLength * scores.size))),
                    scores.size - 1,
                )
            )
        localRadiusIntervals = min(
            localRadiusIntervals,
            maxLocalRadiusIntervals,
        )
    else:
        localRadiusIntervals = None

    rng = np.random.default_rng(randomSeedValue)
    nullOccupancies = np.empty(numBootstrapValue, dtype=np.float64)
    for drawIndex in range(numBootstrapValue):
        draw = np.asarray(
            _wls.stationaryNullBootstrapDraw(
                template,
                blockLength,
                rng,
                maxLocalRadiusIntervals,
            ),
            dtype=np.float64,
        )
        nullOccupancies[drawIndex] = float(np.mean(draw > thresholdOffset))

    nullTailOccupancy = float(np.mean(nullOccupancies))
    nullTailOccupancySD = float(np.std(nullOccupancies, ddof=1))
    nullTailMCSE = float(nullTailOccupancySD / np.sqrt(float(numBootstrapValue)))
    signedTailExcess = float(trackTailOccupancy - nullTailOccupancy)
    budgetFraction = float(max(signedTailExcess, 0.0))

    essMaxLag = int(min(scores.size - 1, blockLength))
    effectiveTotalCount, autocorrelationTime, essLagsUsed = (
        _estimate_effective_sample_size(
            tailTrack.astype(np.float64),
            max_lag=essMaxLag,
        )
    )
    effectiveCount = float(budgetFraction * effectiveTotalCount)

    if not bool(returnDetails):
        return budgetFraction

    details: dict[str, Any] = {
        "bootstrapMethod": "stationary_bootstrap",
        "bootstrapBlockLength": int(blockLength),
        "numBootstrap": int(numBootstrapValue),
        "randomSeed": int(randomSeedValue),
        "stepBP": int(stepBPValue),
        "useLocalBootstrapRadius": bool(useLocalRadius),
        "localRadiusLimitBP": 1_000_000,
        "maxLocalRadiusIntervals": int(maxLocalRadiusIntervals),
        "localRadiusIntervals": localRadiusIntervals,
        "nullCenter": float(nullCenter),
        "nullScale": float(nullScale),
        **templateMetadata,
        "thresholdZ": float(thresholdZValue),
        "tailAlpha": float(tailAlpha),
        "thresholdOffset": float(thresholdOffset),
        "threshold": float(threshold),
        "trackTailOccupancy": float(trackTailOccupancy),
        "nullTailOccupancy": float(nullTailOccupancy),
        "nullTailOccupancySD": float(nullTailOccupancySD),
        "nullTailMCSE": float(nullTailMCSE),
        "signedTailExcess": float(signedTailExcess),
        "budgetFraction": float(budgetFraction),
        "effectiveCount": float(effectiveCount),
        "effectiveTotalCount": float(effectiveTotalCount),
        "autocorrelationTime": float(autocorrelationTime),
        "essMaxLag": int(essMaxLag),
        "essLagsUsed": int(essLagsUsed),
        "numLoci": int(scores.size),
    }
    return budgetFraction, details


def fit_beta_prior_mle(
    successes: np.ndarray,
    totals: np.ndarray,
    init_center: float = 0.05,
    init_strength: float = 10.0,
) -> Tuple[float, float]:
    successes_ = np.asarray(successes, dtype=np.float64)
    totals_ = np.asarray(totals, dtype=np.float64)
    if successes_.shape != totals_.shape:
        raise ValueError("`successes` and `totals` must have the same shape")
    if successes_.size == 0:
        return 1.0, 1.0

    init_center_ = min(max(float(init_center), 1.0e-6), 1.0 - 1.0e-6)
    raw_rates = successes_ / np.maximum(totals_, 1.0)
    pooled_rate = float(
        np.clip(
            np.sum(successes_) / max(np.sum(totals_), 1.0),
            1.0e-6,
            1.0 - 1.0e-6,
        )
    )
    observed_raw_rate_var = (
        float(np.var(raw_rates, ddof=1)) if raw_rates.size > 1 else 0.0
    )
    # Beta-binomial: the minimum structural dispersion is the
    # binomial boundary :math:`\rho = 0`, which gives
    # :math:`\operatorname{Var}(x_c/n_c) = p(1-p)/n_c`, but we use an ESS s.t. this becomes the following floor
    # on the observed variance of raw rates across chromosomes.

    theoretical_min_raw_rate_var = float(
        pooled_rate * (1.0 - pooled_rate) * np.mean(1.0 / np.maximum(totals_, 1.0))
    )
    if observed_raw_rate_var <= theoretical_min_raw_rate_var + 1.0e-12:
        boundary_strength = float(max(1.0e12, 100.0 * np.max(totals_)))
        return (
            pooled_rate * boundary_strength,
            (1.0 - pooled_rate) * boundary_strength,
        )

    def objective(theta: np.ndarray) -> float:
        alpha = float(np.exp(theta[0]))
        beta = float(np.exp(theta[1]))
        loglik = np.sum(
            special.betaln(successes_ + alpha, totals_ - successes_ + beta)
            - special.betaln(alpha, beta)
        )
        return float(-loglik)

    init = np.log(
        np.array(
            [
                init_center_ * float(init_strength),
                (1.0 - init_center_) * float(init_strength),
            ],
            dtype=np.float64,
        )
    )
    result = optimize.minimize(
        objective,
        init,
        method="L-BFGS-B",
    )
    if not result.success:
        logger.warning(
            "Falling back to a weak beta prior while fitting EB budgets: %s",
            result.message,
        )
        return (
            init_center_ * float(init_strength),
            (1.0 - init_center_) * float(init_strength),
        )
    alpha_hat = float(np.exp(result.x[0]))
    beta_hat = float(np.exp(result.x[1]))
    return alpha_hat, beta_hat


def _beta_posterior_budget_quantile(
    successes: float,
    total: float,
    alpha: float,
    beta: float,
    posterior_quantile: float,
    min_budget: float,
    max_budget: float,
) -> float:
    posterior_alpha = float(max(1.0e-12, successes + alpha))
    posterior_beta = float(max(1.0e-12, (total - successes) + beta))
    posterior_quantile_ = float(np.clip(posterior_quantile, 1.0e-6, 1.0 - 1.0e-6))
    posterior_budget = float(
        stats.beta.ppf(
            posterior_quantile_,
            posterior_alpha,
            posterior_beta,
        )
    )
    return float(
        np.clip(
            posterior_budget,
            min_budget,
            max_budget,
        )
    )


def estimate_empirical_bayes_budgets(
    chrom_candidate_counts: Dict[str, float],
    chrom_total_counts: Dict[str, float],
    min_budget: float = 1.0e-4,
    max_budget: float = 0.5,
    init_center: float = 0.05,
    init_strength: float = 10.0,
    posterior_quantile: float = 0.01,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    r"""Estimate per-chromosome budgets with beta-binomial EB shrinkage."""
    chroms = list(chrom_candidate_counts.keys())
    if chroms != list(chrom_total_counts.keys()):
        raise ValueError(
            "`chrom_candidate_counts` and `chrom_total_counts` must share keys in the same order"
        )

    successes = np.array(
        [chrom_candidate_counts[chrom] for chrom in chroms],
        dtype=np.float64,
    )
    totals = np.array(
        [chrom_total_counts[chrom] for chrom in chroms],
        dtype=np.float64,
    )
    raw_budgets = successes / np.maximum(totals, 1.0)
    pooled_rate = float(
        np.clip(
            np.sum(successes) / max(np.sum(totals), 1.0),
            1.0e-6,
            1.0 - 1.0e-6,
        )
    )
    observed_raw_budget_var = (
        float(np.var(raw_budgets, ddof=1)) if raw_budgets.size > 1 else 0.0
    )
    theoretical_min_raw_budget_var = float(
        pooled_rate * (1.0 - pooled_rate) * np.mean(1.0 / np.maximum(totals, 1.0))
    )
    dispersion_at_floor = bool(
        observed_raw_budget_var <= theoretical_min_raw_budget_var + 1.0e-12
    )

    posterior_quantile_ = float(posterior_quantile)
    if not (0.0 < posterior_quantile_ < 1.0):
        raise ValueError("`posterior_quantile` must lie strictly between 0 and 1")

    if len(chroms) <= 1:
        alpha0 = float(init_center) * float(init_strength)
        beta0 = (1.0 - float(init_center)) * float(init_strength)
        shrunk = {
            chrom: float(
                _beta_posterior_budget_quantile(
                    successes[idx],
                    totals[idx],
                    alpha0,
                    beta0,
                    posterior_quantile_,
                    min_budget,
                    max_budget,
                )
            )
            for idx, chrom in enumerate(chroms)
        }
        return shrunk, {
            "alpha": float(alpha0),
            "beta": float(beta0),
            "genomeWideBudgetFraction": float(init_center),
            "priorStrength": float(init_strength),
            "priorDispersion": float(1.0 / (1.0 + alpha0 + beta0)),
            "minimumPriorDispersion": 0.0,
            "observedRawBudgetVariance": float(observed_raw_budget_var),
            "theoreticalMinimumRawBudgetVariance": float(
                theoretical_min_raw_budget_var
            ),
            "priorDispersionAtFloor": bool(False),
            "posteriorSummary": "beta_quantile",
            "posteriorQuantile": float(posterior_quantile_),
            "priorFitMethod": "single_chrom_default",
        }

    if len(chroms) <= 3:
        alpha_hat = float(pooled_rate) * float(init_strength)
        beta_hat = (1.0 - float(pooled_rate)) * float(init_strength)
        shrunk = {
            chrom: _beta_posterior_budget_quantile(
                successes[idx],
                totals[idx],
                alpha_hat,
                beta_hat,
                posterior_quantile_,
                min_budget,
                max_budget,
            )
            for idx, chrom in enumerate(chroms)
        }
        prior_strength = float(alpha_hat + beta_hat)
        return shrunk, {
            "alpha": float(alpha_hat),
            "beta": float(beta_hat),
            "genomeWideBudgetFraction": float(pooled_rate),
            "priorStrength": float(prior_strength),
            "priorDispersion": float(max(0.0, 1.0 / (1.0 + prior_strength))),
            "minimumPriorDispersion": 0.0,
            "observedRawBudgetVariance": float(observed_raw_budget_var),
            "theoreticalMinimumRawBudgetVariance": float(
                theoretical_min_raw_budget_var
            ),
            "priorDispersionAtFloor": bool(
                observed_raw_budget_var <= theoretical_min_raw_budget_var + 1.0e-12
            ),
            "posteriorSummary": "beta_quantile",
            "posteriorQuantile": float(posterior_quantile_),
            "priorFitMethod": "weak_pooled_prior",
        }

    alpha_hat, beta_hat = fit_beta_prior_mle(
        successes,
        totals,
        init_center=init_center,
        init_strength=init_strength,
    )
    shrunk = {}
    for idx, chrom in enumerate(chroms):
        shrunk[chrom] = _beta_posterior_budget_quantile(
            successes[idx],
            totals[idx],
            alpha_hat,
            beta_hat,
            posterior_quantile_,
            min_budget,
            max_budget,
        )
    prior_strength = float(alpha_hat + beta_hat)
    prior_dispersion = float(max(0.0, 1.0 / (1.0 + prior_strength)))

    return shrunk, {
        "alpha": float(alpha_hat),
        "beta": float(beta_hat),
        "genomeWideBudgetFraction": float(alpha_hat / (alpha_hat + beta_hat)),
        "priorStrength": float(prior_strength),
        "priorDispersion": float(prior_dispersion),
        "minimumPriorDispersion": 0.0,
        "observedRawBudgetVariance": float(observed_raw_budget_var),
        "theoreticalMinimumRawBudgetVariance": float(
            theoretical_min_raw_budget_var
        ),
        "priorDispersionAtFloor": bool(dispersion_at_floor),
        "posteriorSummary": "beta_quantile",
        "posteriorQuantile": float(posterior_quantile_),
        "priorFitMethod": "beta_binomial_mle",
    }
