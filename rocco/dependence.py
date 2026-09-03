from __future__ import annotations

import heapq
import math
import warnings
from collections.abc import Mapping
from numbers import Integral, Real
from typing import Any

import numpy as np

_MIN_WINDOWS = 20
_MIN_AUTOSOMES = 4
_MAX_LAG_BP = 50000
MIN_CORRELATION_RADIUS_BP = 2500
_PRIOR_RADIUS_BP = 5000
_ACF_THRESHOLD = 0.1
_ACF_REQUIRED_CROSSINGS = 5
_CLIP_QUANTILES = (0.005, 0.995)


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"`{name}` must be an integer")
    result = int(value)
    if result <= 0:
        raise ValueError(f"`{name}` must be positive")
    return result


def _nonnegative_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"`{name}` must be an integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"`{name}` must be nonnegative")
    return result


def _probability(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"`{name}` must be a real number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0 or result >= 1.0:
        raise ValueError(f"`{name}` must lie in (0, 1)")
    return result


def _prepare_track_acfs(
    track_matrix: np.ndarray,
    max_lag_bins: int,
    clip_quantiles: tuple[float, float],
    min_finite_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_tracks, n_bins = track_matrix.shape
    n_fft = 1 << int((2 * n_bins - 2).bit_length())
    track_acfs = np.full((n_tracks, max_lag_bins + 1), np.nan, dtype=np.float64)
    track_pair_counts = np.full_like(track_acfs, np.nan)
    lower_quantile, upper_quantile = clip_quantiles

    for track_idx in range(n_tracks):
        row = track_matrix[track_idx]
        finite = np.isfinite(row)
        if int(np.sum(finite)) < min_finite_count:
            continue
        finite_values = row[finite]
        lower = float(np.quantile(finite_values, lower_quantile))
        upper = float(np.quantile(finite_values, upper_quantile))
        clipped = np.clip(finite_values, lower, upper)
        clipped -= float(np.mean(clipped))
        values = np.zeros(n_bins, dtype=np.float64)
        values[finite] = clipped
        mask = finite.astype(np.float64, copy=False)

        value_fft = np.fft.rfft(values, n=n_fft)
        mask_fft = np.fft.rfft(mask, n=n_fft)
        numerators = np.fft.irfft(value_fft * np.conjugate(value_fft), n=n_fft)[
            : max_lag_bins + 1
        ]
        pair_counts = np.rint(
            np.maximum(
                np.fft.irfft(mask_fft * np.conjugate(mask_fft), n=n_fft)[
                    : max_lag_bins + 1
                ],
                0.0,
            )
        )
        autocovariance = np.full(max_lag_bins + 1, np.nan, dtype=np.float64)
        np.divide(
            numerators,
            pair_counts,
            out=autocovariance,
            where=pair_counts >= float(min_finite_count),
        )
        variance = float(autocovariance[0])
        if not math.isfinite(variance) or variance <= np.finfo(np.float64).tiny:
            continue
        track_acfs[track_idx] = np.clip(autocovariance / variance, -1.0, 1.0)
        track_pair_counts[track_idx] = pair_counts
    return track_acfs, track_pair_counts


def _median_across_tracks(values: np.ndarray) -> np.ndarray:
    if np.all(np.isfinite(values)):
        return np.median(values, axis=-2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=-2)


def _smooth_acfs(acfs: np.ndarray) -> np.ndarray:
    finite = np.isfinite(acfs)
    padded_values = np.pad(
        np.where(finite, acfs, 0.0), ((0, 0), (2, 2)), mode="constant"
    )
    padded_counts = np.pad(finite.astype(np.int64), ((0, 0), (2, 2)))
    value_prefix = np.pad(
        np.cumsum(padded_values, axis=1), ((0, 0), (1, 0)), mode="constant"
    )
    count_prefix = np.pad(
        np.cumsum(padded_counts, axis=1), ((0, 0), (1, 0)), mode="constant"
    )
    window_sums = value_prefix[:, 5:] - value_prefix[:, :-5]
    window_counts = count_prefix[:, 5:] - count_prefix[:, :-5]
    smoothed = np.full_like(acfs, np.nan)
    np.divide(window_sums, window_counts, out=smoothed, where=window_counts > 0)
    return smoothed


def _crossings_from_pooled_acfs(
    pooled_acfs: np.ndarray,
    step_bp: int,
    threshold: float,
    required_crossings: int,
    include_smoothed: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pooled = np.asarray(pooled_acfs, dtype=np.float64)
    if pooled.ndim == 1:
        pooled = pooled[None, :]
    smoothed = _smooth_acfs(np.abs(pooled))
    lag_values = smoothed[:, 1:]
    below_threshold = np.isfinite(lag_values) & (np.abs(lag_values) < threshold)
    below_prefix = np.pad(
        np.cumsum(below_threshold, axis=1), ((0, 0), (1, 0)), mode="constant"
    )
    run_counts = (
        below_prefix[:, required_crossings:] - below_prefix[:, :-required_crossings]
    )
    qualifying = run_counts == required_crossings
    has_crossing = np.any(qualifying, axis=1)
    crossing_bins = np.argmax(qualifying, axis=1) + 1
    max_lag_bins = pooled.shape[1] - 1
    raw_lag_bp = np.full(pooled.shape[0], max_lag_bins * step_bp, dtype=np.int64)
    raw_lag_bp[has_crossing] = crossing_bins[has_crossing] * step_bp
    valid = np.isfinite(pooled[:, 0]) & (
        np.sum(np.isfinite(lag_values), axis=1) >= required_crossings
    )
    censored = ~has_crossing
    gaussian_factor = 3.0 / (2.0 * math.sqrt(-math.log(threshold)))
    radii_bp = raw_lag_bp.astype(np.float64) * gaussian_factor
    return radii_bp, censored, valid, smoothed if include_smoothed else None


def _periodicity_diagnostics(
    smoothed_acf: np.ndarray,
    step_bp: int,
    crossing_lag_bp: int,
    right_censored: bool,
    required_crossings: int,
) -> tuple[float | None, float, float]:
    values = np.asarray(smoothed_acf[1:], dtype=np.float64)
    finite = np.isfinite(values)
    if int(np.sum(finite)) < 8:
        return None, 0.0, 0.0
    indices = np.arange(values.size, dtype=np.float64)
    design = np.column_stack((indices[finite], np.ones(int(np.sum(finite)))))
    slope, intercept = np.linalg.lstsq(design, values[finite], rcond=None)[0]
    residual = np.zeros_like(values)
    residual[finite] = values[finite] - ((slope * indices[finite]) + intercept)
    powers = np.square(np.abs(np.fft.rfft(residual)))
    if powers.size <= 1:
        dominant_period_bp = None
        oscillation_strength = 0.0
    else:
        powers[0] = 0.0
        total_power = float(np.sum(powers))
        dominant_idx = int(np.argmax(powers))
        if total_power <= np.finfo(np.float64).tiny or dominant_idx == 0:
            dominant_period_bp = None
            oscillation_strength = 0.0
        else:
            dominant_period_bp = float((values.size * step_bp) / float(dominant_idx))
            oscillation_strength = float(powers[dominant_idx] / total_power)

    if right_censored:
        revival_start = max(0, (3 * values.size) // 4)
    else:
        revival_start = min(
            values.size,
            max(0, int(crossing_lag_bp // step_bp) + required_crossings - 1),
        )
    tail = values[revival_start:]
    tail = tail[np.isfinite(tail)]
    post_crossing_revival = float(np.max(np.abs(tail))) if tail.size else 0.0
    return dominant_period_bp, oscillation_strength, post_crossing_revival


def _serializable_float_list(values: np.ndarray) -> list[float | None]:
    return [float(value) if math.isfinite(float(value)) else None for value in values]


def _window_result(
    track_acfs: np.ndarray,
    track_pair_counts: np.ndarray,
    pooled_acf: np.ndarray,
    smoothed_acf: np.ndarray,
    radius_bp: float,
    right_censored: bool,
    step_bp: int,
    threshold: float,
    required_crossings: int,
) -> dict[str, Any]:
    gaussian_factor = 3.0 / (2.0 * math.sqrt(-math.log(threshold)))
    censor_bound_bp = int(round(float(radius_bp) / gaussian_factor))
    raw_crossing_lag_bp = None if right_censored else censor_bound_bp
    pooled_pairs = _median_across_tracks(track_pair_counts)
    finite_pairs = pooled_pairs[np.isfinite(pooled_pairs) & (pooled_pairs > 0.0)]
    period_bp, oscillation_strength, revival = _periodicity_diagnostics(
        pooled_acf,
        step_bp,
        censor_bound_bp,
        right_censored,
        required_crossings,
    )
    return {
        "rawCrossingLagBP": raw_crossing_lag_bp,
        "gaussianEquivalentRadiusBP": float(radius_bp),
        "rightCensored": bool(right_censored),
        "lagsBP": [int(lag * step_bp) for lag in range(pooled_acf.size)],
        "ACF": _serializable_float_list(pooled_acf),
        "smoothedACF": _serializable_float_list(smoothed_acf),
        "finitePairCounts": _serializable_float_list(pooled_pairs),
        "finitePairCountMin": float(np.min(finite_pairs)),
        "finitePairCountMedian": float(np.median(finite_pairs)),
        "finitePairCountMax": float(np.max(finite_pairs)),
        "trackCount": int(track_acfs.shape[0]),
        "validTrackCount": int(np.sum(np.isfinite(track_acfs[:, 0]))),
        "dominantACFPeriodBP": period_bp,
        "oscillationStrength": oscillation_strength,
        "postCrossingACFRevival": revival,
    }


def estimate_dependence_radius_for_window(
    values: np.ndarray,
    step_bp: int,
    *,
    max_lag_bp: int = _MAX_LAG_BP,
    acf_threshold: float = _ACF_THRESHOLD,
    acf_required_crossings: int = _ACF_REQUIRED_CROSSINGS,
    clip_quantiles: tuple[float, float] = _CLIP_QUANTILES,
) -> dict[str, Any]:
    step = _positive_int(step_bp, "step_bp")
    max_lag = _positive_int(max_lag_bp, "max_lag_bp")
    required = _positive_int(acf_required_crossings, "acf_required_crossings")
    threshold = _probability(acf_threshold, "acf_threshold")
    if len(clip_quantiles) != 2:
        raise ValueError("`clip_quantiles` must contain two probabilities")
    lower_quantile = float(clip_quantiles[0])
    upper_quantile = float(clip_quantiles[1])
    if (
        not math.isfinite(lower_quantile)
        or not math.isfinite(upper_quantile)
        or lower_quantile < 0.0
        or upper_quantile > 1.0
        or lower_quantile >= upper_quantile
    ):
        raise ValueError("`clip_quantiles` must be ordered within [0, 1]")

    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("`values` must be a two-dimensional track-by-bin matrix")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("`values` must not be empty")
    max_lag_bins = min(max_lag // step, matrix.shape[1] - 1)
    if max_lag_bins < required:
        raise ValueError("The window contains too few bins for the ACF crossing rule")
    track_acfs, track_pair_counts = _prepare_track_acfs(
        matrix,
        max_lag_bins,
        (lower_quantile, upper_quantile),
        max(20, 2 * required),
    )
    pooled_acf = _median_across_tracks(track_acfs)
    radii_bp, censored, valid, smoothed = _crossings_from_pooled_acfs(
        pooled_acf,
        step,
        threshold,
        required,
        True,
    )
    if not bool(valid[0]):
        raise ValueError("The window lacks a nonconstant finite track")
    if smoothed is None:
        raise RuntimeError("ACF smoothing failed")
    return _window_result(
        track_acfs,
        track_pair_counts,
        pooled_acf,
        smoothed[0],
        float(radii_bp[0]),
        bool(censored[0]),
        step,
        threshold,
        required,
    )


def _is_autosome(chromosome: str) -> bool:
    label = chromosome[3:] if chromosome.startswith("chr") else chromosome
    return label.isdecimal() and int(label) > 0


def _kaplan_meier_quantile(
    durations: np.ndarray, right_censored: np.ndarray, quantile: float
) -> float | None:
    values = np.asarray(durations, dtype=np.float64)
    censored = np.asarray(right_censored, dtype=bool)
    if values.ndim != 1 or censored.shape != values.shape or values.size == 0:
        raise ValueError("Kaplan-Meier inputs must be nonempty aligned vectors")
    order = np.argsort(values, kind="mergesort")
    ordered_values = values[order]
    ordered_censored = censored[order]
    survival = 1.0
    at_risk = int(values.size)
    target_survival = 1.0 - quantile
    for time in np.unique(ordered_values):
        tied = ordered_values == time
        events = int(np.sum(~ordered_censored[tied]))
        removals = int(np.sum(tied))
        if events > 0:
            survival *= 1.0 - (float(events) / float(at_risk))
            if survival <= target_survival:
                return float(time)
        at_risk -= removals
    return None


def _candidate_windows(
    matrices: Mapping[str, np.ndarray],
    coordinates: Mapping[str, np.ndarray],
    step_bp: int,
    window_bp: int,
    window_count: int,
) -> tuple[list[tuple[float, str, int, int, int]], int]:
    window_bins = window_bp // step_bp
    candidates: list[
        tuple[tuple[float, int, int], tuple[float, str, int, int, int]]
    ] = []
    eligible_count = 0
    sorted_chromosomes = sorted(matrices)
    chromosome_ranks = {
        chromosome: rank for rank, chromosome in enumerate(sorted_chromosomes)
    }
    for chromosome in sorted_chromosomes:
        if not _is_autosome(chromosome):
            continue
        source_matrix = matrices[chromosome]
        matrix = np.asarray(source_matrix, dtype=np.float64)
        starts = coordinates[chromosome]
        coordinate_offset = int(starts[0])
        first_tile_start = int(
            ((coordinate_offset + window_bp - 1) // window_bp) * window_bp
        )
        first_offset_bp = first_tile_start - coordinate_offset
        if first_offset_bp % step_bp == 0:
            first_start = first_offset_bp // step_bp
            tile_count = (starts.size - first_start) // window_bins
            tiles_per_chunk = 64
            for first_tile_index in range(0, tile_count, tiles_per_chunk):
                chunk_tile_count = min(
                    tiles_per_chunk,
                    tile_count - first_tile_index,
                )
                chunk_start = first_start + (first_tile_index * window_bins)
                chunk_stop = chunk_start + (chunk_tile_count * window_bins)
                track_chunk = matrix[:, chunk_start:chunk_stop]
                positive_chunk = np.maximum(track_chunk, 0.0)
                positive_chunk[~np.isfinite(positive_chunk)] = 0.0
                chunk_scores = np.sum(
                    positive_chunk.reshape(
                        matrix.shape[0],
                        chunk_tile_count,
                        window_bins,
                    ),
                    axis=(0, 2),
                )
                for tile_offset, raw_score in enumerate(chunk_scores):
                    score = float(raw_score)
                    if score <= 0.0:
                        continue
                    tile_index = first_tile_index + tile_offset
                    tile_start = first_tile_start + (tile_index * window_bp)
                    start = first_start + (tile_index * window_bins)
                    stop = start + window_bins
                    eligible_count += 1
                    candidate = (-score, chromosome, tile_start, start, stop)
                    heap_key = (
                        score,
                        -chromosome_ranks[chromosome],
                        -tile_start,
                    )
                    heap_item = (heap_key, candidate)
                    if len(candidates) < window_count:
                        heapq.heappush(candidates, heap_item)
                    elif heap_key > candidates[0][0]:
                        heapq.heapreplace(candidates, heap_item)
        del matrix
        del source_matrix
    selected = sorted(candidate for _, candidate in candidates)
    return selected, eligible_count


def _prior_result(
    reason: str,
    step_bp: int,
    diagnostics: dict[str, Any],
) -> tuple[int, int, int, dict[str, Any]]:
    intervals = int(math.ceil(_PRIOR_RADIUS_BP / float(step_bp)))
    diagnostics.update(
        {
            "usedPrior": True,
            "status": "prior",
            "reason": reason,
            "insufficientDataReason": reason,
            "estimateBP": float(_PRIOR_RADIUS_BP),
            "lowerBP": float(_PRIOR_RADIUS_BP),
            "upperBP": float(_PRIOR_RADIUS_BP),
            "estimateIntervals": intervals,
            "lowerIntervals": intervals,
            "upperIntervals": intervals,
            "workingSpanBP": float(_PRIOR_RADIUS_BP),
            "workingSpanIntervals": intervals,
            "minimumCorrelationRadiusBP": MIN_CORRELATION_RADIUS_BP,
            "minimumCorrelationRadiusApplied": False,
        }
    )
    return intervals, intervals, intervals, diagnostics


def choose_dependence_span(
    chromosome_matrices: Mapping[str, np.ndarray],
    chromosome_coordinates: Mapping[str, np.ndarray],
    step_bp: int,
    *,
    window_bp: int = 50000,
    window_count: int = 256,
    working_quantile: float = 0.95,
    bootstrap_draws: int = 500,
    random_seed: int = 1729,
) -> tuple[int, int, int, dict[str, Any]]:
    step = _positive_int(step_bp, "step_bp")
    window_size = _positive_int(window_bp, "window_bp")
    requested_windows = _positive_int(window_count, "window_count")
    draws = _positive_int(bootstrap_draws, "bootstrap_draws")
    seed = _nonnegative_int(random_seed, "random_seed")
    working_q = _probability(working_quantile, "working_quantile")
    if working_q <= 0.5:
        raise ValueError("`working_quantile` must lie in (0.5, 1)")
    if requested_windows < _MIN_WINDOWS:
        raise ValueError(f"`window_count` must be at least {_MIN_WINDOWS}")
    if draws < 20:
        raise ValueError("`bootstrap_draws` must be at least 20")
    if window_size % step != 0:
        raise ValueError("`window_bp` must be divisible by `step_bp`")
    if not isinstance(chromosome_matrices, Mapping) or not isinstance(
        chromosome_coordinates, Mapping
    ):
        raise TypeError("Chromosome matrices and coordinates must be mappings")
    matrix_keys = set(chromosome_matrices)
    coordinate_keys = set(chromosome_coordinates)
    if not matrix_keys:
        raise ValueError("`chromosome_matrices` must not be empty")
    if matrix_keys != coordinate_keys:
        raise ValueError("Chromosome matrix and coordinate keys must match exactly")
    if any(not isinstance(chromosome, str) for chromosome in matrix_keys):
        raise TypeError("Chromosome names must be strings")

    first_chromosome = min(matrix_keys)
    first_matrix = np.asarray(chromosome_matrices[first_chromosome])
    if first_matrix.ndim != 2 or first_matrix.shape[0] == 0:
        raise ValueError("Chromosome matrices must be nonempty two-dimensional arrays")
    track_count = int(first_matrix.shape[0])
    del first_matrix
    coordinate_arrays: dict[str, np.ndarray] = {}
    for chromosome in sorted(matrix_keys):
        matrix = np.asarray(chromosome_matrices[chromosome], dtype=np.float64)
        if matrix.ndim != 2 or matrix.shape[0] != track_count:
            raise ValueError("Every chromosome matrix must have the same track count")
        if matrix.shape[1] == 0:
            raise ValueError(f"{chromosome} matrix must contain bins")
        raw_coordinates = np.asarray(chromosome_coordinates[chromosome])
        if raw_coordinates.ndim != 1 or raw_coordinates.size != matrix.shape[1]:
            raise ValueError(f"{chromosome} coordinates must align with matrix columns")
        if np.issubdtype(raw_coordinates.dtype, np.integer):
            if np.issubdtype(raw_coordinates.dtype, np.unsignedinteger) and np.any(
                raw_coordinates > np.iinfo(np.int64).max
            ):
                raise ValueError(f"{chromosome} coordinates exceed int64 range")
            coordinates = np.asarray(raw_coordinates, dtype=np.int64)
        elif np.issubdtype(raw_coordinates.dtype, np.floating):
            numeric_coordinates = np.asarray(raw_coordinates, dtype=np.float64)
            if np.any(~np.isfinite(numeric_coordinates)):
                raise ValueError(f"{chromosome} coordinates must be finite")
            coordinates = numeric_coordinates.astype(np.int64)
            if not np.array_equal(
                numeric_coordinates,
                coordinates.astype(np.float64),
            ):
                raise ValueError(f"{chromosome} coordinates must be integers")
        else:
            raise TypeError(f"{chromosome} coordinates must be numeric")
        coordinate_differences = np.diff(coordinates)
        if np.any(coordinates < 0) or np.any(coordinate_differences <= 0):
            raise ValueError(
                f"{chromosome} coordinates must be nonnegative and increasing"
            )
        if np.any(coordinate_differences != step):
            raise ValueError(f"{chromosome} coordinates must use `step_bp` spacing")
        coordinate_arrays[chromosome] = coordinates
        del matrix

    max_lag_bins = min(_MAX_LAG_BP // step, (window_size // step) - 1)
    if max_lag_bins < _ACF_REQUIRED_CROSSINGS:
        raise ValueError("`window_bp` contains too few bins for the ACF crossing rule")
    diagnostics: dict[str, Any] = {
        "method": "finitePairWindowACF",
        "usedPrior": False,
        "status": "unavailable",
        "reason": "none",
        "stepBP": step,
        "windowBP": window_size,
        "windowCountRequested": requested_windows,
        "windowCountSelected": 0,
        "workingQuantile": working_q,
        "minimumCorrelationRadiusBP": MIN_CORRELATION_RADIUS_BP,
        "bootstrapDraws": draws,
        "bootstrapSeed": seed,
        "trackCount": track_count,
        "autosomeCount": int(sum(_is_autosome(chrom) for chrom in matrix_keys)),
        "selectedCoordinates": [],
        "windowEstimates": [],
    }
    candidates, eligible_window_count = _candidate_windows(
        chromosome_matrices,
        coordinate_arrays,
        step,
        window_size,
        requested_windows,
    )
    diagnostics["eligibleWindowCount"] = eligible_window_count
    if eligible_window_count < _MIN_WINDOWS:
        reason = f"Dependence estimation needs at least {_MIN_WINDOWS} eligible windows"
        return _prior_result(reason, step, diagnostics)

    selected_slots: list[dict[str, Any] | None] = [None] * len(candidates)
    candidate_indices_by_chromosome = {
        chromosome: [
            candidate_index
            for candidate_index, candidate in enumerate(candidates)
            if candidate[1] == chromosome
        ]
        for chromosome in sorted({candidate[1] for candidate in candidates})
    }
    for chromosome, candidate_indices in candidate_indices_by_chromosome.items():
        source_matrix = chromosome_matrices[chromosome]
        matrix = np.asarray(source_matrix, dtype=np.float64)
        for candidate_index in candidate_indices:
            negative_score, _, start_bp, start_idx, stop_idx = candidates[
                candidate_index
            ]
            track_window = matrix[:, start_idx:stop_idx]
            track_acfs, track_pair_counts = _prepare_track_acfs(
                track_window,
                max_lag_bins,
                _CLIP_QUANTILES,
                max(20, 2 * _ACF_REQUIRED_CROSSINGS),
            )
            pooled_acf = _median_across_tracks(track_acfs)
            radii_bp, censored, valid, smoothed = _crossings_from_pooled_acfs(
                pooled_acf,
                step,
                _ACF_THRESHOLD,
                _ACF_REQUIRED_CROSSINGS,
                True,
            )
            if not bool(valid[0]) or smoothed is None:
                continue
            result = _window_result(
                track_acfs,
                track_pair_counts,
                pooled_acf,
                smoothed[0],
                float(radii_bp[0]),
                bool(censored[0]),
                step,
                _ACF_THRESHOLD,
                _ACF_REQUIRED_CROSSINGS,
            )
            selected_slots[candidate_index] = {
                "chromosome": chromosome,
                "startBP": start_bp,
                "endBP": start_bp + window_size,
                "rankingScore": -negative_score,
                "trackACFs": track_acfs,
                "trackPairCounts": track_pair_counts,
                "result": result,
            }
        del matrix
        del source_matrix
    selected = [record for record in selected_slots if record is not None]
    selected_chromosomes = sorted({record["chromosome"] for record in selected})
    diagnostics["windowCountSelected"] = len(selected)
    diagnostics["selectedAutosomeCount"] = len(selected_chromosomes)
    if len(selected) < _MIN_WINDOWS or len(selected_chromosomes) < _MIN_AUTOSOMES:
        reason = (
            f"Dependence estimation needs at least {_MIN_WINDOWS} usable windows "
            f"from at least {_MIN_AUTOSOMES} autosomes"
        )
        return _prior_result(reason, step, diagnostics)

    selected_coordinates = []
    window_estimates = []
    for record in selected:
        result = record["result"]
        selected_coordinates.append(
            {
                "chromosome": record["chromosome"],
                "startBP": int(record["startBP"]),
                "endBP": int(record["endBP"]),
                "rankingScore": float(record["rankingScore"]),
            }
        )
        window_estimates.append(
            {
                "chromosome": record["chromosome"],
                "startBP": int(record["startBP"]),
                "endBP": int(record["endBP"]),
                "rawCrossingLagBP": result["rawCrossingLagBP"],
                "gaussianEquivalentRadiusBP": float(
                    result["gaussianEquivalentRadiusBP"]
                ),
                "rightCensored": bool(result["rightCensored"]),
                "dominantACFPeriodBP": result["dominantACFPeriodBP"],
                "oscillationStrength": float(result["oscillationStrength"]),
                "postCrossingACFRevival": float(result["postCrossingACFRevival"]),
                "finitePairCountMin": float(result["finitePairCountMin"]),
                "finitePairCountMedian": float(result["finitePairCountMedian"]),
                "finitePairCountMax": float(result["finitePairCountMax"]),
            }
        )
    diagnostics["selectedCoordinates"] = selected_coordinates
    diagnostics["windowEstimates"] = window_estimates

    point_durations = np.asarray(
        [record["result"]["gaussianEquivalentRadiusBP"] for record in selected],
        dtype=np.float64,
    )
    point_censored = np.asarray(
        [record["result"]["rightCensored"] for record in selected], dtype=bool
    )
    point_median_bp = _kaplan_meier_quantile(point_durations, point_censored, 0.5)
    point_working_bp = _kaplan_meier_quantile(
        point_durations, point_censored, working_q
    )
    if point_median_bp is None or point_working_bp is None:
        return _prior_result(
            "Window radius median or working quantile is unidentified",
            step,
            diagnostics,
        )

    acf_cube = np.stack([record["trackACFs"] for record in selected], axis=0)
    indices_by_chromosome = {
        chromosome: np.asarray(
            [
                idx
                for idx, record in enumerate(selected)
                if record["chromosome"] == chromosome
            ],
            dtype=np.int64,
        )
        for chromosome in selected_chromosomes
    }
    rng = np.random.default_rng(seed)
    bootstrap_estimate_values: list[float] = []
    bootstrap_working_values: list[float] = []
    chromosome_count = len(selected_chromosomes)
    for draw_idx in range(draws):
        track_sample = rng.integers(0, track_count, size=track_count)
        pooled_draw = _median_across_tracks(acf_cube[:, track_sample, :])
        draw_radii, draw_censored, draw_valid, _ = _crossings_from_pooled_acfs(
            pooled_draw,
            step,
            _ACF_THRESHOLD,
            _ACF_REQUIRED_CROSSINGS,
            False,
        )
        sampled_chromosome_indices = rng.integers(
            0, chromosome_count, size=chromosome_count
        )
        sampled_windows = []
        for chromosome_idx in sampled_chromosome_indices:
            chromosome = selected_chromosomes[int(chromosome_idx)]
            chromosome_windows = indices_by_chromosome[chromosome]
            sampled_windows.append(
                rng.choice(
                    chromosome_windows,
                    size=chromosome_windows.size,
                    replace=True,
                )
            )
        draw_window_indices = np.concatenate(sampled_windows)
        draw_window_indices = draw_window_indices[draw_valid[draw_window_indices]]
        if draw_window_indices.size == 0:
            continue
        draw_estimate = _kaplan_meier_quantile(
            draw_radii[draw_window_indices],
            draw_censored[draw_window_indices],
            0.5,
        )
        draw_working = _kaplan_meier_quantile(
            draw_radii[draw_window_indices],
            draw_censored[draw_window_indices],
            working_q,
        )
        if draw_estimate is None or draw_working is None:
            continue
        bootstrap_estimate_values.append(float(draw_estimate))
        bootstrap_working_values.append(float(draw_working))

    minimum_valid_draws = max(20, int(math.ceil(0.8 * draws)))
    if len(bootstrap_estimate_values) < minimum_valid_draws:
        return _prior_result(
            "Too few bootstrap draws have identified radius quantiles",
            step,
            diagnostics,
        )
    unconstrained_bootstrap_estimates = np.asarray(
        bootstrap_estimate_values,
        dtype=np.float64,
    )
    unconstrained_bootstrap_working_spans = np.asarray(
        bootstrap_working_values,
        dtype=np.float64,
    )
    unconstrained_estimate_bp = float(np.median(unconstrained_bootstrap_estimates))
    unconstrained_working_span_bp = float(
        np.median(unconstrained_bootstrap_working_spans)
    )
    bootstrap_estimates = np.maximum(
        unconstrained_bootstrap_estimates,
        MIN_CORRELATION_RADIUS_BP,
    )
    bootstrap_working_spans = np.maximum(
        unconstrained_bootstrap_working_spans,
        MIN_CORRELATION_RADIUS_BP,
    )

    estimate_bp = float(np.median(bootstrap_estimates))
    lower_bp, upper_bp = (
        float(value) for value in np.quantile(bootstrap_estimates, (0.025, 0.975))
    )
    working_span_bp = float(np.median(bootstrap_working_spans))
    estimate_intervals = int(math.ceil(estimate_bp / float(step)))
    lower_intervals = int(math.ceil(lower_bp / float(step)))
    upper_intervals = int(math.ceil(upper_bp / float(step)))
    working_span_intervals = int(math.ceil(working_span_bp / float(step)))

    periods = np.asarray(
        [
            value
            for value in (
                estimate["dominantACFPeriodBP"] for estimate in window_estimates
            )
            if value is not None
        ],
        dtype=np.float64,
    )
    oscillation_strengths = np.asarray(
        [estimate["oscillationStrength"] for estimate in window_estimates],
        dtype=np.float64,
    )
    revivals = np.asarray(
        [estimate["postCrossingACFRevival"] for estimate in window_estimates],
        dtype=np.float64,
    )
    finite_pair_medians = np.asarray(
        [estimate["finitePairCountMedian"] for estimate in window_estimates],
        dtype=np.float64,
    )
    diagnostics.update(
        {
            "status": "estimated",
            "reason": "none",
            "estimateBP": estimate_bp,
            "lowerBP": lower_bp,
            "upperBP": upper_bp,
            "estimateIntervals": estimate_intervals,
            "lowerIntervals": lower_intervals,
            "upperIntervals": upper_intervals,
            "workingSpanBP": working_span_bp,
            "workingSpanIntervals": working_span_intervals,
            "unconstrainedEstimateBP": unconstrained_estimate_bp,
            "unconstrainedWorkingSpanBP": unconstrained_working_span_bp,
            "minimumCorrelationRadiusApplied": bool(
                unconstrained_estimate_bp < MIN_CORRELATION_RADIUS_BP
                or unconstrained_working_span_bp < MIN_CORRELATION_RADIUS_BP
            ),
            "pointMedianRadiusBP": point_median_bp,
            "pointWorkingSpanBP": point_working_bp,
            "workingSpanLowerBP": float(np.quantile(bootstrap_working_spans, 0.025)),
            "workingSpanUpperBP": float(np.quantile(bootstrap_working_spans, 0.975)),
            "censorFraction": float(np.mean(point_censored)),
            "radiusDistributionBP": point_durations.tolist(),
            "radiusEventObserved": (~point_censored).tolist(),
            "bootstrapEstimateBP": bootstrap_estimates.tolist(),
            "bootstrapWorkingSpanBP": bootstrap_working_spans.tolist(),
            "bootstrapDrawsValid": int(bootstrap_estimates.size),
            "dominantACFPeriodBP": (
                float(np.median(periods)) if periods.size else None
            ),
            "oscillationStrength": float(np.median(oscillation_strengths)),
            "postCrossingACFRevival": float(np.median(revivals)),
            "finitePairCountMedian": float(np.median(finite_pair_medians)),
        }
    )
    return estimate_intervals, lower_intervals, upper_intervals, diagnostics
