r"""
ROCCO: Inference
==================================================================================

Functions for scoring loci and setting chromosome-specific budget and gamma values.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import sys
import warnings
from typing import Any, Dict, Tuple

import numpy as np
from scipy import linalg, ndimage, optimize, signal, special, stats

try:
    from . import _baseline as _baseline_native
except ImportError:
    _baseline_native = None

try:
    from . import _wls as _wls_native
except ImportError:
    _wls_native = None

logger = logging.getLogger(__name__)
_BUDGET_NULL_PROCESS_STATE: dict[str, Any] | None = None


def _robust_scale(values: np.ndarray, floor: float = 1.0e-6) -> float:
    values_ = np.asarray(values, dtype=np.float64)
    if values_.size == 0:
        return float(floor)
    mad = np.median(np.abs(values_ - np.median(values_)))
    return float(max(mad * 1.4826, floor))


def _log_scale_wls_matrix(
    chrom_matrix: np.ndarray,
    pseudocount: float = 1.0,
) -> np.ndarray:
    matrix = np.asarray(chrom_matrix, dtype=np.float64)
    if np.any(~np.isfinite(matrix)):
        raise ValueError("`chrom_matrix` contains non-finite values")
    return np.log2(np.clip(matrix, 0.0, None) + float(pseudocount))


def _resolve_local_baseline_window(
    n_loci: int,
    target_window: int = 101,
) -> int:
    n_loci = int(n_loci)
    if n_loci < 25:
        return 0
    window = int(max(3, target_window))
    if window > n_loci:
        window = n_loci
    if (window % 2) == 0:
        window = window - 1 if window == n_loci else window + 1
    return int(max(0, window))


def _consenrich_whittaker_lambda(block_size: int) -> float:
    r"""Map a smoothing block size to the Whittaker penalty used in Consenrich.

    This follows the same :math:`\texttt{blockSize} \mapsto \lambda` rule used
    by ``clocalBaseline`` in the Consenrich repo:
    https://github.com/nolan-h-hamilton/Consenrich
    """
    block = int(max(3, block_size))
    if (block % 2) == 0:
        block += 1
    w_hat = float(block) * 0.15915494
    return float(7.0 * (w_hat**4))


def _consenrich_whittaker_penalty_bands(
    n_loci: int,
    penalty_lambda: float,
) -> np.ndarray:
    bands = np.zeros((3, int(n_loci)), dtype=np.float64)
    diag = np.full(int(n_loci), 6.0, dtype=np.float64)
    diag[0] = 1.0
    diag[-1] = 1.0
    if n_loci > 1:
        diag[1] = 5.0
        diag[-2] = 5.0
        off1 = np.full(int(n_loci) - 1, -4.0, dtype=np.float64)
        off1[0] = -2.0
        off1[-1] = -2.0
        bands[1, 1:] = float(penalty_lambda) * off1
    if n_loci > 2:
        bands[0, 2:] = float(penalty_lambda)
    bands[2, :] = float(penalty_lambda) * diag
    return bands


def _consenrich_masked_whittaker_baseline(
    y_vals: np.ndarray,
    fit_mask: np.ndarray,
    penalty_bands: np.ndarray,
) -> np.ndarray:
    r"""Solve the masked Whittaker baseline system used by Consenrich.

    This is a small pure-Python/SciPy port of the masked solve
    :math:`(W + \lambda D^\top D)b = Wy` behind
    ``locBaselineMasked_F64`` in the Consenrich repo:
    https://github.com/nolan-h-hamilton/Consenrich
    """
    y_arr = np.asarray(y_vals, dtype=np.float64)
    mask_arr = np.clip(np.asarray(fit_mask, dtype=np.float64), 0.0, None)
    if y_arr.ndim != 1:
        raise ValueError("`y_vals` must be one-dimensional")
    if mask_arr.shape != y_arr.shape:
        raise ValueError("`fit_mask` must match `y_vals`")
    if y_arr.size < 3:
        return y_arr.copy()
    # Add the mask to the main diagonal so the system is
    # :math:`(W + \lambda D^\top D)b = Wy`.
    bands = np.array(penalty_bands, dtype=np.float64, copy=True)
    bands[2, :] += mask_arr
    rhs = mask_arr * y_arr
    return linalg.solveh_banded(
        bands,
        rhs,
        lower=False,
        check_finite=False,
        overwrite_ab=True,
        overwrite_b=True,
    )


def _consenrich_crossfit_whittaker_baseline(
    y_vals: np.ndarray,
    block_size: int = 101,
) -> np.ndarray:
    r"""Estimate a broad local baseline with a Consenrich-style cross-fit smoother.

    This follows the even/odd masked Whittaker trick used by
    ``locBaselineCrossfit2_F64`` in the Consenrich repo, so the final estimate
    is :math:`(\hat b_{\mathrm{even}} + \hat b_{\mathrm{odd}}) / 2`:
    https://github.com/nolan-h-hamilton/Consenrich

    Small inputs return zeros on purpose. That matches the old guard in
    ``clocalBaseline`` and keeps short chromosomes from getting weird.
    """
    y_arr = np.asarray(y_vals, dtype=np.float64)
    if y_arr.ndim != 1:
        raise ValueError("`y_vals` must be one-dimensional")
    n_loci = int(y_arr.size)
    window = _resolve_local_baseline_window(n_loci, target_window=block_size)
    if window == 0:
        return np.zeros_like(y_arr, dtype=np.float64)

    penalty_lambda = _consenrich_whittaker_lambda(window)
    if _baseline_native is not None:
        return np.asarray(
            _baseline_native.crossfit_whittaker_baseline(
                y_arr,
                penalty_lambda=penalty_lambda,
            ),
            dtype=np.float64,
        )
    penalty_bands = _consenrich_whittaker_penalty_bands(
        n_loci,
        penalty_lambda,
    )
    even_mask = ((np.arange(n_loci) % 2) == 0).astype(np.float64)
    odd_mask = 1.0 - even_mask
    even_baseline = _consenrich_masked_whittaker_baseline(
        y_arr,
        even_mask,
        penalty_bands,
    )
    odd_baseline = _consenrich_masked_whittaker_baseline(
        y_arr,
        odd_mask,
        penalty_bands,
    )
    return 0.5 * (even_baseline + odd_baseline)


def _estimate_local_background_matrix(
    centered_matrix: np.ndarray,
    target_window: int = 101,
) -> tuple[np.ndarray, int, float]:
    matrix = np.asarray(centered_matrix, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("`centered_matrix` must be two-dimensional")
    n_samples, n_loci = matrix.shape
    window = _resolve_local_baseline_window(n_loci, target_window=target_window)
    if window == 0:
        return np.zeros_like(matrix, dtype=np.float64), 0, 0.0

    penalty_lambda = _consenrich_whittaker_lambda(window)
    if _baseline_native is not None:
        # Apply the same cross-fit baseline to each sample track.
        local_baselines = np.asarray(
            _baseline_native.crossfit_whittaker_baseline(
                matrix,
                penalty_lambda=penalty_lambda,
            ),
            dtype=np.float64,
        )
        if not np.all(np.isfinite(local_baselines)):
            raise ValueError("Local baseline fit produced non-finite values")
        return local_baselines, window, penalty_lambda
    penalty_bands = _consenrich_whittaker_penalty_bands(n_loci, penalty_lambda)
    local_baselines = np.empty_like(matrix, dtype=np.float64)
    even_mask = ((np.arange(n_loci) % 2) == 0).astype(np.float64)
    odd_mask = 1.0 - even_mask
    for sample_idx in range(n_samples):
        even_baseline = _consenrich_masked_whittaker_baseline(
            matrix[sample_idx],
            even_mask,
            penalty_bands,
        )
        odd_baseline = _consenrich_masked_whittaker_baseline(
            matrix[sample_idx],
            odd_mask,
            penalty_bands,
        )
        local_baselines[sample_idx] = 0.5 * (even_baseline + odd_baseline)
    if not np.all(np.isfinite(local_baselines)):
        raise ValueError("Local baseline fit produced non-finite values")
    return local_baselines, window, penalty_lambda


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
    if _wls_native is None:
        raise RuntimeError("Make sure native C extensions are built and available")
    precision_floor_ratio_ = float(max(precision_floor_ratio, 0.0))
    # run EB munc + WLS to get per-locus scores and details for the constrained optimization (DP)
    result_native = _wls_native.score_centered_wls(
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
    result = (
        np.asarray(scores_arr, dtype=np.float64),
        {
            "mean": mean,
            "raw_variance": np.asarray(raw_var_arr, dtype=np.float64),
            "prior_variance": np.asarray(prior_var_arr, dtype=np.float64),
            "moderated_variance": np.asarray(moderated_var_arr, dtype=np.float64),
            "standard_error": se,
            "z_scores": mean / np.maximum(se, 1.0e-8),
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
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    low_memory: bool = False,
    return_details: bool = False,
) -> np.ndarray | Tuple[np.ndarray, Dict[str, Any]]:
    r"""Score loci with an EB-moderated summary on baseline-corrected log signal.

    We model the log-scaled signal at each locus as
    :math:`y_{ij} = \log(x_{ij}) = b_{ij} + \mu_j + e_{ij}` where
    :math:`b_{ij}` is a broad local background term from a cross-fit Whittaker
    smoother (see Consenrich api). On the resulting ``m x n`` centered matrix, ROCCO follows the
    Consenrich pattern track by track: each sample gets a local variance track
    from rolling AR(1) innovation variances, a separate global monotone prior
    variance trend as a function of absolute signal level, and an EB-shrunk
    posterior variance track based on the prior variance. Final locus estimate and standard error are
    then computed by WLS over the full data matrix and the full posterior
    variance matrix. The default score for the constrained optimization is the 'moderated'
    standardized effect :math:`t_j = \mu_j / \mathrm{se}_j`.
    """
    matrix = _log_scale_wls_matrix(chrom_matrix)
    if matrix.ndim != 2:
        raise ValueError("`chrom_matrix` must be two-dimensional")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("`chrom_matrix` must be non-empty")

    # Use a robust pilot offset before baseline fitting so the smoother targets
    # :math:`b_{ij}` instead of burning effort on a sample-level shift.
    baseline_init = np.median(matrix, axis=1, keepdims=True)
    global_centered = matrix - baseline_init
    local_baselines, local_window, local_lambda = _estimate_local_background_matrix(
        global_centered
    )
    centered = global_centered - local_baselines
    del matrix
    del global_centered
    del local_baselines
    scores, core_details = _score_centered_wls_matrix(
        centered,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
    )
    if not np.all(np.isfinite(scores)):
        raise ValueError("Locus scoring produced non-finite values")
    centered_out = centered.astype(
        np.float32 if low_memory else np.float64,
        copy=False,
    )
    del centered

    details = {
        "input_scale": "log2p1",
        "local_baseline_window": int(local_window),
        "local_baseline_lambda": float(local_lambda),
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
        "centered_matrix": centered_out,
    }
    if return_details:
        return scores.astype(np.float64), details
    return scores.astype(np.float64)


def benjamini_hochberg(
    p_values: np.ndarray,
    fdr: float = 0.01,
) -> np.ndarray:
    p_values_ = np.asarray(p_values, dtype=np.float64)
    if p_values_.ndim != 1:
        raise ValueError("`p_values` must be one-dimensional")
    m = p_values_.shape[0]
    if m == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(p_values_)
    ranked = p_values_[order]
    thresholds = float(fdr) * (np.arange(1, m + 1) / float(m))
    passing = ranked <= thresholds
    mask = np.zeros(m, dtype=bool)
    if np.any(passing):
        cutoff = np.max(np.where(passing)[0])
        mask[order[: cutoff + 1]] = True
    return mask


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
    # Geyer's IPS keeps :math:`\tau_{\mathrm{int}}` on the valid side by
    # summing adjacent autocorrelation pairs until the first nonpositive pair.
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


def _resolve_budget_ess_max_lag(
    n_loci: int,
    dependence_lag_hint: int | None = None,
) -> int:
    r"""Resolve the ESS autocorrelation lag cap from a broad background scale."""
    n_loci_ = int(max(1, n_loci))
    if dependence_lag_hint is None:
        return int(min(n_loci_ - 1, max(16, 4 * min(n_loci_, 101))))
    return int(
        min(
            n_loci_ - 1,
            max(16, 4 * max(1, min(n_loci_, int(dependence_lag_hint)))),
        )
    )


def _resolve_budget_bootstrap_bandwidth(
    n_loci: int,
    dependence_lag_hint: int | None = None,
) -> int:
    r"""Resolve the dependent-multiplier bandwidth for the budget null (DWB)."""
    n_loci_ = int(max(1, n_loci))
    if n_loci_ <= 1:
        return 1
    if dependence_lag_hint is None:
        return int(min(n_loci_ - 1, max(8, round(n_loci_ ** (1.0 / 3.0)))))
    return int(min(n_loci_ - 1, max(8, int(dependence_lag_hint))))


def _build_budget_bootstrap_kernel(
    bandwidth: int,
) -> np.ndarray:
    r"""Build the Bartlett kernel used by the dependent wild bootstrap to induce short-range dependence."""
    bandwidth_ = int(max(1, bandwidth))
    support = np.arange(-bandwidth_, bandwidth_ + 1, dtype=np.float64)
    kernel = np.maximum(1.0 - (np.abs(support) / float(bandwidth_ + 1)), 0.0)
    kernel /= np.sqrt(np.sum(kernel * kernel))
    return kernel.astype(np.float64, copy=False)


def _generate_dependent_wild_weights(
    n_loci: int,
    kernel: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    r"""Draw a short-range dependent multiplier process for the wild bootstrap.

    Independent innovations are smoothed by a precomputed *Bartlett kernel* so
    the multiplier field has mean zero, unit variance, and a fixed dependence
    scale across all draws (i.e., the dependent wild bootstrap).
    """
    n_loci_ = int(max(1, n_loci))
    if n_loci_ == 1:
        return np.ones(1, dtype=np.float64)

    kernel_ = np.asarray(kernel, dtype=np.float64)
    innovations = rng.standard_normal(n_loci_ + kernel_.size - 1)
    weights = signal.fftconvolve(innovations, kernel_, mode="valid")
    weights = np.asarray(weights, dtype=np.float64)
    weights -= float(np.mean(weights))
    weight_scale = float(np.std(weights))
    if not np.isfinite(weight_scale) or weight_scale <= 1.0e-8:
        fallback = rng.choice(np.array([-1.0, 1.0], dtype=np.float64), size=n_loci_)
        fallback -= float(np.mean(fallback))
        weight_scale = float(np.std(fallback))
        return fallback / max(weight_scale, 1.0e-6)
    return weights / weight_scale


def _update_running_moments(
    count: int,
    mean: float,
    m2: float,
    new_value: float,
) -> tuple[int, float, float]:
    # Welford's algorithm for numerically stable online mean and variance updates
    count_ = int(count) + 1
    delta = float(new_value) - float(mean)
    mean_ = float(mean) + (delta / float(count_))
    delta2 = float(new_value) - mean_
    m2_ = float(m2) + (delta * delta2)
    return count_, mean_, m2_


def _budget_null_stable_enough(
    count: int,
    mean: float,
    m2: float,
    min_draws: int,
    abs_tol: float,
    rel_tol: float,
) -> bool:
    if int(count) < int(max(2, min_draws)):
        return False
    sample_var = float(max(m2 / float(max(count - 1, 1)), 0.0))
    stderr = float(np.sqrt(sample_var / float(max(count, 1))))
    target = float(max(abs_tol, rel_tol * max(abs(mean), 1.0e-6)))
    return bool(stderr <= target)


def _choose_budget_pool_context() -> mp.context.BaseContext:
    start_methods = mp.get_all_start_methods()
    if "fork" in start_methods:
        return mp.get_context("fork")
    return mp.get_context(start_methods[0])


def _init_budget_null_process(
    residual_template: np.ndarray,
    lower_bound_z: float,
    prior_df: float,
    min_effect: float | None,
    precision_floor_ratio: float,
    null_center: float,
    null_soft_scale: float,
    null_threshold: float,
    kernel: np.ndarray,
    base_seed: int,
) -> None:
    global _BUDGET_NULL_PROCESS_STATE
    _BUDGET_NULL_PROCESS_STATE = {
        "residual_template": np.asarray(residual_template, dtype=np.float64),
        "lower_bound_z": float(lower_bound_z),
        "prior_df": float(prior_df),
        "min_effect": None if min_effect is None else float(max(min_effect, 0.0)),
        "precision_floor_ratio": float(max(precision_floor_ratio, 0.0)),
        "null_center": float(null_center),
        "null_soft_scale": float(null_soft_scale),
        "null_threshold": float(null_threshold),
        "kernel": np.asarray(kernel, dtype=np.float64),
        "base_seed": int(base_seed),
    }


def _compute_budget_null_draw(draw_index: int) -> tuple[float, float, float, float]:
    state = _BUDGET_NULL_PROCESS_STATE
    if state is None:
        raise RuntimeError("Budget null process state is not initialized")

    residual_template = np.asarray(state["residual_template"], dtype=np.float64)
    lower_bound_z = float(state["lower_bound_z"])
    prior_df = float(state["prior_df"])
    min_effect = state["min_effect"]
    precision_floor_ratio = float(state["precision_floor_ratio"])
    null_center = float(state["null_center"])
    null_soft_scale = float(state["null_soft_scale"])
    null_threshold = float(state["null_threshold"])
    kernel = np.asarray(state["kernel"], dtype=np.float64)
    base_seed = int(state["base_seed"])
    n_samples, n_loci = residual_template.shape
    rng = np.random.default_rng(base_seed + (104729 * (int(draw_index) + 1)))
    bootstrap_centered = np.empty_like(residual_template, dtype=np.float64)

    for sample_idx in range(n_samples):
        wild_weights = _generate_dependent_wild_weights(
            n_loci,
            kernel=kernel,
            rng=rng,
        )
        bootstrap_centered[sample_idx] = residual_template[sample_idx] * wild_weights

    bootstrap_scores, _ = _score_centered_wls_matrix(
        bootstrap_centered,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
    )
    bootstrap_residual_scores = (
        np.asarray(bootstrap_scores, dtype=np.float64) - null_center
    )
    bootstrap_positive = np.clip(bootstrap_residual_scores, 0.0, None)
    return (
        float(np.mean(bootstrap_positive)),
        float(np.mean(bootstrap_positive / null_soft_scale)),
        float(np.mean(bootstrap_positive > 0.0)),
        float(np.mean(bootstrap_scores > null_threshold)),
    )


def _fit_budget_null_residual_template(
    centered_matrix: np.ndarray,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""Fit the one-sided null residual template used by the budget bootstrap.

    ROCCO scores the centered matrix as :math:`y_{ij} = \mu_j + e_{ij}` after
    removing :math:`a_i + b_{ij}`. The budget null sets
    :math:`\mu_j^{(0)} = \max(\hat\mu_j, 0)` to zero out positive consensus
    signal, then uses

    .. math::

       \tilde e_{ij} = y_{ij} - \max(\hat\mu_j, 0)

    as the residual template for the bootstrap.
    """
    observed_scores, score_details = _score_centered_wls_matrix(
        centered_matrix,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
    )
    mu_hat = np.asarray(score_details["mean"], dtype=np.float64)
    positive_consensus = np.clip(mu_hat, 0.0, None)
    residual_template = (
        np.asarray(centered_matrix, dtype=np.float64) - positive_consensus[None, :]
    )
    return residual_template, observed_scores.astype(np.float64), positive_consensus


def _estimate_wild_bootstrap_score_null(
    centered_matrix: np.ndarray,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    observed_scores: np.ndarray | None = None,
    dependence_lag_hint: int | None = None,
    num_null_draws: int = 25,
    random_seed: int = 0,
    progress_label: str | None = None,
    num_processes: int = 1,
    min_null_draws: int | None = None,
    stability_abs_tol: float = 5.0e-3,
    stability_rel_tol: float = 5.0e-2,
) -> dict[str, float | int | str | np.ndarray]:
    r"""Estimate the chromosome score null by a dependent wild residual bootstrap.

    The fitted null starts from the residual template
    :math:`\tilde e_{ij}`. For each bootstrap draw, ROCCO generates a
    short-range dependent multiplier field :math:`W_{ij}` and forms

    .. math::

       y_{ij}^{*} = \tilde e_{ij} W_{ij}.

    The WLS score is recomputed on :math:`y_{ij}^{*}`. The null center and null
    scale come from the score field of the fitted null template, while the
    bootstrap draws estimate the positive-tail score mass expected under that
    fitted null.
    """
    centered = np.asarray(centered_matrix, dtype=np.float64)
    residual_template, fitted_scores, positive_consensus = (
        _fit_budget_null_residual_template(
            centered,
            lower_bound_z=lower_bound_z,
            prior_df=prior_df,
            min_effect=min_effect,
            precision_floor_ratio=precision_floor_ratio,
        )
    )
    if observed_scores is None:
        observed_scores_ = fitted_scores
    else:
        observed_scores_ = np.asarray(observed_scores, dtype=np.float64)
        if observed_scores_.shape[0] != centered.shape[1]:
            raise ValueError(
                "`observed_scores` must have the same number of loci as `centered_matrix`"
            )

    null_reference_scores, _ = _score_centered_wls_matrix(
        residual_template,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
    )
    finite_null_scores = np.asarray(null_reference_scores, dtype=np.float64)
    null_center = float(np.median(finite_null_scores))
    null_reference_residuals = finite_null_scores - null_center
    negative_reference = null_reference_residuals[null_reference_residuals <= 0.0]
    if negative_reference.size == 0:
        negative_magnitudes = np.abs(null_reference_residuals)
    else:
        negative_magnitudes = -negative_reference
    if negative_magnitudes.size == 0:
        negative_magnitudes = np.array([0.0], dtype=np.float64)
    mirrored_reference = np.concatenate((-negative_magnitudes, negative_magnitudes))
    null_scale = float(_robust_scale(mirrored_reference))
    if not np.isfinite(null_center) or not np.isfinite(null_scale):
        raise ValueError("Budget null fit produced non-finite values")
    null_soft_scale = float(max(null_scale, 1.0e-6))
    # Conservative score threshold at :math:`\mu_0 + 2\sigma_0`, where sigma_0 is null scale
    null_threshold = float(null_center + (2.0 * null_scale))

    n_samples, n_loci = centered.shape
    bandwidth = _resolve_budget_bootstrap_bandwidth(
        n_loci,
        dependence_lag_hint=dependence_lag_hint,
    )
    kernel = _build_budget_bootstrap_kernel(bandwidth)
    num_draws = int(max(1, num_null_draws))
    min_draws = int(
        min(num_draws, max(4, 8 if min_null_draws is None else min_null_draws))
    )
    process_count = int(min(max(1, num_processes), num_draws))
    batch_size = int(max(1, process_count))
    draws_used = 0
    mean_mass = 0.0
    m2_mass = 0.0
    mean_units = 0.0
    m2_units = 0.0
    mean_fraction = 0.0
    m2_fraction = 0.0
    mean_tail_occupancy = 0.0
    m2_tail_occupancy = 0.0

    if process_count <= 1:
        _init_budget_null_process(
            residual_template,
            lower_bound_z,
            prior_df,
            min_effect,
            precision_floor_ratio,
            null_center,
            null_soft_scale,
            null_threshold,
            kernel,
            int(random_seed),
        )
        for batch_start in range(0, num_draws, batch_size):
            draw_ids = range(batch_start, min(num_draws, batch_start + batch_size))
            for draw_id in draw_ids:
                null_mass, null_units, null_fraction, null_tail_occupancy = (
                    _compute_budget_null_draw(draw_id)
                )
                draws_used, mean_mass, m2_mass = _update_running_moments(
                    draws_used,
                    mean_mass,
                    m2_mass,
                    null_mass,
                )
                _, mean_units, m2_units = _update_running_moments(
                    draws_used - 1,
                    mean_units,
                    m2_units,
                    null_units,
                )
                _, mean_fraction, m2_fraction = _update_running_moments(
                    draws_used - 1,
                    mean_fraction,
                    m2_fraction,
                    null_fraction,
                )
                _, mean_tail_occupancy, m2_tail_occupancy = _update_running_moments(
                    draws_used - 1,
                    mean_tail_occupancy,
                    m2_tail_occupancy,
                    null_tail_occupancy,
                )
                if progress_label:
                    sys.stderr.write(f"\r{progress_label}: {draws_used}/{num_draws}")
                    sys.stderr.flush()
            if _budget_null_stable_enough(
                draws_used,
                mean_units,
                m2_units,
                min_draws=min_draws,
                abs_tol=stability_abs_tol,
                rel_tol=stability_rel_tol,
            ):
                break
    else:
        context = _choose_budget_pool_context()
        with context.Pool(
            processes=process_count,
            initializer=_init_budget_null_process,
            initargs=(
                residual_template,
                lower_bound_z,
                prior_df,
                min_effect,
                precision_floor_ratio,
                null_center,
                null_soft_scale,
                null_threshold,
                kernel,
                int(random_seed),
            ),
        ) as pool:
            for batch_start in range(0, num_draws, batch_size):
                draw_ids = list(
                    range(batch_start, min(num_draws, batch_start + batch_size))
                )
                batch_results = pool.map(_compute_budget_null_draw, draw_ids)
                for (
                    null_mass,
                    null_units,
                    null_fraction,
                    null_tail_occupancy,
                ) in batch_results:
                    draws_used, mean_mass, m2_mass = _update_running_moments(
                        draws_used,
                        mean_mass,
                        m2_mass,
                        null_mass,
                    )
                    _, mean_units, m2_units = _update_running_moments(
                        draws_used - 1,
                        mean_units,
                        m2_units,
                        null_units,
                    )
                    _, mean_fraction, m2_fraction = _update_running_moments(
                        draws_used - 1,
                        mean_fraction,
                        m2_fraction,
                        null_fraction,
                    )
                    _, mean_tail_occupancy, m2_tail_occupancy = _update_running_moments(
                        draws_used - 1,
                        mean_tail_occupancy,
                        m2_tail_occupancy,
                        null_tail_occupancy,
                    )
                    if progress_label:
                        sys.stderr.write(
                            f"\r{progress_label}: {draws_used}/{num_draws}"
                        )
                        sys.stderr.flush()
                if _budget_null_stable_enough(
                    draws_used,
                    mean_units,
                    m2_units,
                    min_draws=min_draws,
                    abs_tol=stability_abs_tol,
                    rel_tol=stability_rel_tol,
                ):
                    break

    if progress_label:
        sys.stderr.write("\n")
        sys.stderr.flush()

    null_units_sd = float(np.sqrt(max(m2_units / float(max(draws_used - 1, 1)), 0.0)))
    null_units_stderr = float(
        np.sqrt(
            max(m2_units / float(max(draws_used - 1, 1)), 0.0)
            / float(max(draws_used, 1))
        )
    )
    null_tail_occupancy_sd = float(
        np.sqrt(max(m2_tail_occupancy / float(max(draws_used - 1, 1)), 0.0))
    )
    null_tail_occupancy_stderr = float(
        np.sqrt(
            max(m2_tail_occupancy / float(max(draws_used - 1, 1)), 0.0)
            / float(max(draws_used, 1))
        )
    )

    return {
        "observed_scores": observed_scores_.astype(np.float64),
        "null_center": float(null_center),
        "null_scale": float(null_scale),
        "null_positive_mass": float(mean_mass),
        "null_positive_units": float(mean_units),
        "null_positive_fraction": float(mean_fraction),
        "null_positive_units_sd": float(null_units_sd),
        "null_positive_units_stderr": float(null_units_stderr),
        "null_threshold": float(null_threshold),
        "null_tail_occupancy": float(mean_tail_occupancy),
        "null_tail_occupancy_sd": float(null_tail_occupancy_sd),
        "null_tail_occupancy_stderr": float(null_tail_occupancy_stderr),
        "negative_support_size": int(negative_magnitudes.size),
        "negative_fraction": float(
            negative_magnitudes.size / max(int(finite_null_scores.size), 1)
        ),
        "num_null_draws": int(draws_used),
        "max_null_draws": int(num_draws),
        "adaptive_stop": bool(draws_used < num_draws),
        "wild_bandwidth": int(bandwidth),
        "wild_process": "bartlett_multiplier",
        "null_method": "dependent_wild_residual_bootstrap",
        "null_reference_mean_positive_consensus": float(np.mean(positive_consensus)),
        "null_reference_max_positive_consensus": float(np.max(positive_consensus)),
    }


def estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
    centered_matrix: np.ndarray,
    observed_scores: np.ndarray | None = None,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    dependence_lag_hint: int | None = None,
    num_null_draws: int = 25,
    random_seed: int = 0,
    progress_label: str | None = None,
    num_processes: int = 1,
    return_details: bool = False,
) -> float | Tuple[float, Dict[str, Any]]:
    r"""Estimate a conservative enriched fraction from a wild-bootstrap null.

    We fit a null residual template by subtracting the positive
    part of the fitted locus effect, then bootstrap that template with a
    short-range dependent multiplier process. With :math:`\mu_0` and
    :math:`\sigma_0` estimated from the fitted null score field, the observed
    soft count is

    .. math::

       m_{\mathrm{obs}} =
       \frac{1}{n}\sum_j \frac{(s_j - \mu_0)_+}{\max(\sigma_0, 10^{-6})},

    and the null contribution is the bootstrap average

    .. math::

       m_0 =
       \mathbb{E}_{*}\left[
       \frac{(S_j^{*} - \mu_0)_+}{\max(\sigma_0, 10^{-6})}
       \right].

    Let :math:`t_0 = \mu_0 + 2\sigma_0`. we estimate the raw chromosome
    budget from above that fitted null threshold. If
    :math:`\hat p_{\mathrm{obs}} = n^{-1}\sum_j I(s_j > t_0)` and
    :math:`\hat p_0` is the bootstrap average of the same statistic under the
    fitted null, with standard deviation :math:`\hat \sigma_{p_0}`,
    then ROCCO uses

    .. math::

       \hat \pi_1 =
       \operatorname{clip}(
       \hat p_{\mathrm{obs}} - \hat p_0 - \hat \sigma_{p_0},
       0,
       1)

    as the raw enriched fraction for the budget EB step.
    """
    centered = np.asarray(centered_matrix, dtype=np.float64)
    if centered.ndim == 1:
        centered = centered[np.newaxis, :]
    if centered.ndim != 2:
        raise ValueError("`centered_matrix` must be one- or two-dimensional")

    _, n_loci = centered.shape
    if n_loci <= 0:
        raise ValueError("`centered_matrix` must contain at least one locus")

    null_meta = _estimate_wild_bootstrap_score_null(
        centered,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
        observed_scores=observed_scores,
        dependence_lag_hint=dependence_lag_hint,
        num_null_draws=num_null_draws,
        random_seed=random_seed,
        progress_label=progress_label,
        num_processes=num_processes,
    )
    observed_scores_ = np.asarray(null_meta["observed_scores"], dtype=np.float64)
    null_center = float(null_meta["null_center"])
    null_scale = float(null_meta["null_scale"])
    null_soft_scale = float(max(null_scale, 1.0e-6))
    residual_scores = observed_scores_ - null_center
    # Positive-part score mass above the fitted bootstrap null:
    # :math:`(s_j - \mu_0)_+ / \sigma_0`.
    observed_excess = np.clip(residual_scores, 0.0, None)
    observed_negative = np.clip(-residual_scores, 0.0, None)
    observed_soft_counts = observed_excess / null_soft_scale
    observed_positive_fraction = float(np.mean(observed_excess > 0.0))
    observed_negative_fraction = float(np.mean(observed_negative > 0.0))
    observed_excess_mass = float(np.mean(observed_excess))
    observed_excess_units = float(np.mean(observed_soft_counts))
    null_excess_mass = float(null_meta["null_positive_mass"])
    null_excess_units = float(null_meta["null_positive_units"])
    null_excess_units_sd = float(null_meta["null_positive_units_sd"])
    null_threshold = float(null_meta["null_threshold"])
    observed_tail_occupancy = float(np.mean(observed_scores_ > null_threshold))
    null_tail_occupancy = float(null_meta["null_tail_occupancy"])
    null_tail_occupancy_sd = float(null_meta["null_tail_occupancy_sd"])
    ess_max_lag = _resolve_budget_ess_max_lag(
        n_loci,
        dependence_lag_hint=dependence_lag_hint,
    )
    effective_total_count, tau_int, ess_lags_used = _estimate_effective_sample_size(
        observed_soft_counts,
        max_lag=ess_max_lag,
    )
    nonnull_fraction = float(
        np.clip(
            observed_tail_occupancy - null_tail_occupancy - null_tail_occupancy_sd,
            0.0,
            1.0,
        )
    )
    if (
        not np.isfinite(nonnull_fraction)
        or not np.isfinite(effective_total_count)
        or not np.isfinite(tau_int)
    ):
        raise ValueError("Budget initialization produced non-finite values")
    details = {
        "observed_positive_fraction": float(observed_positive_fraction),
        "observed_negative_fraction": float(observed_negative_fraction),
        "null_positive_fraction": float(null_meta["null_positive_fraction"]),
        "observed_excess_mass": float(observed_excess_mass),
        "null_excess_mass": float(null_excess_mass),
        "observed_excess_units": float(observed_excess_units),
        "null_excess_units": float(null_excess_units),
        "null_excess_units_sd": float(null_excess_units_sd),
        "null_excess_units_stderr": float(null_meta["null_positive_units_stderr"]),
        "null_threshold": float(null_threshold),
        "observed_tail_occupancy": float(observed_tail_occupancy),
        "null_tail_occupancy": float(null_tail_occupancy),
        "null_tail_occupancy_sd": float(null_tail_occupancy_sd),
        "null_tail_occupancy_stderr": float(null_meta["null_tail_occupancy_stderr"]),
        "null_center": float(null_center),
        "null_scale": float(null_scale),
        "nonnull_fraction": float(nonnull_fraction),
        "effective_count": float(nonnull_fraction * effective_total_count),
        "effective_total_count": float(effective_total_count),
        "autocorrelation_time": float(tau_int),
        "ess_max_lag": float(ess_max_lag),
        "ess_lags_used": float(ess_lags_used),
        "num_loci": float(n_loci),
        "negative_support_size": float(null_meta["negative_support_size"]),
        "negative_fraction": float(null_meta["negative_fraction"]),
        "num_null_draws": float(null_meta["num_null_draws"]),
        "max_null_draws": float(null_meta["max_null_draws"]),
        "adaptive_stop": bool(null_meta["adaptive_stop"]),
        "wild_bandwidth": float(null_meta["wild_bandwidth"]),
        "wild_process": str(null_meta["wild_process"]),
        "null_method": str(null_meta["null_method"]),
        "null_reference_mean_positive_consensus": float(
            null_meta["null_reference_mean_positive_consensus"]
        ),
        "null_reference_max_positive_consensus": float(
            null_meta["null_reference_max_positive_consensus"]
        ),
    }
    if return_details:
        return nonnull_fraction, details
    return nonnull_fraction


def _estimate_wild_bootstrap_direct_score_null(
    score_track: np.ndarray,
    dependence_lag_hint: int | None = None,
    num_null_draws: int = 25,
    random_seed: int = 0,
    progress_label: str | None = None,
    min_null_draws: int | None = None,
    stability_abs_tol: float = 5.0e-3,
    stability_rel_tol: float = 5.0e-2,
) -> dict[str, float | int | str | np.ndarray]:
    r"""Draw nulls from generic score tracks by a dependent wild bootstrap"""
    scores = np.asarray(score_track, dtype=np.float64)
    if scores.ndim != 1:
        raise ValueError("`score_track` must be one-dimensional")
    if scores.size == 0:
        raise ValueError("`score_track` must contain at least one locus")

    observed_scores = scores.astype(np.float64, copy=False)
    positive_consensus = np.clip(observed_scores, 0.0, None)
    residual_template = observed_scores - positive_consensus
    null_reference_scores = residual_template
    null_center = float(np.median(null_reference_scores))
    null_reference_residuals = null_reference_scores - null_center
    negative_reference = null_reference_residuals[null_reference_residuals <= 0.0]
    if negative_reference.size == 0:
        negative_magnitudes = np.abs(null_reference_residuals)
    else:
        negative_magnitudes = -negative_reference
    if negative_magnitudes.size == 0:
        negative_magnitudes = np.array([0.0], dtype=np.float64)
    mirrored_reference = np.concatenate((-negative_magnitudes, negative_magnitudes))
    null_scale = float(_robust_scale(mirrored_reference))
    if not np.isfinite(null_center) or not np.isfinite(null_scale):
        raise ValueError("Direct-score budget null fit produced non-finite values")
    null_soft_scale = float(max(null_scale, 1.0e-6))
    null_threshold = float(null_center + (2.0 * null_scale))

    bandwidth = _resolve_budget_bootstrap_bandwidth(
        observed_scores.size,
        dependence_lag_hint=dependence_lag_hint,
    )
    kernel = _build_budget_bootstrap_kernel(bandwidth)
    num_draws = int(max(1, num_null_draws))
    min_draws = int(
        min(num_draws, max(4, 8 if min_null_draws is None else min_null_draws))
    )
    draws_used = 0
    mean_mass = 0.0
    m2_mass = 0.0
    mean_units = 0.0
    m2_units = 0.0
    mean_fraction = 0.0
    m2_fraction = 0.0
    mean_tail_occupancy = 0.0
    m2_tail_occupancy = 0.0
    rng = np.random.default_rng(int(random_seed))

    # For each draw, generate a short-range dependency-inducing conv. kernel for the wild bootstrap
    for draw_id in range(num_draws):
        wild_weights = _generate_dependent_wild_weights(
            observed_scores.size,
            kernel=kernel,
            rng=rng,
        )
        # convolve w/ template
        bootstrap_scores = residual_template * wild_weights
        # center null draws
        bootstrap_residual_scores = bootstrap_scores - null_center
        bootstrap_positive = np.clip(bootstrap_residual_scores, 0.0, None)
        null_mass = float(np.mean(bootstrap_positive))
        null_units = float(np.mean(bootstrap_positive / null_soft_scale))
        null_fraction = float(np.mean(bootstrap_positive > 0.0))
        null_tail_occupancy = float(np.mean(bootstrap_scores > null_threshold))

        draws_used, mean_mass, m2_mass = _update_running_moments(
            draws_used,
            mean_mass,
            m2_mass,
            null_mass,
        )
        _, mean_units, m2_units = _update_running_moments(
            draws_used - 1,
            mean_units,
            m2_units,
            null_units,
        )
        _, mean_fraction, m2_fraction = _update_running_moments(
            draws_used - 1,
            mean_fraction,
            m2_fraction,
            null_fraction,
        )
        _, mean_tail_occupancy, m2_tail_occupancy = _update_running_moments(
            draws_used - 1,
            mean_tail_occupancy,
            m2_tail_occupancy,
            null_tail_occupancy,
        )
        if progress_label:
            sys.stderr.write(f"\r{progress_label}: {draws_used}/{num_draws}")
            sys.stderr.flush()
        if _budget_null_stable_enough(
            draws_used,
            mean_units,
            m2_units,
            min_draws=min_draws,
            abs_tol=stability_abs_tol,
            rel_tol=stability_rel_tol,
        ):
            break

    if progress_label:
        sys.stderr.write("\n")
        sys.stderr.flush()

    null_units_sd = float(np.sqrt(max(m2_units / float(max(draws_used - 1, 1)), 0.0)))
    null_units_stderr = float(
        np.sqrt(
            max(m2_units / float(max(draws_used - 1, 1)), 0.0)
            / float(max(draws_used, 1))
        )
    )
    null_tail_occupancy_sd = float(
        np.sqrt(max(m2_tail_occupancy / float(max(draws_used - 1, 1)), 0.0))
    )
    null_tail_occupancy_stderr = float(
        np.sqrt(
            max(m2_tail_occupancy / float(max(draws_used - 1, 1)), 0.0)
            / float(max(draws_used, 1))
        )
    )

    # build metadata dict with null estimates and diagnostics
    return {
        "observed_scores": observed_scores.astype(np.float64),
        "null_center": float(null_center),
        "null_scale": float(null_scale),
        "null_positive_mass": float(mean_mass),
        "null_positive_units": float(mean_units),
        "null_positive_fraction": float(mean_fraction),
        "null_positive_units_sd": float(null_units_sd),
        "null_positive_units_stderr": float(null_units_stderr),
        "null_threshold": float(null_threshold),
        "null_tail_occupancy": float(mean_tail_occupancy),
        "null_tail_occupancy_sd": float(null_tail_occupancy_sd),
        "null_tail_occupancy_stderr": float(null_tail_occupancy_stderr),
        "negative_support_size": int(negative_magnitudes.size),
        "negative_fraction": float(
            negative_magnitudes.size / max(int(null_reference_scores.size), 1)
        ),
        "num_null_draws": int(draws_used),
        "max_null_draws": int(num_draws),
        "adaptive_stop": bool(draws_used < num_draws),
        "wild_bandwidth": int(bandwidth),
        "wild_process": "bartlett_multiplier",
        "null_method": "dependent_wild_score_bootstrap",
        "null_reference_mean_positive_consensus": float(np.mean(positive_consensus)),
        "null_reference_max_positive_consensus": float(np.max(positive_consensus)),
    }


def estimate_budget_nonnull_fraction_from_score_track(
    score_track: np.ndarray,
    dependence_lag_hint: int | None = None,
    num_null_draws: int = 25,
    random_seed: int = 0,
    progress_label: str | None = None,
    num_processes: int = 1,
    return_details: bool = False,
) -> float | Tuple[float, Dict[str, Any]]:
    r"""Estimate a conservative enriched fraction directly from a score track."""
    _ = int(max(1, num_processes))
    scores = np.asarray(score_track, dtype=np.float64)
    if scores.ndim != 1:
        raise ValueError("`score_track` must be one-dimensional")
    if scores.size == 0:
        raise ValueError("`score_track` must contain at least one locus")

    null_meta = _estimate_wild_bootstrap_direct_score_null(
        scores,
        dependence_lag_hint=dependence_lag_hint,
        num_null_draws=num_null_draws,
        random_seed=random_seed,
        progress_label=progress_label,
    )
    observed_scores = np.asarray(null_meta["observed_scores"], dtype=np.float64)
    null_center = float(null_meta["null_center"])
    null_scale = float(null_meta["null_scale"])
    null_soft_scale = float(max(null_scale, 1.0e-6))
    residual_scores = observed_scores - null_center
    observed_excess = np.clip(residual_scores, 0.0, None)
    observed_negative = np.clip(-residual_scores, 0.0, None)
    observed_soft_counts = observed_excess / null_soft_scale
    observed_positive_fraction = float(np.mean(observed_excess > 0.0))
    observed_negative_fraction = float(np.mean(observed_negative > 0.0))
    observed_excess_mass = float(np.mean(observed_excess))
    observed_excess_units = float(np.mean(observed_soft_counts))
    null_excess_mass = float(null_meta["null_positive_mass"])
    null_excess_units = float(null_meta["null_positive_units"])
    null_excess_units_sd = float(null_meta["null_positive_units_sd"])
    null_threshold = float(null_meta["null_threshold"])
    observed_tail_occupancy = float(np.mean(observed_scores > null_threshold))
    null_tail_occupancy = float(null_meta["null_tail_occupancy"])
    null_tail_occupancy_sd = float(null_meta["null_tail_occupancy_sd"])
    ess_max_lag = _resolve_budget_ess_max_lag(
        scores.size,
        dependence_lag_hint=dependence_lag_hint,
    )
    effective_total_count, tau_int, ess_lags_used = _estimate_effective_sample_size(
        observed_soft_counts,
        max_lag=ess_max_lag,
    )
    nonnull_fraction = float(
        np.clip(
            observed_tail_occupancy - null_tail_occupancy - null_tail_occupancy_sd,
            0.0,
            1.0,
        )
    )
    if (
        not np.isfinite(nonnull_fraction)
        or not np.isfinite(effective_total_count)
        or not np.isfinite(tau_int)
    ):
        raise ValueError(
            "Direct-score budget initialization produced non-finite values"
        )

    details = {
        "observed_positive_fraction": float(observed_positive_fraction),
        "observed_negative_fraction": float(observed_negative_fraction),
        "null_positive_fraction": float(null_meta["null_positive_fraction"]),
        "observed_excess_mass": float(observed_excess_mass),
        "null_excess_mass": float(null_excess_mass),
        "observed_excess_units": float(observed_excess_units),
        "null_excess_units": float(null_excess_units),
        "null_excess_units_sd": float(null_excess_units_sd),
        "null_excess_units_stderr": float(null_meta["null_positive_units_stderr"]),
        "null_threshold": float(null_threshold),
        "observed_tail_occupancy": float(observed_tail_occupancy),
        "null_tail_occupancy": float(null_tail_occupancy),
        "null_tail_occupancy_sd": float(null_tail_occupancy_sd),
        "null_tail_occupancy_stderr": float(null_meta["null_tail_occupancy_stderr"]),
        "null_center": float(null_center),
        "null_scale": float(null_scale),
        "nonnull_fraction": float(nonnull_fraction),
        "effective_count": float(nonnull_fraction * effective_total_count),
        "effective_total_count": float(effective_total_count),
        "autocorrelation_time": float(tau_int),
        "ess_max_lag": float(ess_max_lag),
        "ess_lags_used": float(ess_lags_used),
        "num_loci": float(scores.size),
        "negative_support_size": float(null_meta["negative_support_size"]),
        "negative_fraction": float(null_meta["negative_fraction"]),
        "num_null_draws": float(null_meta["num_null_draws"]),
        "max_null_draws": float(null_meta["max_null_draws"]),
        "adaptive_stop": bool(null_meta["adaptive_stop"]),
        "wild_bandwidth": float(null_meta["wild_bandwidth"]),
        "wild_process": str(null_meta["wild_process"]),
        "null_method": str(null_meta["null_method"]),
        "null_reference_mean_positive_consensus": float(
            null_meta["null_reference_mean_positive_consensus"]
        ),
        "null_reference_max_positive_consensus": float(
            null_meta["null_reference_max_positive_consensus"]
        ),
    }
    if return_details:
        return nonnull_fraction, details
    return nonnull_fraction


def estimate_budget_nonnull_fraction_from_empirical_null(
    centered_matrix: np.ndarray,
    observed_scores: np.ndarray | None = None,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    dependence_lag_hint: int | None = None,
    num_null_draws: int = 25,
    random_seed: int = 0,
    progress_label: str | None = None,
    num_processes: int = 1,
    return_details: bool = False,
) -> float | Tuple[float, Dict[str, Any]]:
    r"""Wrapper for the wild-bootstrap budget estimator."""
    return estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
        centered_matrix,
        observed_scores=observed_scores,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
        dependence_lag_hint=dependence_lag_hint,
        num_null_draws=num_null_draws,
        random_seed=random_seed,
        progress_label=progress_label,
        num_processes=num_processes,
        return_details=return_details,
    )


def estimate_budget_nonnull_fraction_from_resampled_null(
    centered_matrix: np.ndarray,
    observed_scores: np.ndarray | None = None,
    lower_bound_z: float = 1.0,
    prior_df: float = 5.0,
    min_effect: float | None = None,
    precision_floor_ratio: float = 0.01,
    num_null_draws: int = 25,
    mean_block_length: int | None = None,
    null_threshold_scale: float = 1.0,
    random_seed: int = 0,
    progress_label: str | None = None,
    num_processes: int = 1,
    return_details: bool = False,
) -> float | Tuple[float, Dict[str, Any]]:
    r"""Wrapper for the wild-bootstrap budget estimator."""
    _ = (null_threshold_scale,)
    return estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
        centered_matrix,
        observed_scores=observed_scores,
        lower_bound_z=lower_bound_z,
        prior_df=prior_df,
        min_effect=min_effect,
        precision_floor_ratio=precision_floor_ratio,
        dependence_lag_hint=mean_block_length,
        num_null_draws=num_null_draws,
        random_seed=random_seed,
        progress_label=progress_label,
        num_processes=num_processes,
        return_details=return_details,
    )


def _estimate_local_peak_width(
    y_vals: np.ndarray,
    peak_idx: int,
    peak_search_radius: int = 2,
) -> float | None:
    n = int(y_vals.size)
    if n < 3:
        return None

    left_idx = max(0, int(peak_idx) - int(max(0, peak_search_radius)))
    right_idx = min(
        n,
        int(peak_idx) + int(max(0, peak_search_radius)) + 1,
    )
    if right_idx - left_idx < 1:
        return None

    local_peak_idx = left_idx + int(np.argmax(y_vals[left_idx:right_idx]))
    if local_peak_idx <= 0 or local_peak_idx >= (n - 1):
        return None

    peak_height = float(y_vals[local_peak_idx])
    if not np.isfinite(peak_height):
        return None

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        widths, _, _, _ = signal.peak_widths(
            y_vals,
            np.array([local_peak_idx], dtype=np.int64),
            rel_height=0.5,
        )

    if widths.size == 0:
        return None

    width_hat = float(widths[0])
    if not np.isfinite(width_hat) or width_hat <= 0.0:
        return None
    return float(max(1.0, width_hat))


def _estimate_local_peak_log_width_variance(
    y_vals: np.ndarray,
    peak_idx: int,
    half_window: int,
) -> tuple[float | None, float | None]:
    r"""Estimate local log-width and its uncertainty from a small smoothing ensemble.

    The width itself is measured with the heuristic ``scipy.signal.peak_widths`` after a local
    peak re-centering step. We then evaluate the same local peak on a few nearby
    smoothings of the local signal and use the spread of the resulting log-widths
    as a heteroskedastic uncertainty estimate for the EB summary downstream.
    """
    n = int(y_vals.size)
    if n < 8:
        return None, None

    # local window around the peak for width estimation and bootstrapping, with a minimum size of 5
    local_left = max(0, int(peak_idx) - int(max(2, half_window)))
    local_right = min(n, int(peak_idx) + int(max(2, half_window)) + 1)
    local_y = np.asarray(y_vals[local_left:local_right], dtype=np.float64)
    local_peak_idx = int(peak_idx) - local_left
    width_hat = _estimate_local_peak_width(
        local_y,
        local_peak_idx,
        peak_search_radius=2,
    )
    if width_hat is None:
        return None, None

    smooth_size = int(max(3, min(local_y.size - (1 - (local_y.size % 2)), 9)))
    if smooth_size % 2 == 0:
        smooth_size -= 1
    smooth_size = max(3, smooth_size)
    if smooth_size >= local_y.size:
        smooth_size = local_y.size - 1 if local_y.size % 2 == 0 else local_y.size
    smooth_size = max(3, smooth_size)
    candidate_windows = sorted(
        {
            1,
            smooth_size,
            min(
                local_y.size if (local_y.size % 2) == 1 else local_y.size - 1,
                smooth_size + 2,
            ),
        }
    )
    log_widths: list[float] = []
    for window_length in candidate_windows:
        if window_length <= 1:
            local_variant = local_y
        else:
            local_variant = signal.savgol_filter(
                local_y,
                window_length=window_length,
                polyorder=min(2, window_length - 1),
                mode="interp",
            )
        width_variant = _estimate_local_peak_width(
            np.asarray(local_variant, dtype=np.float64),
            local_peak_idx,
            peak_search_radius=2,
        )
        if width_variant is None:
            continue
        # log-scale widths st rare/wide peaks do not unduly dominate the variance estimate for more typical peaks
        log_widths.append(float(np.log(width_variant)))

    if len(log_widths) == 0:
        return None, None
    log_width_arr = np.asarray(log_widths, dtype=np.float64)
    log_width = float(np.median(log_width_arr))
    if log_width_arr.size < 2:
        return log_width, 1.0e-4
    return log_width, float(max(1.0e-4, np.var(log_width_arr, ddof=1)))


def _select_context_features(
    vals: np.ndarray,
    min_span: int | None = 3,
    max_span: int | None = 64,
    max_order: int = 5,
) -> Dict[str, Any]:
    r"""Detect stable local features for width and gamma initialization."""
    y = np.asarray(vals, dtype=np.float64)
    n = int(y.size)
    if n < 100:
        raise ValueError(
            "Input signal is too small for width-based context estimation."
        )

    min_span_ = 3 if min_span is None else int(min_span)
    max_span_ = (
        int(max(10, min(50, np.floor(np.log2(n + 1) * 2))))
        if max_span is None
        else int(max_span)
    )
    if max_span_ <= 0:
        raise ValueError("`max_span` must be positive.")

    y_pos = np.clip(y, 0.0, None)
    y_log = np.log1p(y_pos)
    pos_log = y_log[y > 0]
    if pos_log.size <= max(1, int(max_order)):
        raise ValueError(
            "Insufficient positive signal for width-based context estimation."
        )

    smooth_size = min(int(max(1, min_span_)), int(max_span_ / 2))
    # smooth out log signal before prominence-based peaks in scipy.find_peaks
    y_log_smooth = ndimage.uniform_filter1d(
        y_log,
        size=smooth_size,
        mode="nearest",
    )
    min_feature_count = int(max(1, (2 * np.log2(n + 1))))
    thr_log = float(np.mean(pos_log))
    start_order = int(max(1, max_order))
    best_order = 1
    best_features = np.array([], dtype=np.int64)
    best_score = -1.0

    for order in range(start_order, 0, -1):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            warnings.filterwarnings("ignore", category=UserWarning)
            features, props = signal.find_peaks(
                y_log_smooth,
                distance=int(max(1, order * max(1, min_span_))),
                prominence=1.0e-4,
            )

        prominences = props.get("prominences")
        prom_masked: np.ndarray | None = None
        if features.size:
            feature_mask = y_log[features] > thr_log
            if prominences is not None:
                feature_mask &= prominences > 0.0
            widths_half, _, _, _ = signal.peak_widths(
                y_log_smooth,
                features.astype(np.int64, copy=False),
                rel_height=0.5,
            )
            feature_mask &= np.isfinite(widths_half)
            feature_mask &= widths_half >= float(max(1, min_span_))
            if max_span_ > 0:
                feature_mask &= widths_half <= float(max_span_ * 2)
            if prominences is not None:
                prom_masked = prominences[feature_mask]
            features = features[feature_mask]

        if features.size == 0:
            score = -np.inf
        elif prom_masked is None:
            logger.warning(
                "Prominences missing during context-size estimation; using feature count as fallback score."
            )
            score = float(features.size)
        else:
            top_k = int(min(1000, prom_masked.size))
            if top_k <= 0:
                score = float(features.size)
            else:
                top = np.partition(prom_masked, -top_k)[-top_k:]
                score = float(np.sum(np.log1p(top)) / np.sqrt(float(top_k)))

        if score > best_score:
            best_score = score
            best_order = order
            best_features = features.astype(np.int64, copy=False)

    chosen_features = np.unique(best_features.astype(np.int64))
    chosen_features.sort()
    logger.info(
        "Chose peak distance order=%s with num_features=%s (score=%.3f).",
        best_order,
        int(chosen_features.size),
        best_score,
    )
    if chosen_features.size == 0:
        raise ValueError(
            "Could not identify stable features for width-based context estimation."
        )

    base_q = 0.05
    feature_baselines = np.empty(chosen_features.size, dtype=np.float64)
    for i, idx in enumerate(chosen_features):
        left = max(0, idx - max_span_)
        right = min(n - 1, idx + max_span_)
        left_q = float(np.quantile(y_log[left : idx + 1], base_q))
        right_q = float(np.quantile(y_log[idx : right + 1], base_q))
        feature_baselines[i] = max(left_q, right_q)

    feature_scores = y_log[chosen_features] - feature_baselines
    keep_mask = feature_scores > 0.0
    if np.any(keep_mask):
        chosen_features = chosen_features[keep_mask]
        feature_scores = feature_scores[keep_mask]

    keep_count = int(
        min(
            1000,
            feature_scores.size,
            max(min_feature_count, n // max(8, max_span_)),
        )
    )
    if keep_count <= 0:
        raise ValueError(
            "No informative features remained for width-based context estimation."
        )

    keep = np.argpartition(-feature_scores, keep_count - 1)[:keep_count]
    feature_index_array = np.unique(chosen_features[keep].astype(np.int64))
    feature_index_array.sort()
    return {
        "y_log": y_log.astype(np.float64, copy=False),
        "max_span": int(max_span_),
        "feature_indices": feature_index_array.astype(np.int64, copy=False),
        "feature_scores_log": feature_scores.astype(np.float64, copy=False),
        "num_features": int(feature_index_array.size),
        "best_order": int(best_order),
    }


def estimate_context_size(
    vals: np.ndarray,
    min_span: int | None = 3,
    max_span: int | None = 64,
    band_z: float = 1.0,
    max_order: int = 5,
    return_details: bool = False,
) -> tuple[int, int, int] | tuple[int, int, int, Dict[str, Any]]:
    r"""Estimate a characteristic peak width from local half-height widths.

    We detect prominent local features on a smoothed log-scale track, estimate a
    half-height width for each feature with ``scipy.signal.peak_widths``, derive
    a per-feature log-width uncertainty from a small smoothing ensemble, and then
    shrink widths with a simple normal-normal EB approach on log scale.
    """
    feature_meta = _select_context_features(
        vals,
        min_span=min_span,
        max_span=max_span,
        max_order=max_order,
    )
    y_log = np.asarray(feature_meta["y_log"], dtype=np.float64)
    max_span_ = int(feature_meta["max_span"])
    feature_index_array = np.asarray(
        feature_meta["feature_indices"],
        dtype=np.int64,
    )

    noise_window = int(min(max_span_, 32))
    s_hat_list: list[float] = []
    sigma2_list: list[float] = []
    for peak_idx in feature_index_array:
        log_width, sigma_s2 = _estimate_local_peak_log_width_variance(
            y_vals=y_log,
            peak_idx=int(peak_idx),
            half_window=noise_window,
        )
        if log_width is None or sigma_s2 is None:
            continue
        s_hat_list.append(log_width)
        sigma2_list.append(sigma_s2)

    if len(s_hat_list) == 0:
        raise ValueError("Failed to estimate widths from the detected features.")

    s_hat_arr = np.asarray(s_hat_list, dtype=np.float64)
    sigma2_arr = np.asarray(sigma2_list, dtype=np.float64)
    n_feat = int(s_hat_arr.size)
    var_s = float(np.var(s_hat_arr, ddof=1)) if n_feat > 1 else 0.0
    mean_sigma2 = float(np.mean(sigma2_arr))
    tau2_mom = float(max(0.0, var_s - mean_sigma2))
    tau2_max = float(max(1.0e-6, var_s + mean_sigma2) * 10.0)

    def _negative_log_likelihood(between_feature_var: float) -> float:
        tau_sq = float(max(0.0, between_feature_var))
        total_var = np.maximum(sigma2_arr + tau_sq, 1.0e-12)
        inv_total_var = 1.0 / total_var
        mu_hat = float(np.sum(inv_total_var * s_hat_arr) / np.sum(inv_total_var))
        residuals = s_hat_arr - mu_hat
        return 0.5 * float(
            np.sum(np.log(total_var)) + np.sum((residuals * residuals) * inv_total_var)
        )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        res = optimize.minimize_scalar(
            _negative_log_likelihood,
            bounds=(0.0, tau2_max),
            method="bounded",
            options={"xatol": 1.0e-4, "maxiter": 50},
        )

    if getattr(res, "success", False):
        tau_sq_hat = float(res.x)
        logger.info("Context-size tau^2 MLE plugin = %.6f", tau_sq_hat)
    else:
        tau_sq_hat = tau2_mom
        logger.warning(
            "Failed to optimize context-size tau^2; using MoM estimate tau^2=%.6f.",
            tau2_mom,
        )

    v_hat = np.maximum(sigma2_arr + tau_sq_hat, 1.0e-12)
    w_hat = 1.0 / v_hat
    mu_hat = float(np.sum(w_hat * s_hat_arr) / np.sum(w_hat))
    mu_var = float(1.0 / np.sum(w_hat))
    pred_std = float(np.sqrt(max(0.0, tau_sq_hat + mu_var)))
    log_lower = float(mu_hat - band_z * pred_std)
    log_upper = float(mu_hat + band_z * pred_std)
    point_estimate = float(np.exp(mu_hat))
    width_lower = max(1.0, float(np.exp(log_lower)))
    width_upper = max(1.0, float(np.exp(log_upper)), width_lower)

    if max_span_ > 0:
        point_estimate = min(point_estimate, float(max_span_))
        width_lower = min(width_lower, float(max_span_))
        width_upper = min(width_upper, float(max_span_))

    details = {
        "feature_indices": feature_index_array.astype(np.int64, copy=False),
        "feature_scores_log": np.asarray(
            feature_meta["feature_scores_log"],
            dtype=np.float64,
        ),
        "num_features": int(feature_index_array.size),
        "tau_sq_hat": float(tau_sq_hat),
        "context_log_mean": float(mu_hat),
        "context_log_mean_var": float(mu_var),
        "context_max_span": int(max_span_),
        "best_order": int(feature_meta["best_order"]),
    }
    logger.info(
        "Estimated context size on natural scale: point=%.4f lower=%.4f upper=%.4f",
        point_estimate,
        width_lower,
        width_upper,
    )
    if return_details:
        return (
            int(point_estimate),
            int(width_lower),
            int(width_upper),
            details,
        )
    return int(point_estimate), int(width_lower), int(width_upper)


def _estimate_gamma_signal_scale(
    scores: np.ndarray,
    feature_indices: np.ndarray,
    context_width: int,
) -> float:
    r"""Estimate the local positive score scale used for ``gamma_hat`` initialization."""
    scores_ = np.asarray(scores, dtype=np.float64)
    positive_scores = np.clip(scores_, 0.0, None)
    feature_indices_ = np.asarray(feature_indices, dtype=np.int64)
    local_mean_excess = []
    half_width = max(1, int(context_width))
    for idx in feature_indices_:
        left = max(0, int(idx) - half_width)
        right = min(scores_.size, int(idx) + half_width + 1)
        local_scores = positive_scores[left:right]
        if local_scores.size == 0:
            continue
        # Use :math:`\mathbb{E}[(s - q_{0.1})_+]` to measure local positive mass.
        baseline = float(np.quantile(local_scores, 0.1))
        mean_excess = float(np.mean(np.maximum(local_scores - baseline, 0.0)))
        if np.isfinite(mean_excess) and mean_excess > 0.0:
            local_mean_excess.append(mean_excess)

    if local_mean_excess:
        local_mean_excess_arr = np.asarray(
            local_mean_excess,
            dtype=np.float64,
        )
        signal_scale = float(np.median(local_mean_excess_arr))
        return signal_scale

    positive = positive_scores[positive_scores > 0]
    if positive.size == 0:
        raise ValueError("Cannot estimate gamma from scores without positive signal.")
    signal_scale = float(np.quantile(positive, 0.75))
    return signal_scale


def estimate_gamma_from_scores(
    scores: np.ndarray,
    gamma_scale: float = 0.5,
    min_span: int | None = 3,
    max_span: int | None = 64,
    band_z: float = 1.0,
    max_order: int = 5,
    context_size_summary: tuple[int, int, int] | None = None,
) -> tuple[float, Dict[str, float | int]]:
    r"""Initialize ``gamma`` from peak width and local positive score scale.

    The raw estimate is
    ``gamma_hat = gamma_scale * context_size_point * signal_scale`` where
    ``context_size_point`` is the point width estimate and
    ``signal_scale`` is the median local positive score excess around detected
    features.
    """
    scores_ = np.asarray(scores, dtype=np.float64)
    positive_scores = np.clip(scores_, 0.0, None)
    if context_size_summary is None:
        context_point, width_lower, width_upper, context_meta = estimate_context_size(
            positive_scores,
            min_span=min_span,
            max_span=max_span,
            band_z=band_z,
            max_order=max_order,
            return_details=True,
        )
        feature_indices = np.asarray(
            context_meta["feature_indices"],
            dtype=np.int64,
        )
    else:
        context_point = int(context_size_summary[0])
        width_lower = int(context_size_summary[1])
        width_upper = int(context_size_summary[2])
        feature_meta = _select_context_features(
            positive_scores,
            min_span=min_span,
            max_span=max_span,
            max_order=max_order,
        )
        feature_indices = np.asarray(
            feature_meta["feature_indices"],
            dtype=np.int64,
        )
        context_meta = {
            "best_order": int(feature_meta["best_order"]),
            "num_features": int(feature_indices.size),
        }

    signal_scale = _estimate_gamma_signal_scale(
        scores_,
        feature_indices,
        context_point,
    )

    gamma_hat = float(
        max(1.0e-6, float(gamma_scale) * float(max(1, context_point)) * signal_scale)
    )
    return gamma_hat, {
        "gamma_scale": float(gamma_scale),
        "signal_scale": float(signal_scale),
        "context_size_point": int(context_point),
        "context_size_lower": int(width_lower),
        "context_size_upper": int(width_upper),
        "num_features": int(context_meta["num_features"]),
        "feature_detection_order": int(context_meta["best_order"]),
    }


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
            "genome_wide_budget": float(init_center),
            "prior_strength": float(init_strength),
            "prior_dispersion": float(1.0 / (1.0 + alpha0 + beta0)),
            "min_prior_dispersion": 0.0,
            "observed_raw_budget_var": float(observed_raw_budget_var),
            "theoretical_min_raw_budget_var": float(theoretical_min_raw_budget_var),
            "prior_dispersion_at_floor": bool(False),
            "posterior_summary": "beta_quantile",
            "posterior_quantile": float(posterior_quantile_),
            "prior_fit_method": "single_chrom_default",
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
            "genome_wide_budget": float(pooled_rate),
            "prior_strength": float(prior_strength),
            "prior_dispersion": float(max(0.0, 1.0 / (1.0 + prior_strength))),
            "min_prior_dispersion": 0.0,
            "observed_raw_budget_var": float(observed_raw_budget_var),
            "theoretical_min_raw_budget_var": float(theoretical_min_raw_budget_var),
            "prior_dispersion_at_floor": bool(
                observed_raw_budget_var <= theoretical_min_raw_budget_var + 1.0e-12
            ),
            "posterior_summary": "beta_quantile",
            "posterior_quantile": float(posterior_quantile_),
            "prior_fit_method": "weak_pooled_prior",
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
        "genome_wide_budget": float(alpha_hat / (alpha_hat + beta_hat)),
        "prior_strength": float(prior_strength),
        "prior_dispersion": float(prior_dispersion),
        "min_prior_dispersion": 0.0,
        "observed_raw_budget_var": float(observed_raw_budget_var),
        "theoretical_min_raw_budget_var": float(theoretical_min_raw_budget_var),
        "prior_dispersion_at_floor": bool(dispersion_at_floor),
        "posterior_summary": "beta_quantile",
        "posterior_quantile": float(posterior_quantile_),
        "prior_fit_method": "beta_binomial_mle",
    }
