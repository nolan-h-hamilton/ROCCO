from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

try:
    from . import _chain_dp
except ImportError:
    _chain_dp = None


def objective_value(
    solution: np.ndarray,
    scores: np.ndarray,
    switch_costs: np.ndarray | float,
) -> float:
    solution_ = np.asarray(solution, dtype=np.float64)
    scores_ = np.asarray(scores, dtype=np.float64)
    if np.isscalar(switch_costs):
        switch_costs_ = np.full(
            max(solution_.shape[0] - 1, 0),
            float(switch_costs),
            dtype=np.float64,
        )
    else:
        switch_costs_ = np.asarray(switch_costs, dtype=np.float64)
    penalty = 0.0
    if solution_.shape[0] > 1:
        penalty = float(switch_costs_ @ np.abs(np.diff(solution_, 1)))
    return float(-(scores_ @ solution_) + penalty)


def build_switch_costs(
    scores: np.ndarray,
    gamma: float = 1.0,
) -> np.ndarray:
    scores_ = np.asarray(scores, dtype=np.float64)
    if scores_.ndim != 1:
        raise ValueError("`scores` must be a one-dimensional array")
    if scores_.shape[0] <= 1:
        return np.zeros(0, dtype=np.float64)
    return np.full(scores_.shape[0] - 1, float(gamma), dtype=np.float64)


def solve_penalized_chain(
    scores: np.ndarray,
    switch_costs: np.ndarray,
    selection_penalty: float,
) -> Tuple[np.ndarray, float, int]:
    r"""Solve the penalized binary chain problem for one chromosome.

    This dynamic program solves the penalized binary chain problem in linear time
    and returns integral solutions.

    .. math::

       \max_{z \in \{0,1\}^n}
       \sum_{j=1}^n (s_j - \lambda) z_j
       \;-\;
       \sum_{j=1}^{n-1} c_j |z_{j+1} - z_j|,

    where :math:`s_j` are locus scores, :math:`c_j` are boundary penalties,
    and :math:`\lambda` is the direct penalty on selecting a locus.

    Similar problems are treated in  Johnson (2013, JCGS, doi:10.1080/10618600.2012.681238) and
    Madrid Padilla et al. (2017, JMLR 18). For the broader chain-DP context,
    see Forney (1973, Proc. IEEE, doi:10.1109/PROC.1973.9030).
    """
    if _chain_dp is None:
        raise RuntimeError("Make sure native C extensions are built and available")
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs_ = np.ascontiguousarray(switch_costs, dtype=np.float64)
    solution, penalized_objective, selected_count = _chain_dp.solve_penalized_chain(
        scores_,
        switch_costs_,
        float(selection_penalty),
    )
    return (
        np.asarray(solution, dtype=np.uint8),
        float(penalized_objective),
        int(selected_count),
    )


def solve_chrom_exact(
    scores: np.ndarray,
    budget: Optional[float] = None,
    gamma: float = 0.25,
    selection_penalty: Optional[float] = None,
    return_details: bool = False,
) -> Tuple[np.ndarray, float] | Tuple[np.ndarray, float, Dict[str, float]]:
    r"""Solve one chromosome with the exact penalized-chain dynamic program.

    If ``selection_penalty`` is not supplied and ``budget`` is supplied, use a
    quantile-derived soft penalty :math:`\lambda`.

    If ``selection_penalty`` is supplied, skip that calibration step and
    solve the penalized chain directly with the supplied value.
    """
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs = build_switch_costs(
        scores_,
        gamma=gamma,
    )
    budget_mode = "unpenalized"
    budget_ = None
    target_count = None
    if selection_penalty is None and budget is not None:
        budget_ = float(budget)
        if not np.isfinite(budget_) or budget_ < 0.0 or budget_ > 1.0:
            raise ValueError("`budget` must be finite and lie in [0, 1]")
        target_count = int(np.floor(len(scores_) * budget_))
        penalty_quantile = float(np.clip(1.0 - budget_, 0.0, 1.0))
        selection_penalty_ = float(
            max(0.0, np.quantile(scores_, penalty_quantile))
        )
        budget_mode = "soft_selection_penalty"
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs,
            selection_penalty_,
        )
    elif selection_penalty is None:
        selection_penalty_ = 0.0
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs,
            selection_penalty_,
        )
    else:
        selection_penalty_ = float(selection_penalty)
        if not np.isfinite(selection_penalty_):
            raise ValueError("`selection_penalty` must be finite")
        budget_mode = "manual_selection_penalty"
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs,
            selection_penalty_,
        )

    objective = objective_value(solution, scores_, switch_costs)
    if not return_details:
        return solution.astype(np.uint8, copy=False), objective
    details = {
        "penalized_objective": float(penalized_objective),
        "selected_count": int(selected_count),
        "selected_fraction": float(selected_count / len(scores_)),
        "selection_penalty": float(selection_penalty_),
        "budget_mode": budget_mode,
        "soft_budget_penalty": float(selection_penalty_),
    }
    if budget_ is not None:
        details["budget"] = float(budget_)
        details["budget_target_count"] = int(target_count)
    return (
        solution.astype(np.uint8, copy=False),
        objective,
        details,
    )
