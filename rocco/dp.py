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

    Note this replaces the lp-based approach for a slightly relaxed problem and achieves linear time complexity and integral solutions.
    We still have to bisect for a suitable selection penalty if we want to meet a budget constraint, but that is cheap.

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


def calibrate_selection_penalty(
    scores: np.ndarray,
    switch_costs: np.ndarray,
    target_count: int,
    max_iter: int = 60,
) -> Tuple[float, np.ndarray, float, int]:
    r"""Find a selection penalty that yields a solution with an acceptable number of selected loci by bisection on DP solutions."""
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs_ = np.ascontiguousarray(switch_costs, dtype=np.float64)
    n = scores_.shape[0]
    if n == 0:
        raise ValueError("`scores` cannot be empty")
    target_count_ = int(max(0, min(target_count, n)))
    if target_count_ == n:
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs_,
            0.0,
        )
        return 0.0, solution, penalized_objective, selected_count

    lower = float(np.min(scores_) - np.sum(switch_costs_) - 1.0)
    upper = float(np.max(scores_) + np.sum(switch_costs_) + 1.0)

    lower_solution, lower_value, lower_count = solve_penalized_chain(
        scores_,
        switch_costs_,
        lower,
    )
    # find bounds before bisecting
    while lower_count <= target_count_:
        lower -= max(1.0, abs(lower))
        lower_solution, lower_value, lower_count = solve_penalized_chain(
            scores_,
            switch_costs_,
            lower,
        )
    # skip
    best_solution, best_value, best_count = solve_penalized_chain(
        scores_,
        switch_costs_,
        upper,
    )
    while best_count > target_count_:
        upper += max(1.0, abs(upper))
        best_solution, best_value, best_count = solve_penalized_chain(
            scores_,
            switch_costs_,
            upper,
        )

    # bisection
    for _ in range(max_iter):
        # solve @ midpoint
        midpoint = (lower + upper) / 2.0
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs_,
            midpoint,
        )

        # ---pick side---

        if selected_count > target_count_:
            lower = midpoint
            lower_solution = solution
            lower_value = penalized_objective
            lower_count = selected_count

        else:
            upper = midpoint
            best_solution = solution
            best_value = penalized_objective
            best_count = selected_count

    return upper, best_solution, best_value, best_count


def solve_chrom_exact(
    scores: np.ndarray,
    budget: Optional[float] = None,
    gamma: float = 1.0,
    selection_penalty: Optional[float] = None,
    return_details: bool = False,
) -> Tuple[np.ndarray, float] | Tuple[np.ndarray, float, Dict[str, float]]:
    r"""Solve one chromosome with the exact penalized-chain dynamic program.

    If ``selection_penalty`` is not supplied and ``budget`` is supplied, we
    find a penalty :math:`\lambda` that forces feasibility. Note that we
    do not have to saturate the budget entirely, only stay below it.

    If ``selection_penalty`` is supplied, skip that calibration step and
    solve the penalized chain directly with the supplied value.
    """
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs = build_switch_costs(
        scores_,
        gamma=gamma,
    )
    if selection_penalty is None:
        if budget is None:
            selection_penalty_ = 0.0
            solution, penalized_objective, selected_count = solve_penalized_chain(
                scores_,
                switch_costs,
                selection_penalty_,
            )
        else:
            target_count = int(np.floor(len(scores_) * float(budget)))
            (
                selection_penalty_,
                solution,
                penalized_objective,
                selected_count,
            ) = calibrate_selection_penalty(
                scores_,
                switch_costs,
                target_count=target_count,
            )
    else:
        selection_penalty_ = float(selection_penalty)
        solution, penalized_objective, selected_count = solve_penalized_chain(
            scores_,
            switch_costs,
            selection_penalty_,
        )

    objective = objective_value(solution, scores_, switch_costs)
    if not return_details:
        return solution.astype(np.uint8, copy=False), objective
    return (
        solution.astype(np.uint8, copy=False),
        objective,
        {
            "penalized_objective": float(penalized_objective),
            "selected_count": int(selected_count),
            "selected_fraction": float(selected_count / len(scores_)),
            "selection_penalty": float(selection_penalty_),
        },
    )
