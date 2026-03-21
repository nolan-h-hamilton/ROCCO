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
    scale_gamma: bool = False,
    beta: Optional[float] = None,
    denom_const: Optional[float] = None,
) -> np.ndarray:
    scores_ = np.asarray(scores, dtype=np.float64)
    if scores_.ndim != 1:
        raise ValueError("`scores` must be a one-dimensional array")
    if scores_.shape[0] <= 1:
        return np.zeros(0, dtype=np.float64)
    if not scale_gamma:
        return np.full(scores_.shape[0] - 1, float(gamma), dtype=np.float64)
    beta_ = 0.5 if beta is None else float(beta)
    denom_const_ = 1.0 if denom_const is None else float(denom_const)
    return float(gamma) / ((np.abs(np.diff(scores_, 1)) + denom_const_) ** beta_)


def _solve_penalized_chain_python(
    scores: np.ndarray,
    switch_costs: np.ndarray,
    selection_penalty: float,
) -> Tuple[np.ndarray, float, int]:
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs_ = np.ascontiguousarray(switch_costs, dtype=np.float64)
    if scores_.ndim != 1:
        raise ValueError("`scores` must be one-dimensional")
    n = scores_.shape[0]
    if n == 0:
        raise ValueError("`scores` cannot be empty")
    if n > 1 and switch_costs_.shape[0] != n - 1:
        raise ValueError("`switch_costs` must have length len(scores) - 1")

    # at locus i, we keep the best penalized objective ending in state 0 or 1.
    # The value recursion only needs the previous locus, so the objective and
    # selected-count updates stay O(1) in memory. We still keep one backpointer
    # array per state so we can reconstruct the final binary chain.
    bt0 = np.zeros(n, dtype=np.uint8)
    bt1 = np.zeros(n, dtype=np.uint8)
    prev0_val = 0.0
    prev0_count = 0
    prev1_val = float(scores_[0] - selection_penalty)
    prev1_count = 1

    # Move left to right and update the best penalized objective for each end state.
    # Note, tie-breaks favor solutions with _fewer_ selected loci.
    for i in range(1, n):

        switch_cost = float(switch_costs_[i - 1])
        stay0_val = prev0_val
        stay0_count = prev0_count
        # implicit: if previous state was 1, we pay to go to 0
        switch0_val = prev1_val - switch_cost
        switch0_count = prev1_count
        if switch0_val > stay0_val or (
            switch0_val == stay0_val and switch0_count < stay0_count
        ):
            new0_val = switch0_val
            new0_count = switch0_count
            bt0[i] = 1
        else:
            new0_val = stay0_val
            new0_count = stay0_count

        # add score_i - selection_penalty, then decide whether
        # that 1-segment is continuing or starts exactly at this locus.
        # If it starts here,
        # we have to pay the fragmentation/TV cost. Otherwise, we just
        # add the score and selection penalty to the previous state-1 value.
        stay1_val = prev1_val + float(scores_[i] - selection_penalty)
        stay1_count = prev1_count + 1
        switch1_val = prev0_val - switch_cost + float(scores_[i] - selection_penalty)
        switch1_count = prev0_count + 1

        # store which previous state gave the best state-1 value, breaking ties toward fewer selected loci.
        # backpointer array reconstructs the solution after the forward pass.
        if switch1_val > stay1_val or (
            switch1_val == stay1_val and switch1_count < stay1_count
        ):
            new1_val = switch1_val
            new1_count = switch1_count
            bt1[i] = 0
        else:
            new1_val = stay1_val
            new1_count = stay1_count
            bt1[i] = 1

        prev0_val = new0_val
        prev0_count = new0_count
        prev1_val = new1_val
        prev1_count = new1_count

    # final state
    if prev1_val > prev0_val or (prev1_val == prev0_val and prev1_count < prev0_count):
        best_val = prev1_val
        best_count = prev1_count
        state = 1
    else:
        best_val = prev0_val
        best_count = prev0_count
        state = 0

    # viterbi backtracking yields the solution
    solution = np.zeros(n, dtype=np.uint8)
    solution[-1] = state
    for i in range(n - 1, 0, -1):
        state = int(bt0[i]) if state == 0 else int(bt1[i])
        solution[i - 1] = state
    return solution, float(best_val), int(best_count)


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
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs_ = np.ascontiguousarray(switch_costs, dtype=np.float64)
    # use pure python if extension didn't build
    if _chain_dp is None:
        return _solve_penalized_chain_python(
            scores_,
            switch_costs_,
            float(selection_penalty),
        )
    # fast backend, no python
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
    beta: Optional[float] = None,
    denom_const: Optional[float] = None,
    scale_gamma: bool = False,
    selection_penalty: Optional[float] = None,
    return_details: bool = False,
) -> Tuple[np.ndarray, float] | Tuple[np.ndarray, float, Dict[str, float]]:
    r"""Solve one chromosome with the exact penalized-chain DP.

    If ``selection_penalty`` is not supplied and ``budget`` is supplied, we
    find a penalty :math:`\lambda` that forces feasibility. Note that we
    do not have to saturate the budget entirely, only stay below it.

    If ``selection_penalty`` is supplied, ROCCO skips that calibration step and
    solves the penalized chain directly with the supplied value.
    """
    scores_ = np.ascontiguousarray(scores, dtype=np.float64)
    switch_costs = build_switch_costs(
        scores_,
        gamma=gamma,
        scale_gamma=scale_gamma,
        beta=beta,
        denom_const=denom_const,
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
