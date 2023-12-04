r"""
# loci.py
ROCCO creates a Locus object (see `locus.py`) for each
$\ell_i, \forall i = 1 \ldots n$. These Locus objects are then linked
together as a Loci object, defined in this code.

Scoring: `Loci.score_loci()`

Optimization: `Loci.rocco_lp()`

Integral solutions: `Loci.run_rr()`

#### [Project Homepage](https://github.com/nolan-h-hamilton/ROCCO/)
"""

import copy
import math
import random
import numpy as np
import cvxpy as cp
from scipy import stats
from .locus import Locus

class LociIterator:
    """
    Iterator for traversing the Loci object

    Args:
        head (Locus): first Locus object in the Loci object

    Returns:
        LociIterator: an iterator object

    Raises:
        StopIteration: raised if end of genomic region is reached during iteration
    """
    def __init__(self, head):
        self.current = head

    def __iter__(self):
        return self

    def __next__(self):
        if not self.current:
            raise StopIteration

        item = copy.copy(self.current)
        self.current = self.current.right
        return item

class Loci:

    def __init__(self, head=None, tail=None):
        self.head = head
        self.tail = tail

    def __iter__(self):
        return LociIterator(self.head)

    def is_empty(self) -> bool:
        """
        Check if the Loci object is empty

        Returns:
            bool: True if Loci object is empty--False otherwise.
        """
        if self.head is None and self.tail is None:
            return True
        return False

    def append_locus(self,locus):
        """
        Append a new locus as the tail of the Loci object.

        Args:
            locus (Locus): The locus to be appended.

        Returns:
            None
        """
        if  self.is_empty():
            self.head = locus
            self.tail = locus
            locus.left=None
            locus.right=None
        else:
            tail_cpy = self.tail
            self.tail = locus
            locus.right = None
            locus.left = tail_cpy
            tail_cpy.right=locus

    def combine(self,head_,tail_):
        """Combine contiguous accessible regions

        This function takes the starting and ending node
        for a sequence of contiguous loci (i.e., a sequence
        such that all loci satisfy `Locus.accessible = 1`) and
        merges them into a single locus object. It defines
        the new locus's signal vector as the sum
        of signal vectors among loci in the contiguous
        sequence.

        Args:
            head_ (Locus): first locus object in sequence of contiguous
              selected loci
            tail_ (Locus): last locus object in sequence of contiguous
              selected loci
        """
        cpy = head_
        accessible_ = 0
        size_ = 0
        sig_matrix = []
        subloc = Loci()

        accessible_= 0
        if head_.accessible == 1:
            accessible_ = 1

        while head_ is not None and head_.position <= tail_.position:
            size_ += head_.size
            sig_matrix.append(head_.sig_data)
            subloc.append_locus(copy.copy(head_))
            head_ = head_.right

        sig_matrix = np.array([np.array(sig) for sig in sig_matrix]).transpose()
        sig_data_ = np.sum(sig_matrix,axis=1)

        # create new locus
        new_loc = Locus(position=cpy.position,left=cpy.left,right=head_,
                        size=size_, parent=None, subloci=subloc,
                        sig_data=sig_data_,accessible=accessible_)

        if head_ is None:
            self.tail = new_loc
        else:
            head_.left = new_loc
            cpy.right = new_loc

        if cpy.left is not None:
            cpy.left.right = new_loc
        else:
            self.head = new_loc
        s_head_ = subloc.head
        while s_head_ is not None and s_head_.position <= tail_.position:
            s_head_.parent=new_loc
            s_head_ = s_head_.right

    def as_list(self) -> list:
        """
        Get a list representation of the Loci object.

        Returns:
            list: list representation of the Loci object, where the Locus elements
                are ordered as in the Loci object.
        """
        head_ = self.head
        list_repr = []
        while head_ is not None:
            list_repr.append(copy.copy(head_))
            head_ = head_.right
        return list_repr

    def combine_selected(self):
        r"""
        Find adjacent accessible regions and combine them into a single locus

        This function identifies sequences of locus objects, $\ell_i$, satisfying
        $\ell_i$`.accessible > 0` and combines them into a single locus via `Loci.combine()`.


        Notes:
            - `Loci.combine()` defines signal data for the created locus as the vector sum of
                signal data for each Locus object comprising the contiguous accessible region.

        Returns:
            None
        """
        head_ = self.head
        while head_ is not None:
            to_combine = []
            head_cpy = head_
            while head_cpy is not None and head_cpy.accessible > 0 :
                to_combine.append(head_cpy)
                head_cpy = head_cpy.right
            if len(to_combine) > 1:
                self.combine(to_combine[0], to_combine[-1])
            head_ = head_cpy
            if head_ is not None:
                head_ = head_.right

    def get_sig_mat(self) -> np.ndarray:
        r"""
        Builds the signal matrix $\mathbf{S}_{chr} \in \mathbb{R}^{K \times n}$ as defined in the paper.

        Loop through each of the $n$ Locus object in `self` (type: Loci),
        where each `Locus.sig_data` is a $K$-length list containing
        the signal values for each of the $K$ samples at the respective Locus.
        Append to `sig_mat` and then take the transpose to return a $K \times n$
        np.ndarray.

        Returns:
            S_chr (np.ndarray): $K \times n$ sample-by-locus track matrix.
        """
        sig_mat = []
        head_ = self.head
        while head_ is not None:
            sig_mat.append(head_.sig_data)
            head_ = head_.right
        S_chr = np.array(sig_mat).T
        return S_chr

    def score_loci(self, tau: float = 0.0, c1: float = 1.0, c2: float = 1.0, c3: float = 1.0, eps_l=1e-4) -> np.ndarray:
        r"""
        compute locus score vector, i.e., $\mathcal{S}(i), \forall i = 1 \ldots n$
        given the $K \times n$ signal matrix $\mathbf{S}_{chr}$.

        Args:
            tau (float): median signal threshold (see $\tau$ in paper)
            c1 (float): weight for the first term in  $\mathcal{S}(i)$,
                (median signal) $g_1$
            c2 (float): weight for second term in $\mathcal{S}(i)$,
                (MAD dispersion) $g_2$
            c3 (float): weight for third term in  $\mathcal{S}(i)$, $g_3$

        Returns:
            np.ndarray: the $n$-dimensional score vector $\mathcal{S}(i), \forall i = 1 \ldots n$
        """
        def g3_vec(vvec) -> np.ndarray:
            g3_vals = np.zeros(len(vvec))
            for i,Loc in enumerate(vvec):
                if i == 0:
                    vel = (abs(vvec[i] - vvec[i+1]))
                if i == len(vvec)-1:
                    vel = abs(vvec[i] - vvec[i-1])
                else:
                    vel = max(abs(vvec[i] - vvec[i+1]),
                        abs(vvec[i] - vvec[i-1]))

                vel /= vvec[i]+1
                g3_vals[i] = vel
            return g3_vals

        sig_mat = self.get_sig_mat().T
        med_vec = np.median(sig_mat,axis=1) # g_1
        mad_vec = stats.median_abs_deviation(sig_mat,axis=1) # g_2
        g3_vec_ = g3_vec(med_vec) # g_3
        s_vec = c1*med_vec - c2*mad_vec + c3*g3_vec_

        for i, val in enumerate(s_vec):
            if med_vec[i] <= tau:
                s_vec[i] = 0
        return s_vec - eps_l


    def run_rr(self, lp_sol, N, loci_scores, budget, gam, eps = 1e-8, sumsq_penalty=None) -> np.ndarray:
        r"""
        Carry out the $\texttt{RR}$ rounding procedure to return a
        'good' integral solution with reference to the LP solution

        Args:
            lp_sol (numpy.ndarray): The LP solution
            N (int): Number of RR solutions to evaluate. If $N \leq 0$, the
                floor-epsilon protocol is applied.
            loci_scores (np.ndarray): $\mathcal{S}$
            budget (float): budget parameter
            gam (float): weight for $\sum_{i=1}^{i=n-1} |\ell_i - \ell_{i+1}|$ (discontig.) penalty
            eps (float): `init_sol = np.floor(lp_sol + eps)`. Decreased iteratively if initial `eps` does
                not yield a feasible solution.
            sumsq_penalty: Experimental. `None` by default. If not `None`, add $$\mathsf{sumsq\_penalty}\cdot \mathbf{\ell}^{T}\mathbf{\ell}$$ to
                the objective function. In this case, the $f(\mathbf{\ell})$ becomes strongly convex and the
                relaxed solution is guaranteed unique. Can be viewed as a penalty on the number of selections
                in supplement to the "hard" budget constraint.

        Returns:
        np.ndarray: the integral $\texttt{RR}$ solution, OR $\texttt{floor\_eps}$ solution if $N \leq 0$.
    """
        def obj(sol, loci_scores, gam, sumsq_penalty=sumsq_penalty):
            """
            Return numeric value of objective function given solution `sol`
            """
            if sumsq_penalty is None:
                return (-loci_scores@sol
                       + gam*np.sum(np.abs(np.diff(sol,1))))
            return (-loci_scores@sol
                       + gam*np.sum(np.abs(np.diff(sol,1))) + sumsq_penalty*(sol@sol.T))

        n = len(loci_scores)
        eps_cpy = eps
        # initialize as floor_eps solution
        init_sol = np.floor(lp_sol + eps)
        while np.sum(init_sol) > np.floor(n*budget):
            # loop guarantees the feasibility of solutions
            eps_cpy = eps_cpy/2
            init_sol = np.floor(lp_sol + eps_cpy)
        if N <= 0:
            return init_sol

        init_score = obj(sol=init_sol, loci_scores=loci_scores, gam=gam, sumsq_penalty=sumsq_penalty)

        nonint_loci = [i for i in range(len(lp_sol)) if lp_sol[i] > 0 and lp_sol[i] < 1]
        rr_sol = init_sol
        best_score = init_score
        for j in range(N):
            ell_rand_n = copy.copy(lp_sol)
            # for efficiency, only round `nonint_loci`
            for idx in nonint_loci:
                if random.random() <= lp_sol[idx]:
                    ell_rand_n[idx] = 1
                else:
                    ell_rand_n[idx] = 0

            score = obj(sol=ell_rand_n, loci_scores=loci_scores, gam=gam, sumsq_penalty=sumsq_penalty)

            is_feas = (np.sum(ell_rand_n) <= math.floor(n*budget))
            if is_feas and score < best_score:
                rr_sol = ell_rand_n
                best_score = score
        return rr_sol

    def rocco_lp(self, budget: float = .035, tau: float = 0, gam: float = 1, c1: float = 1,
             c2: float = 1, c3: float = 1, verbose_: bool = False, solver: str = "CLARABEL",
             N: int = 50, solver_reltol: float = 1.0e-8, solver_maxiter: float = 10000,
             solver_feastol: float = 1e-8, sumsq_penalty=None) -> cp.Problem:
        r"""
        Solve the LP-relaxation and apply the randomization procedure,
        $\texttt{RR}$ ($N > 0$), or $\texttt{floor\_eps}$ ($N \leq 0$) for an
        integral solution $\mathbf{\ell}^{*}$.

        Args:
            budget (float): determines budget constraint $\sum_i \ell_i \leq nb$
            tau (float): $\tau$ in $\mathcal{S}$
            gam (float): $\gamma$ in obj. function $f(\mathbf{\ell})$
            c1 (float): $c_1$ value in $\mathcal{S}$.
            c2 (float): $c_2$ value in $\mathcal{S}$
            c3 (float): $c_3$ value in $\mathcal{S}$
            N (int): $\texttt{RR}$ procedure iterations. See `Loci.run_rr()`
            verbose_ (bool): Verbosity flag for the solver
            solver (str): Defaults to `CLARABEL`. `ECOS, PDLP, MOSEK` are also supported.
            solver_reltol (float): acceptable relative optimality gap for solver termination.
            solver_feastol (float): the solver will consider solutions within `solver_feastol`
                of the constraint boundaries as feasible
            solver_maxiter (int): bounds the number of solving iterations before terminating
            sumsq_penalty: Experimental. If not `None`, add $$\mathsf{sumsq\_penalty}\cdot \mathbf{\ell}^{T}\mathbf{\ell}$$ to
                the objective function. In this case, the $f(\mathbf{\ell})$ becomes strongly convex and the
                relaxed *solution* is thus guaranteed unique. Can be viewed as a soft penalty on the number of selections
                in supplement to the "hard" budget constraint.
        Returns:
            cp.Problem: a CVXPY problem object
    """
        loci_scores = self.score_loci(tau, c1, c2, c3) # S_i, i=1...n
        n = len(loci_scores)
        # define problem in CVXPY
        ell = cp.Variable(n)
        z = cp.Variable(n-1) # aux. variables
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n),
                          z >= cp.diff(ell,1),
                          z >= -1*cp.diff(ell,1)]
        if sumsq_penalty is None:
            problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(z)),
                             constraints)
        else:
            problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(z) + sumsq_penalty*cp.sum_squares(ell)),
                             constraints)

        # Refer to `Solve method options` at
        # https://www.cvxpy.org/tutorial/advanced/index.html
        if solver.lower() not in [x.lower() for x in cp.installed_solvers()]:
            raise cp.error.SolverError(
                f'\nSolver {solver} is not available.\
                \nInstalled solvers: {cp.installed_solvers()}\
                \nhttps://www.cvxpy.org/tutorial/advanced/index.html\n')

        elif solver.lower() == 'clarabel':
            problem.solve(solver=cp.CLARABEL,
                          tol_gap_rel=solver_reltol,
                          tol_feas=solver_feastol,
                          max_iter=solver_maxiter,
                          verbose=verbose_)

        elif solver.lower() == "ecos":
            problem.solve(solver=cp.ECOS,
                          reltol=solver_reltol,
                          max_iters=solver_maxiter,
                          feastol=solver_feastol,
                          verbose=verbose_)

        elif solver.lower() == "mosek":
            # https://docs.mosek.com/latest/pythonapi/parameters.html
            MOSEK_OPTS = {'MSK_IPAR_NUM_THREADS': 1,
                        'MSK_IPAR_PRESOLVE_LINDEP_ABS_WORK_TRH': 10}
            try:
                problem.solve(cp.MOSEK, mosek_params=MOSEK_OPTS, eps=solver_reltol, bfs=True, verbose=verbose_)
            except cp.error.SolverError as ex:
                print("Ensure a valid MOSEK license is available: https://docs.mosek.com/latest/licensing/quickstart.html")
                raise

        elif solver.lower() == "pdlp":
            try:
                problem.solve(solver=cp.PDLP, verbose=verbose_)
            except cp.error.SolverError as ex:
                print("Ensure a supported version of ortools is installed for cvxpy")
                raise

        p_stat = problem.status
        if p_stat is None or p_stat in ['infeasible','unbounded'] or problem.variables()[0].value is None:
            raise cp.error.SolverError(f'\nFailed to obtain optimal solution.\
                \nProblem status: {problem.status}.\n')

        lp_sol = problem.variables()[0].value
        sol_rr  = self.run_rr(lp_sol, N, loci_scores, budget, gam, sumsq_penalty=sumsq_penalty)

        head_ = self.head
        i = 0
        while head_ is not None:
            head_.accessible = sol_rr[i]
            head_ = head_.right
            i += 1

        return problem
