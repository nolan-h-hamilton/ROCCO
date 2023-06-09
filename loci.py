r"""
Loci object: ROCCO creates a Locus object (see `locus.py`) for each
$\ell_i, \forall i = 1 \ldots n$. These Locus objects are then linked
together as a Loci object, defined in this code.

The necessary functionality to construct the signal matrix,
$\mathbf{S}_{chr}$, is implemented by way of the Loci object.
"""
import copy
import math
import random
import numpy as np
import cvxpy as cp
from scipy import stats
from locus import Locus


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
        Builds the signal matrix $\mathbf{S}_{chr}$ as defined in the paper.

        Each Locus object's `sig_data` attribute is a $K \times 1$ column
        vector, so $\mathbf{S}_{chr}$ is a $K \times n$ matrix.

        Returns:
            None
        """
        sig_mat = []
        head_ = self.head
        while head_ is not None:
            sig_mat.append(head_.sig_data)
            head_ = head_.right
        return np.array(sig_mat)

    def g3_vec(self, metric_vec) -> np.ndarray:
        r"""
        Compute locus score term $g_3$ for each locus.

        Args:
            metric_vec (np.ndarray): (by default) a $1 \times n$ vector of
                median signal values at each locus.

        Returns:
            np.ndarray: A vector $g_3(i), \forall i=1 \ldots n$
        """
        v_vec = []
        for i,Loc in enumerate(metric_vec):
            if i == 0:
                vel = (abs(metric_vec[i] - metric_vec[i+1]))
            if i == len(metric_vec)-1:
                vel = abs(metric_vec[i] - metric_vec[i-1])
            else:
                vel = max(abs(metric_vec[i] - metric_vec[i+1]),
                      abs(metric_vec[i] - metric_vec[i-1]))

            vel /= metric_vec[i]+1
            v_vec.append(vel)
        return np.array(v_vec)

    def score_loci(self, tau: float = 0.0, c1: float = 1.0, c2: float = 1.0, c3: float = 1.0) -> np.ndarray:
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
        sig_mat = self.get_sig_mat()
        med_vec = np.median(sig_mat,axis=1) # g_1
        mad_vec = stats.median_abs_deviation(sig_mat,axis=1) # g_2
        vel_vec = self.g3_vec(med_vec) # g_3
        s_vec = c1*med_vec - c2*mad_vec + c3*vel_vec

        for i, val in enumerate(s_vec):
            if med_vec[i] <= tau:
                s_vec[i] = 0
        return s_vec

    def run_rr(self, lp_sol, N, loci_scores, budget, gam, eps=1e-5) -> np.ndarray:
        r"""
        Carry out the $\texttt{RR}$ rounding procedure to return a
        'good' integral solution with reference to the LP solution

        Args:
            lp_sol (numpy.ndarray): The LP solution
            N (int): Number of RR solutions to evaluate. If $N \leq 0$, the
                $\texttt{floor\_eps}$ rounding protocol is applied in which
                any decision variable satisfying $\ell_i + \epsilon < 1$
                (`lp_sol[i] + eps < 1`) is rounded down to zero.
            loci_scores (list): mathcal{S} for each locus
            budget (float): budget parameter
            gam (float): gamma
            eps (float): epsilon value for "floor" rounding (Default: eps=1e-5)

    Returns:
        np.ndarray: the integral $\texttt{RR}$ solution, OR $\texttt{floor\_eps}$ sol if $N \leq 0$.
    """
        n = len(loci_scores)

        if N <= 0:
            print("loci.run_rr(): N <= 0 --> returning floor_eps solution")
            return np.floor(lp_sol[:n] + eps)

        rr_sol = None
        best_score = 1e6
        for j in range(N):
            ell_rand_n = []
            for i, ell_i in enumerate(lp_sol[:n]):
                if random.random() <= ell_i:
                    ell_rand_n.append(1)
                else:
                    ell_rand_n.append(0)

            ell_rand_n = np.array(ell_rand_n)
            score = (-loci_scores@ell_rand_n
                       + gam*np.sum(np.abs(np.diff(ell_rand_n,1))))
            is_feas = (np.sum(ell_rand_n) <= math.floor(n*budget))
            if is_feas and score < best_score:
                rr_sol = ell_rand_n
                best_score = score

        if rr_sol is None:
            print("loci.run_rr(): no good solution found via RR --> returning floor solution")
            return np.floor(lp_sol[:n] + eps)

        return np.array(rr_sol)


    def rocco_lp(self, budget: float = .035, tau: float = 0, gam: float = 1, c1: float = 1,
             c2: float = 1, c3: float = 1, verbose_: bool = False, solver: str = "ECOS",
             N: int = 50) -> cp.Problem:
        r"""
        Solve relaxed problem as an LP and assign binary accessibility
        prediction to each locus using the $\texttt{RR}$ or $\texttt{floor\_eps}$
        procedure.

        Args:
            budget (float): budget constraint
            tau (float): tau in $\mathcal{S}(i)$
            gam (float): gamma
            c1 (float): c1 value in $\mathcal{S}(i)$
            c2 (float): c2 value in $\mathcal{S}(i)$
            c3 (float): c3 value in $\mathcal{S}(i)$
            N (int): RR procedure iterations. If N <= 0,
                \texttt{floor\_eps} procedure is applied...
                See `Loci.run_rr()`
            verbose_ (bool): Verbosity flag for the solver
            solver (str): the solver to use: either "ECOS" or "MOSEK"

        Returns:
            cp.Problem: a CVXPY problem object
    """
        loci_scores = self.score_loci(tau, c1, c2, c3) # S_i, i=1...n
        n = len(loci_scores)
        # define problem in CVXPY
        ell = cp.Variable((n,1)) # decision variables
        z = cp.Variable((n-1,1)) # auxiliary variables
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n),
                          z >= cp.diff(ell,1),
                          z >= -1*cp.diff(ell,1)]
        problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(z)),
                             constraints)
        if solver == "ECOS":
            problem.solve(solver=cp.ECOS,verbose=verbose_,max_iters=10000)

        if solver == "MOSEK":
            try:
                problem.solve(cp.MOSEK, verbose=verbose_)
            except Exception as ex:
                print("Ensure a valid `MOSEK` license file is\
                available in your home directory and that the\
                mosek-python interface is installed, e.g.,\
                via 'pip install mosek'")
                raise ex

        lp_sol = np.hstack([[x.value[0] for x in problem.variables()[0]],
                            [x.value[0] for x in problem.variables()[1]]])
        sol_rr  = self.run_rr(lp_sol, N, loci_scores, budget, gam)

        head_ = self.head
        i = 0
        while head_ is not None:
            head_.accessible = sol_rr[i]
            head_ = head_.right
            i += 1

        return problem

    def rocco_ip(self, budget: float = .035, tau: float = 0, gam: float = 1, c1: float = 1,
             c2: float = 1, c3: float = 1, verbose_: bool = False, solver: str = "ECOS_BB") -> cp.Problem:
        r"""
        Solve the unrelaxed optimization problem as an integer program.

        Args:
            budget (float): budget constraint
            tau (float): tau in $\mathcal{S}(i)$
            gam (float): gamma
            c1 (float): c1 value in $\mathcal{S}(i)$
            c2 (float): c2 value in $\mathcal{S}(i)$
            c3 (float): c3 value in $\mathcal{S}(i)$
            verbose_ (bool): Verbosity flag for the solver
            solver (str): the solver to use: either "ECOS" or "MOSEK"

        Returns:
            cp.Problem: a CVXPY problem object

        Notes:
            - There is no efficiency guarantee for this version of
                the problem. Solving may be intractable in some
                cases. MOSEK, or another commerical-grade
                integer optimization package is recommended.
    """
        loci_scores = self.score_loci(tau, c1, c2, c3)
        n = len(loci_scores)
        # define problem in CVXPY
        ell = cp.Variable((n,1),integer=True)
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n)]
        problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(cp.abs(cp.diff(ell,1)))),
                             constraints)

        if solver == "ECOS" or solver =="ECOS_BB":
            problem.solve(cp.ECOS_BB, verbose=verbose_,max_iters=500,
                          feastol=1e-5, abstol=1e-5,reltol=1e-5)
        if solver == "MOSEK":
            try:
                problem.solve(cp.MOSEK, verbose=verbose_)
            except Exception as ex:
                print("Ensure a valid `MOSEK` license file is\
                available in your home directory and that the\
                mosek-python interface is installed, e.g.,\
                via 'pip install mosek'")
                raise ex

        head_ = self.head
        ip_sol = [x.value[0] for x in problem.variables()[0]]
        i = 0
        while head_ is not None:
            head_.accessible = round(ip_sol[i])
            head_ = head_.right
            i += 1


        return problem
