"""
A Loci object is a doubly-linked list of Locus objects with some functions
and attributes useful for running ROCCO.

ROCCO creates a Locus object for each entry in the input
`.wig` files and links them together using a Loci object. 

ROCCO_chrom.py defines the signal data for the i-th Locus object
as a list comprised of the signal values at the i-th entries of the
input `.wig` files. In this sense, the Loci object for a given set of
samples corresponds to the signal matrix mathbf{S}_{mathscr{L}} defined
in the paper.
"""
import copy
import math
import random
import numpy as np
import cvxpy as cp
from scipy import stats
from locus import Locus


class LociIterator:
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
    """A collection of locus objects, doubly-linked"""
    def __init__(self, head=None, tail=None):
        self.head = head
        self.tail = tail

    def __iter__(self):
        return LociIterator(self.head)

    def is_empty(self):
        if self.head is None and self.tail is None:
            return True
        return False

    def append_locus(self,locus):
        """Add locus as the new tail of the Loci object"""
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

        # new signal data is a $K \times 1$ vector
        # obtained as the sum signal at each locus
        # comprising the contiguous sequence.
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
        """Get a list representation of the Loci object

        Returns:
            A list representation of the Loci object: i.e.,
              a list with Locus elements ordered as in the
              Loci object.
        """
        head_ = self.head
        list_repr = []
        while head_ is not None:
            list_repr.append(copy.copy(head_))
            head_ = head_.right
        return list_repr

    def combine_selected(self):
        """Find ajacent accessible regions and combine them

        This function identifies contiguous selected regions
        and calls `combine()` to merge them.
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

    def get_sig_mat(self):
        """Build mathbf{S}_{chr} as defined in the paper

        Each locus's sig_data propety is a K times 1 column
        vector, so mathbf{S}_{chr} is a K times n matrix.
        """
        sig_mat = []
        head_ = self.head
        while head_ is not None:
            sig_mat.append(head_.sig_data)
            head_ = head_.right
        return np.array(sig_mat)

    def g3_vec(self, metric_vec):
        """Compute a list of g_3 values for each locus

        Returns:
            a locus-ordered $1 \times n$ numpy array of $g_3$ values
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

    def score_loci(self,tau=1,c1=1,c2=1,c3=1):
        """Compute scores at each locus using mathcal{S}(i)"""
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
        """
        Carry out the RR rounding procedure to find a good integral solution.

        Note: if `N<=0`, the `floor_eps` protocol is applied
        in which any decision variable satisfying `l_i + eps < 1`
        is rounded down to zero.

        Args:
            lp_sol (numpy.ndarray): The LP solution
            N (int): Number of RR solutions to evaluate
            loci_scores (list): mathcal{S} for each locus
            budget (float): budget
            gam (float): gamma
            eps (float, optional): epsilon value for "floor" rounding (Default: eps=1e-5)

    Returns:
        list: an integral vector of decision variables for each locus.

    """
        n = len(loci_scores)

        if N <= 0:
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
            print("loci.run_rr(): returning floor solution")
            return np.floor(lp_sol[:n] + eps)

        return np.array(rr_sol)


    def rocco_lp(self, budget=.035, tau=1, gam=1, c1=1, c2=1, c3=1,
                 verbose_=False, solver="ECOS", N=50) -> cp.Problem:
        """
        Solve LP as defined in the paper and assign accessibility
        prediction tag `0,1` to each locus.

        Args:
            budget (float): budget constraint
            tau (float): tau in mathcal{S}(i)
            gam (float): gamma in LP
            c1 (float): c1 value in mathcal{S}(i)
            c2 (float): c2 value in mathcal{S}(i)
            c3 (float): c3 value in mathcal{S}(i)
            N (int): RR procedure iters. `N<=0` --> floor_eps procedure applied
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

    def rocco_ip(self, budget=.10, tau=1, gam=1,
                 c1=1, c2=1, c3=1,
                 verbose_=False,
                 solver='ECOS_BB'):
        """
        Solve integer program (unrelaxed version) as defined in paper.

        more precise, but slow--It is recommended to use a commerical 
        grade solver, e.g., MOSEK if solving this unrelaxed version
        of the problem

        Parameters and return value same as `Loci.rocco_lp()`
        """

        loci_scores = self.score_loci(tau, c1, c2, c3) # S_i, i=1...n
        n = len(loci_scores)
        # define problem in CVXPY
        ell = cp.Variable((n,1),integer=True) # decision variables-integers
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n)]
        problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(cp.abs(cp.diff(ell,1)))),
                             constraints)

        if solver == "ECOS" or solver =="ECOS_BB":
            problem.solve(cp.ECOS_BB, verbose=verbose_,max_iters=100,
                          feastol=1e-4, abstol=1e-4,reltol=1e-4)
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
