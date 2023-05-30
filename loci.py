"""defines Loci object and associated functions/attributes"""
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
        for a sequence of contiguous loci, all of which
        satisfy `Locus.accessible = 1` and merges them
        into a single locus object (peak). It defines
        the created locus's signal vector with the sum
        of enrichment values among loci in the contiguous
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
        # obtained as the mean signal at each locus
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

            
    def as_list(self):
        """Get a list representation of the Loci object

        Returns:
            An ordered list of all locus objects
        """
        head_ = self.head
        l = []
        while head_ is not None:
            l.append(copy.copy(head_))
            head_ = head_.right
        return l

    
    def combine_selected(self):
        """Find ajacent accessible regions and combine them

        This function identifies contiguous selected regions
        and calls `combine()` to merge them.
        """
        loci_list = self.as_list()
        i = 0
        while i < len(loci_list)-1:
            j = 0
            to_combine = []
            while loci_list[i+j].accessible > 0 and i+j < len(loci_list)-1:
                to_combine.append(loci_list[i+j])
                j += 1
                
            if len(to_combine) > 0:
                self.combine(to_combine[0],to_combine[-1])
                
            i += len(to_combine)+1

            
    def get_sig_mat(self):
        """Build $\mathbf{S}_{chr}$ as defined in the paper
        
        Each locus's sig_data propety is a $K \times 1$ column
        vector, so $\mathbf{S}_{chr}$ is a $K \times n$ matrix.
        """
        sig_mat = []
        head_ = self.head
        while head_ is not None:
            sig_mat.append(head_.sig_data)
            head_ = head_.right
        return np.array(sig_mat)

    
    def get_vel_vec(self, metric_vec):
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
        """Compute scores at each locus using $\mathcal{S}(i)$"""
        sig_mat = self.get_sig_mat()
        med_vec = np.median(sig_mat,axis=1) # g_1
        mad_vec = stats.median_abs_deviation(sig_mat,axis=1) # g_2
        vel_vec = self.get_vel_vec(med_vec) # g_3
        s_vec = c1*med_vec - c2*mad_vec + c3*vel_vec
        
        for i, val in enumerate(s_vec):
            if med_vec[i] <= tau:
                s_vec[i] = 0
        return s_vec

    def run_rr(self, lp_sol, N, loci_scores, budget, gam, eps=1e-5):
        """Carry out the RR rounding procedure"""
        n = len(loci_scores)
        
        if N <= 0:
            return [math.floor(x+eps) for x in lp_sol[:n]]
        
        rr_sol = None
        best_score = 0
        for j in range(N):
            ell_rand_n = []
            ell_rand_aux = []
            for i, ell_i in enumerate(lp_sol):
                if random.random() <= ell_i:
                    if i >= n:
                        ell_rand_aux.append(1)
                    else:
                        ell_rand_n.append(1)
                else:
                    if i >= n:
                        ell_rand_aux.append(0)
                    else:
                        ell_rand_n.append(0)

            ell_rand_n = np.array(ell_rand_n)
            ell_rand_aux = np.array(ell_rand_aux)
            score = (-loci_scores@ell_rand_n
                       + gam*np.sum(np.abs(np.diff(ell_rand_n,1))))
            is_feas = (np.sum(ell_rand_n) <= math.floor(n*budget))
            if is_feas and score < best_score:
                rr_sol = np.concatenate((ell_rand_n,ell_rand_aux),
                                        axis=None)
                best_score = score
                
        if rr_sol is None:
            print("loci.run_rr(): returning floor solution")
            return [math.floor(x+eps) for x in lp_sol[:n]]

        return rr_sol[:n]


    def rocco_lp(self, budget=.035, tau=1, gam=1,
                 c1=1, c2=1, c3=1,
                 verbose_=False,
                 solver="ECOS",
                 N=50):
        """Solve LP as defined in the paper""" 
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
            problem.solve(solver=cp.ECOS,verbose=verbose_,max_iters=1000)
            
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
        """Solve integer program (unrelaxed version) as defined in paper.

        It is recommended to use a commerical grade solver, e.g.,
        MOSEK if solving this unrelaxed version of the problem
        """
        loci_scores = self.score_loci(tau, c1, c2, c3) # S_i, i=1...n
        n = len(loci_scores)
        # define problem in CVXPY
        ell = cp.Variable((n,1),integer=True) # decision variables-integers
        z = cp.Variable((n-1,1)) # auxiliary variables
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n),
                          z >= cp.diff(ell,1),
                          z >= -1*cp.diff(ell,1)]
        problem = cp.Problem(cp.Minimize(-loci_scores@ell + gam*cp.sum(z)),
                             constraints)

        # this may be very slow using the open-source default solver, ECOS
        # consider using the MOSEK solver instead for integral problems
        if solver == "ECOS_BB":
            problem.solve(cp.ECOS_BB, verbose=verbose_,max_iters=500)
        if solver == "MOSEK":
            try:
                problem.solve(cp.MOSEK, verbose=verbose_)
            except Exception as ex:
                print("Ensure a valid `MOSEK` license file is\
                available in your home directory and that the\
                mosek-python interface is installed, e.g.,\
                via 'pip install mosek'")
                raise ex

        dv_vals = []
        for dec_var in problem.variables():
            for dvar_ in dec_var:
                dv_vals.append(dvar_.value[0])
        head_ = self.head
        i = 0

        while head_ is not None:
            head_.accessible = round(dv_vals[i])
            head_ = head_.right
            i += 1
            
        return problem
