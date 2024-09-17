# -*- coding: utf-8 -*-
r"""
ROCCO: (R)obust (O)pen (C)hromatin Detection via (C)onvex (O)ptimization
==================================================================================

.. image:: logo.png
   :width: 400px
   :align: center

What
----

ROCCO is an algorithm for efficient identification of "consensus peaks" in multiple HTS data samples (namely, ATAC-seq), where read densities are consistently enriched across samples or particularly strong enrichment is observed in a nontrivial subset of samples.

Example Behavior
~~~~~~~~~~~~~~~~~~

In the image below, ROCCO is run on a set of ten heterogeneous ATAC-seq samples (lymphoblast) from independent donors (ENCODE). The samples' tracks are colored gray.

* ROCCO consensus peaks (default parameters) are shown in blue
* MACS2 (pooled, `q=.01`) consensus peaks are shown in red.
* ENCODE cCREs are included as a rough reference of potentially active regions, but note that these regions are not specific to the data samples used in this analysis, nor are they derived from the same cell type or assay.

.. image:: example_behavior.png
   :width: 600px
   :align: center

How
---

ROCCO models consensus peak calling as a constrained optimization problem with an upper-bound on the total proportion of the genome selected as open/accessible and a 'total variation' or 'fragmentation' penalty to promote spatial consistency in active regions and sparsity elsewhere.

Why
---

ROCCO offers several attractive features:

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals
2. **Scaling to large sample sizes** with an asymptotic time complexity independent of sample size
3. **No required training data** or a heuristically determined set of initial candidate peak regions
4. **No rigid thresholds** on the minimum number/width of supporting samples/replicates
5. **Mathematically tractable model** permitting worst-case analysis of runtime and performance


Paper/Citation
--------------

If using ROCCO in your research, please cite the `original paper <https://doi.org/10.1093/bioinformatics/btad725>`_ in *Bioinformatics* (DOI: `btad725`)

.. code-block:: text

   Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
   Bioinformatics, Volume 39, Issue 12, December 2023

Documentation
---------------

ROCCO's documentation is available at `https://nolan-h-hamilton.github.io/ROCCO/ <https://nolan-h-hamilton.github.io/ROCCO/>`_


Installation:
-------------

**PyPI (pip)**:

To install ROCCO via PyPI/pip, you can run one of the following commands:

.. code-block:: bash

    pip install rocco --upgrade # most recent version

.. code-block:: bash

    pip install --pre rocco # most recent release candidate (if available)


**Build ROCCO from Source**

To build ROCCO from source:

1. Clone or download the repository:

.. code-block:: bash

    git clone https://github.com/nolan-h-hamilton/ROCCO.git
    cd ROCCO


2. Build and install the package:

.. code-block:: bash

    python setup.py sdist bdist_wheel
    pip install -e .

System-Level Dependencies:
----------------------------

ROCCO utilizes the popular bioinformatics software Samtools (http://www.htslib.org) and bedtools (https://bedtools.readthedocs.io/en/latest/).

If these are not already available, you can install system-wide with a package manager e.g., 

For Homebrew (MacOS):

.. code-block:: bash

    brew install samtools
    brew install bedtools


For Ubuntu/Debian:

.. code-block:: bash

    sudo apt-get install samtools
    sudo apt-get install bedtools

Conda:

.. code-block:: bash

    conda install bioconda::bedtools
    conda install bioconda::samtools


Both `bedtools` and `samtools` can easily be built from source, too (See  Samtools (http://www.htslib.org) and bedtools (https://bedtools.readthedocs.io/en/latest/).)


Input/Output
-------------

* Input: BAM alignments or BigWig tracks from multiple data samples
   * If BigWig tracks are used as input, no preprocessing can be performed at the alignment level.

* Output: BED file of consensus peak regions


Usage
------

**Run with defaults using BAM input files**:

.. code-block:: bash

        rocco -i sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam -g hg38

**Run on a subset of chromosomes**:
 
* Can be used for quicker debugging
 
.. code-block:: bash

        rocco -i sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam -g hg38 --chroms chr21 chr22

**Run with parametric-sigmoid  transformation of scores**:

* Useful to promote integrality in the LP relaxation and create separation between lower/higher scores.

* See :func:`parsig`

.. code-block:: bash

        rocco -i sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam -g hg38 --use_parsig

**Run ROCCO on 'local-ratio' transformed data**:

* See :func:`rocco.readtracks.apply_transformation`

.. code-block:: bash

        rocco -i sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam -g hg38 --transform_local_ratio

**Run ROCCO with dynamic fragmentation penalties, γ**:

* See :func:`solve_relaxation_chrom_pdlp` or :func:`solve_relaxation_chrom_glop`.

.. code-block:: bash

        rocco -i sample1.bam sample2.bam sample3.bam sample4.bam sample5.bam -g hg38 --scale_gamma

In the above examples, BigWig files can be used too, but as mentioned, no preprocessing can be performed on this data at the *alignment* level.

Usage: Visualization
=====================

See below for a visualization of the effects of several of ROCCO's fundamental options for preprocessing, scoring, optimization, etc.

.. image:: rocco_options.png
   :width: 600px
   :align: center

.. note::

    Using the module(s) themselves programmatically instead of the command-line interface will allow greater flexibility in terms of input/output, data preprocessing/postprocessing, and parallelization but will require more coding on the user's end.

Documentation is provided for the primary :mod:`rocco.rocco` module directly below.

"""

#!/usr/bin/env python
import argparse
import copy
import logging
import multiprocessing
import os
import random
import uuid
from pprint import pformat
from typing import Tuple
import sys

import numpy as np
import pandas as pd
import pysam
import pybedtools
import pyBigWig as pbw

import scipy.stats as stats
import scipy.signal as signal

from google.protobuf import text_format
from ortools.linear_solver import pywraplp
import ortools.glop.parameters_pb2 as parameters_pb2
from ortools.pdlp import solvers_pb2

from rocco.constants import GENOME_DICT
from rocco.readtracks import *

# set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(filename)s: %(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING,
                    format='%(filename)s: %(asctime)s - %(levelname)s - %(message)s')


def _objective_function(sol:np.ndarray, scores:np.ndarray, gamma:float) -> float:
    return (-scores@sol
                + gamma*np.sum(np.abs(np.diff(sol,1))))


def _check_integrality(sol:np.ndarray, int_tol:float=1e-6) -> bool:
    return np.all(np.isin(sol, [int_tol, 1-int_tol]))


def _get_input_type(input_file:str) -> str:
    r"""Determine if `input_file` is a BAM or BigWig file

    The file type is determined by the file extension: '.bam', '.bw', etc. and is not robust
    to incorrectly labelled files.

    :raises ValueError: If file extension is not supported
    :return: a string (extension) representing the file type
    :rtype: str

    """
    file_type = None
    file_ext = str(os.path.splitext(input_file.lower())[1][1:]).lower()

    if file_ext in ['bam']:
        file_type = 'bam'
    elif file_ext in ['bw','bigwig']:
        file_type = 'bw'
    elif file_ext in ['bed', 'bedgraph', 'bg', 'wig', 'wiggle']:
        bedgraph_notice = (
            """\nBedGraph (or other general, non-binary, wiggle-like) input is not supported: Input file must be a BAM alignment file or BigWig file. Please convert to BigWig or use the original BAM files as input. BedGraph/Wiggle files can be converted to BigWig files using UCSC's binary utilities `bedGraphToBigWig`, `wigToBigWig`, respectively, available at <http://hgdownload.soe.ucsc.edu/admin/exe/> or Conda/Mamba <https://anaconda.org/bioconda/ucsc-wigtobigwig>, <https://anaconda.org/bioconda/ucsc-bedgraphtobigwig>.\n""")
        raise ValueError(bedgraph_notice)
    if file_type is None:
        raise ValueError('Input file must be a BAM alignment file or BigWig file')

    return file_type


def minmax_scale_scores(vector:np.ndarray, min_val:float, max_val:float) -> np.ndarray:
    r"""Scale a vector (scores) to the range [min_val,max_val]

    Note that even if data is already in the range [min_val,max_val], this function will still scale the data
    such that the minimum value is `min_val` and the maximum value is `max_val`.

    :param vector: Vector to scale
    :type vector: np.ndarray
    :param min_val: minimum value of the scaled vector (should be near-zero for gamma-scaling)
    :type min_val: float
    :param max_val: maximum value of the scaled vector
    :type max_val: float

    :return: Scaled vector
    :rtype: np.ndarray

    """
    if min_val > max_val:
        raise ValueError('`min_val` must be less than `max_val`')
    if min_val == max_val:
        return np.ones_like(vector)*min_val
    return min_val + (max_val - min_val)*(vector - np.min(vector))/(np.max(vector) - np.min(vector))


def score_central_tendency_chrom(chrom_matrix, method='quantile', quantile=0.50, tprop=0.05, power=1.0) -> np.ndarray:
    r"""Calculate the central tendency of read measures across samples (rows) at each genomic position (column) in a given chrom.
    
    :param chrom_matrix: Matrix of values for a given chromosome
    :type chrom_matrix: np.ndarray
    :param method: Method to calculate central tendency. Options are 'quantile', 'tmean' (trimmed-mean), and 'mean'. Default is quantile, q=.50 (median)
    :type method: str
    :param quantile: Quantile to use if `method` is 'quantile'
    :type quantile: float
    :param tprop: Proportion of values to trim if `method` is 'tmean'
    :type tprop: float
    :param power: Power to raise the central tendency value to
    :type power: float
    :return: Central tendency score
    :rtype: np.ndarray
    
    """
    if get_shape(chrom_matrix)[0] == 1:
        return chrom_matrix[0,:]**power

    method_ = clean_string(method)
    central_tendency_stats = None

    if method_ == 'quantile':
        central_tendency_stats = np.quantile(chrom_matrix, quantile, axis=0, method='nearest')

    if method_ == "tmean":
        # Calculate lower and upper limits for trimming
        llim = np.quantile(chrom_matrix, tprop, axis=0, method='nearest')
        ulim = np.quantile(chrom_matrix, 1-tprop, axis=0, method='nearest')
        
        # Apply trimmed mean column-wise
        central_tendency_stats = np.array([
            stats.tmean(chrom_matrix[:, i], limits=(llim[i], ulim[i]), inclusive=(True, True))
            for i in range(chrom_matrix.shape[1])
        ])

    if method_ == "mean":
        central_tendency_stats = np.mean(chrom_matrix, axis=0)

    if central_tendency_stats is None:
        raise ValueError(f"Central tendency method not recognized: {method}")

    return central_tendency_stats**power


def score_dispersion_chrom(chrom_matrix:np.ndarray, method:str='mad', rng=(25,75), tprop=.05, power:float=1.0) -> np.ndarray:
    r"""Calculate the dispersion of read measures across samples (rows) at each genomic position (column) in a given chrom.
    
    :param chrom_matrix: Matrix of values for a given chromosome
    :type chrom_matrix: np.ndarray
    :param method: Method to calculate dispersion. Options are 'mad', 'iqr', 'tstd' (trimmed std.), and 'std'. Default is 'mad'
    :type method: str
    :param tprop: Proportion of values to trim if `method` is 'tstd'
    :type tprop: float
    :param rng: Range of quantiles to use if `method` is 'iqr'
    :type rng: Tuple[int,int]
    :param power: Power to raise the dispersion value to
    :type power: float
    :return: Dispersion scores
    :rtype: np.ndarray
    
    """

    if get_shape(chrom_matrix)[0] == 1:
        return np.zeros_like(chrom_matrix)[0,:]**power
    
    method_ = clean_string(method)
    dispersion_stats = None

    if method_ == 'mad':
        dispersion_stats = stats.median_abs_deviation(chrom_matrix, axis=0)
    if method_ == 'iqr':
        return stats.iqr(chrom_matrix, rng=rng, axis=0)
    if method_ == 'std':
        dispersion_stats = np.std(chrom_matrix, axis=0)
    if method_ == "tstd":
        llim = np.quantile(chrom_matrix, tprop, method='nearest', axis=0)
        ulim = np.quantile(chrom_matrix, 1-tprop, method='nearest', axis=0)
        dispersion_stats = stats.tstd(chrom_matrix, limits=(llim, ulim), inclusive=(True, True), axis=0)

    if dispersion_stats is None:
        raise ValueError(f"Dispersion method not recognized or could not successfully execute: {method}")

    return dispersion_stats**power


def score_boundary_chrom(signal_vector: np.ndarray, denom:float=1.0, power:float=1.0) -> np.ndarray:
    r"""Calculate the boundary stats for the chromosome score (maximum absolute difference in either direction)
    
    :param signal_vector: (Assumed: central tendency stats) vector for a given chromosome
    :type signal_vector: np.ndarray
    :param denom: Denominator for boundary stats
    :type denom: float
    :param power: Power to raise the boundary stats to
    :type power: float
    :return: Boundary scores
    :rtype: np.ndarray
    
    """
    boundary_stats = np.zeros_like(signal_vector, dtype=float)
    for i in range(len(signal_vector)):
        if i == 0:
            boundary_stats[i] = abs(signal_vector[i] - signal_vector[i+1])/(denom + abs(signal_vector[i]))
        elif i == len(signal_vector) - 1:
            boundary_stats[i] = abs(signal_vector[i] - signal_vector[i-1])/(denom + abs(signal_vector[i]))
        else:
            boundary_stats[i] = max(abs(signal_vector[i] - (signal_vector[i-1])), abs(signal_vector[i] - signal_vector[i+1]))/(denom + abs(signal_vector[i]))

    return boundary_stats**power


def score_chrom_linear(central_tendency_vec:np.ndarray, dispersion_vec:np.ndarray, boundary_vec:np.ndarray, gamma=None, c_1=1.0, c_2=-1.0, c_3=1.0, eps_neg=-1.0e-3, parsig_B=None, parsig_M=None, parsig_R=None) -> np.ndarray:
    r"""Return scores :math:`(\mathbf{G}\mathbf{c})^{\top}` where :math:`\mathbf{G}` is the :math:`n \times 3` matrix of central tendency scores, dispersion scores, and boundary scores for a given chromosome and :math:`\mathbf{c}` is the 3D vector of coefficients.
    
    This is the default scoring function for ROCCO, but various alternatives (e.g., log2fc) have successfully been
    applied as well.
    
    :param central_tendency_vec: Central tendency scores for a given chromosome
    :type central_tendency_vec: np.ndarray
    :param dispersion_vec: Dispersion scores for a given chromosome
    :type dispersion_vec: np.ndarray
    :param boundary_vec: Boundary scores for a given chromosome
    :type boundary_vec: np.ndarray
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere. The ideal value of `gamma` is for a given dataset is dependent on the user's preference. Increase to encourage solutions that are 'block-like'.
    :type  gamma: float
    :param c_1: Coefficient for central tendency scores (:math:`c_1 g_1(i),~i=1,\ldots,n` in the paper)
    :type c_1: float
    :param c_2: Coefficient for dispersion scores (:math:`c_2 g_2(i),~i=1,\ldots,n` in the paper)
    :type c_2: float
    :param c_3: Coefficient for boundary scores (:math:`c_3 g_3(i),~i=1,\ldots,n` in the paper)
    :type c_3: float
    :param use_parsig: If True, apply the custom sigmoid transformation to the scores
    :type use_parsig: False
    :param parsig_M: Upper bound (sup.) of scores under the transformation
    :type parsig_M: float
    :param parsig_B: An inflection point occurs in the transformation at `quantile(scores, parsig_B)/2`. For a given `parsig_R`, larger `parsig_B` value will result in more scores being pushed toward the minimum value of the transformed scores.
    :type parsig_B: float
    :param parsig_R: Scales the rate of change around the inflection point at `quantile(scores, parsig_B)/2`. Higher values of `parsig_R` will result in a steeper sigmoid curve, approaching a step function in the limit with two distinct values.
    :type parsig_R: float
    :param eps_neg: Negative epsilon value to add to the scores
    :type eps_neg: float
    :return: chromosome scores
    :rtype: np.ndarray

    """
    if c_1 < 0: 
        logger.info('Central tendency score coefficient is negative. In the default implementation, this may reward weak signals.')
    if c_2 > 0:
        logger.info('Dispersion score coefficient is positive. In the default implementation, this may reward regions with inconsistent signals among samples.')
        
    chrom_scores = (c_1*central_tendency_vec
                    + c_2*dispersion_vec
                    + c_3*boundary_vec)

    if parsig_B is not None or parsig_R is not None or parsig_M is not None:
        chrom_scores = parsig(chrom_scores, gamma, parsig_B=parsig_B, parsig_M=parsig_M, parsig_R=parsig_R)
    chrom_scores += eps_neg
    return chrom_scores


def parsig(scores, gamma=None, parsig_B=None, parsig_M=None, parsig_R=None, scale_quantile:float=0.50) -> np.ndarray:
    r"""'Parametric sigmoid function' that can be used to transform scores, amplifying the reduced costs/(negative dual price in the primal formulation) values of the most appealing loci/dvars. This can be used to further promote integrality.
    
    :param scores: Scores for each genomic position within a given chromosome. Assumed to be scores generated with :func:`score_chrom_linear`, but can work for others as well.
    :type scores: np.ndarray
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere. The ideal value of `gamma` is for a given dataset is dependent on the user's preference. Increase to encourage solutions that are 'block-like'.
    :type gamma: float
    :param parsig_M: Upper bound (sup.) of scores under the transformation
    :type parsig_M: float
    :param parsig_B: An inflection point occurs in the transformation at `quantile(scores, parsig_B)/2`. For a given `parsig_R`, larger `parsig_B` value will result in more scores being pushed toward the minimum value of the transformed scores.
    :type parsig_B: float
    :param parsig_R: Scales the rate of change around the inflection point at `quantile(scores, parsig_B)/2`. Higher values of `parsig_R` will result in a steeper sigmoid curve, approaching a step function in the limit with two distinct values.
    :type parsig_R: float
    :param scale_quantile: Quantile of unique values in initial scores used to extend the new range of scores under the transformation
    :type scale_quantile: float
    :return: transformed scores
    :rtype: np.ndarray

    """

    if gamma is None:
        gamma = 1.0

    if parsig_B is None:
        parsig_B = 0.95
    elif parsig_B < 0 or parsig_B > 1:
        raise ValueError('`parsig_B` must be in the range (0,1)')
    
    if parsig_R is None:
        parsig_R= 2.0
    elif parsig_R < 0:
        raise ValueError('`parsig_R` must be greater than zero')

    scores_ = np.array(scores)
    scores_unique = np.array(sorted(list(set(list(scores)))))
    scale_ = abs(np.quantile(scores_unique, q=scale_quantile, method='nearest'))
    M_ = (2*gamma +1) + scale_ if parsig_M is None else parsig_M
    if M_ < 0:
        raise ValueError('`parsig_M` must be greater than zero')

    B_ = np.quantile(scores_, parsig_B, method='nearest')
    parsig_values = M_ / (1 + np.exp(-(parsig_R * ((scores - (B_ / 2))))))
    parsig_values = parsig_values - np.min(parsig_values) # in case first inflection occurs below zero
    logger.info(f'Parsig transformation applied with M={M_}, B={B_}, R={parsig_R}')
    return parsig_values


def get_warm_idx(scores, budget, gamma, warm_thresh=None) -> Tuple[np.ndarray, np.ndarray, float]:
    r"""Prior to solving, identify 'warm' indices--those with scores greater than the worst-case fragmentation penalty that could be incurred by their selection assuming integrality (`2*gamma`). Note that this loses meaning if gamma is scaled. 

    :param scores: Scores for each genomic position within a given chromosome
    :type scores: np.ndarray
    :param budget: :math:`b` upper bounds the proportion of the chromosome that can be selected as open/accessible
    :type budget: float
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere. The ideal value of `gamma` is for a given dataset is dependent on the user's preference. Increase to encourage solutions that are 'block-like'.
    :type gamma: float
    :param warm_thresh: Threshold for warm indices. If None, the threshold is set to 2*gamma
    :type warm_thresh: float

    :return: Warm indices, warm scores, and the proportion of warm indices to the maximum number of selections
    :rtype: Tuple[np.ndarray, np.ndarray, float]
    
    """
    n = len(scores)
    max_selections = np.floor(n*budget)
    warm_idx = []
    warm_scores = []
    if warm_thresh is None:
        warm_thresh = 2*gamma + 1
    for i in range(n):
        if scores[i] > warm_thresh:
            warm_idx.append(i)
            warm_scores.append(scores[i])
    warm_idx = np.array([int(x) for _, x in sorted(zip(warm_scores, warm_idx), reverse=True, key=lambda pair: pair[0])], dtype=int)
    
    if len(warm_idx) == 0:
        return np.array([]), np.array([]), 0.0

    return warm_idx, scores[warm_idx], len(warm_idx)/max_selections


def solve_relaxation_chrom_pdlp(scores,
                    budget:float=0.035,
                    gamma:float=1.0,
                    beta:float=None,
                    denom_const:float=None,
                    pdlp_proto=None,
                    pdlp_presolve_options_use_glop:bool=True,
                    pdlp_termination_criteria_eps_optimal_absolute:float=1.0e-8,
                    pdlp_termination_criteria_eps_optimal_relative:float=1.0e-8,
                    pdlp_use_low_precision:bool=False,
                    scale_gamma:bool=False,
                    hint=None,
                    threads:int=0,
                    verbose:bool=False,
                    save_model:str=None) -> Tuple[np.ndarray, float]:
    r"""Solve the relaxation for a specific chromosome using the *first-order* method, pdlp
    
    See the  `full paper for pdlp <https://proceedings.neurips.cc/paper/2021/hash/a8fbbd3b11424ce032ba813493d95ad7-Abstract.html>`_ for a complete technical exposition.
    
    `OR-tools linear programming resources and documentation <https://developers.google.com/optimization/lp>`_
    
    :param scores: Scores for each genomic pisition within a given chromosome
    :type scores: np.ndarray
    :param budget: :math:`b` upper bounds the proportion of the chromosome that can be selected as open/accessible
    :type budget: float
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere. The ideal value of `gamma` is for a given dataset is dependent on the user's preference. Increase to encourage solutions that are 'block-like'.
    :type  gamma: float
    :param beta: Exponent for the denominator in the fragmentation penalty. If None, defaults to 1/2 when `scale_gamma` is True.
    :type beta: float
    :param denom_const: Constant to add to the denominator in the fragmentation penalty. If None, defaults to 1.0 when `scale_gamma` is True.
    :type denom_const: float
    :param scale_gamma: If `True`, Scale the fragmentation penalty (γ) positionally based on the difference between adjacent scores. This will yield a more nuanced fragmentation penalty that does not discourage differences in adjacent decision variables if their corresponding scores are dissimilar and thus reflect a true change in state that should be reflected in the solution. See `beta` and `denom_const` for more details.
    :type scale_gamma: bool
    :param pdlp_proto: pdlp-specific protocol buffer. If this is not None, the explicit solver arguments in this function definition are ignored. See `<https://protobuf.dev>`_ for more information on protocol buffers and `solvers.proto <https://github.com/google/or-tools/blob/2c333f58a37d7c75d29a58fd772c9b3f94e2ca1c/ortools/pdlp/solvers.proto>`_ for Google's pdlp-specific protocol buffer.
    :type pdlp_proto: solvers_pb2.PrimalDualHybridGradientParams
    :param pdlp_presolve_options_use_glop: Use glop's presolve routines but solve with pdlp. Recommended for most cases unless the user is confident in the problem's structure and the solver's behavior and computational resources are limited.
    :type pdlp_presolve_options_use_glop: bool
    :param pdlp_termination_criteria_eps_optimal_absolute: Appealing to strong duality for LPs, the duality gap (difference between the primal objective function and the dual objective function at a given iteration) must be less than this value *plus a scaling of the primal/dual objective values* (See `solvers.proto` linked above for exact details).  If computational resources are limited, consider using `1.0e-4` per the `ortools <https://developers.google.com/optimization/lp/lp_advanced>`_ documentation.
    :type pdlp_termination_criteria_eps_optimal_absolute: float
    :param pdlp_termination_criteria_eps_optimal_relative: Relative termination criterion for the pdlp solver. If computational resources are limited, consider using `1.0e-4` per the `ortools <https://developers.google.com/optimization/lp/lp_advanced>`_ documentation.
    :type pdlp_termination_criteria_eps_optimal_relative: float
    :param pdlp_use_low_precision: Use "loose solve"/low precision mode. This will override the termination arguments to weakened criteria (1.0e-4).
    :type pdlp_use_low_precision: bool
    :param threads: Number of threads to use for optimization. Default is 0 (use default number of threads defined in `solvers.proto`). If threads is negative, the number of threads used will be the maximum of the number of available CPUs minus 1, or 1 in the case of a single-core machine.
    :type threads: int
    :param verbose: Enable verbose output
    :type verbose: bool
    :param save_model: Save the model as an MPS file
    :type save_model: str
    :return: Solution vector and optimal value
    :rtype: Tuple[np.ndarray, float]
    
    """

    solver = pywraplp.Solver.CreateSolver('PDLP')
    pdlp_parameters = solvers_pb2.PrimalDualHybridGradientParams()
    pdlp_parameters.presolve_options.use_glop = pdlp_presolve_options_use_glop
    pdlp_parameters.termination_criteria.simple_optimality_criteria.eps_optimal_absolute = pdlp_termination_criteria_eps_optimal_absolute
    pdlp_parameters.termination_criteria.simple_optimality_criteria.eps_optimal_relative = pdlp_termination_criteria_eps_optimal_relative
    if pdlp_use_low_precision:
        logger.info(f'Using "loose solve"/low precision mode: {pdlp_use_low_precision}')
        logger.info('Overriding termination arguments to weakened criteria (1.0e-4)')
        pdlp_parameters.termination_criteria.simple_optimality_criteria.eps_optimal_absolute = 1.0e-4
        pdlp_parameters.termination_criteria.simple_optimality_criteria.eps_optimal_relative = 1.0e-4
    if pdlp_proto is None:
        proto = text_format.MessageToString(pdlp_parameters)
    else:
        proto = text_format.MessageToString(pdlp_proto)
    solver.SetSolverSpecificParametersAsString(proto)

    if threads >= 1:
        solver.SetNumThreads(threads)
    elif threads < 0:
        solver.SetNumThreads(max(multiprocessing.cpu_count()-1,1))

    if verbose:
        solver.EnableOutput()

    n = len(scores)
    logger.info(f'Building problem with {n} variables')
    ell = [solver.NumVar(0, 1, f'ell_{i}') for i in range(n)]
    ell_aux = [solver.NumVar(0, 1, f'ell_aux_{j}') for j in range(n-1)]
    logger.info(f'Setting constraints')
    solver.Add(sum(ell) <= np.floor(budget * n)*1.0, name='budget')
    for i in range(n - 1):
        # log every 10% of the way
        if i % (n//10) == 0:
            logger.info(f'Constraint ({i}/{n})')
        solver.Add(ell_aux[i] >= ell[i] - ell[i + 1])
        solver.Add(ell_aux[i] >= -(ell[i] - ell[i + 1]))
    objective = solver.Objective()
    logger.info(f'Setting coefficients in the objective')
    for i in range(n):
        # log every 10% of the way
        if i % (n//10) == 0:
            logger.info(f'Coefficient ({i}/{n})')
        objective.SetCoefficient(ell[i], -1.0*scores[i])
    denom = None
    if scale_gamma:
        if denom_const is None:
            denom_const = 1.0
        if beta is None:
            beta = 1.0/2.0
        logger.info(f'Scaling γ positionally: `γ/(|Score_i - Score_i+1| + ε)^β`, with β={beta} and ε={denom_const}')
        denom = (np.abs(np.diff(scores,1)) + denom_const)**beta
    
    for j in range(n-1):
        # log every 10% of the way
        if j % (n//10) == 0:
            logger.info(f'Coefficient (aux.) ({j}/{n})')
        if scale_gamma:
            objective.SetCoefficient(ell_aux[j], (1.0*gamma)/denom[j])
        elif denom is None:
            objective.SetCoefficient(ell_aux[j], (1.0*gamma))
    objective.SetMinimization()
    
    if hint is not None:
        try:
            solver.SetHint(ell, hint)
        except Exception as e:
            logger.info(f'Could not apply hint solution:\n{e}')

    # threads = 0: use default number of threads
    if threads > 0:
        solver.SetNumThreads(threads)
    elif threads < 0:
        solver.SetNumThreads(max(multiprocessing.cpu_count()-1,1))

    logger.info('Solving...')
    solver.Solve()

    ell_arr = np.zeros(n)
    for i in range(n):
        ell_arr[i] = ell[i].solution_value()
    optimal_value = solver.Objective().Value()

    summary = {
    'iterations': solver.iterations(),
    'time': solver.wall_time(),
    }

    if verbose:
        logger.info(f'Solver summary:\n{summary}')

    if save_model is not None:
        try:
            with open(save_model, 'w') as f:
                f.write(solver.ExportModelAsMpsFormat(fixed_format=True, obfuscated=False))
        except Exception as e:
            logger.info(f'Could not save model:\n{e}')
    return ell_arr, optimal_value


def solve_relaxation_chrom_glop(scores,
                     budget:float=0.035,
                     gamma:float=1.0, 
                     beta:float=None,
                     denom_const:float=None,
                     glop_parameters=None,
                     glop_dual_feasibility_tolerance:float=1.0e-8,
                     glop_primal_feasibility_tolerance:float=1.0e-8,
                     glop_use_scaling:bool=True,
                     glop_initial_basis:str='TRIANGULAR',
                     glop_push_to_vertex: bool=True,
                     glop_use_dual_simplex=None,
                     glop_allow_simplex_algorithm_change:bool=False,
                     scale_gamma:bool=False,
                     hint=None,
                     threads:int=0,
                     verbose:bool=False,
                     save_model:str=None) -> Tuple[np.ndarray, float]:
    r"""Solve the relaxation for a specific chromosome using the *simplex-based* method, glop

    `OR-tools linear programming resources and documentation <https://developers.google.com/optimization/lp>`_
    
    A simplex-based metho that yields *corner-point feasible* (vertex) solutions that are attractive for a variety of technical reasons, especially in ROCCO's setting. In practice, however, for sufficiently strict termination criteria, PDLP (:func:`solve_relaxation_chrom_pdlp`) yields nearly identical solutions to glop and scales better to large problems.
    Particularly after the randomized rounding step, the difference in solutions is negligible. 
    
    :param scores: Scores for each genomic pisition within a given chromosome
    :type scores: np.ndarray
    :param budget: :math:`b` upper bounds the proportion of the chromosome that can be selected as open/accessible
    :type budget: float
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere. The ideal value of `gamma` is for a given dataset is dependent on the user's preference. Increase to encourage solutions that are 'block-like'.
    :type  gamma: float
    :param beta: Exponent for the denominator in the fragmentation penalty. If None, defaults to 1/2 when `scale_gamma` is True.
    :type beta: float
    :param denom_const: Constant to add to the denominator in the fragmentation penalty. If None, defaults to 0.5 when `scale_gamma` is True.
    :type denom_const: float
    :param scale_gamma: If `True`, Scale the fragmentation penalty (γ) positionally based on the difference between adjacent scores. This will yield a more nuanced fragmentation penalty that does not discourage differences in adjacent decision variables if their corresponding scores are dissimilar and thus reflect a true change in state that should be reflected in the solution. See `beta` and `denom_const` for more details.
    :type scale_gamma: bool
    :param glop_parameters: glop-specific protocol buffer. If this is not None, the explicit solver arguments in this function definition are ignored. See `<https://protobuf.dev>`_ for more information on protocol buffers and `parameters.proto <https://github.com/google/or-tools/blob/stable/ortools/glop/parameters.proto>`_ for Google's glop-specific protocol buffer.
    :type glop_parameters: parameters_pb2.GlopParameters
    :param glop_dual_feasibility_tolerance: Dual feasibility tolerance for glop.
    :type glop_dual_feasibility_tolerance: float
    :param glop_primal_feasibility_tolerance: Primal feasibility tolerance for glop.
    :type glop_primal_feasibility_tolerance: float
    :param glop_use_scaling: Use scaling for glop. If `True`, the scaling method is set to `EQUILIBRATION` and the cost scaling is set to `CONTAIN_ONE_COST_SCALING`. Recommended for most cases unless the user is confident the problem is not poorly-scaled. Future releases of ROCCO should support additional options.
    :type glop_use_scaling: bool
    :param glop_initial_basis: Determines which of the glop-supported heuristics to identify an initial basis is executed. Options are `'NONE'`, `'TRIANGULAR'`, and `'MAROS'`. `'TRIANGULAR'` is the default and recommended option.
    :type glop_initial_basis: str
    :param threads: Number of threads to use for optimization. Default is 0 (use default number of threads defined in `solvers.proto`). If threads is negative, the number of threads used will be the maximum of the number of available CPUs minus 1, or 1 in the case of a single-core machine.
    :type threads: int
    :param verbose: Enable verbose output
    :type verbose: bool
    :param save_model: Save the model as an MPS file
    :type save_model: str
    :return: Solution vector and optimal value
    :rtype: Tuple[np.ndarray, float]
    
    """
    solver = pywraplp.Solver.CreateSolver('GLOP')
    if glop_parameters is None:
        glop_parameters = parameters_pb2.GlopParameters()
        glop_parameters.dual_feasibility_tolerance = glop_dual_feasibility_tolerance
        glop_parameters.primal_feasibility_tolerance = glop_primal_feasibility_tolerance
        glop_parameters.push_to_vertex = glop_push_to_vertex
        if glop_initial_basis is not None:
            if glop_initial_basis.upper() == 'NONE':
                glop_parameters.initial_basis = parameters_pb2.GlopParameters.NONE
            if glop_initial_basis.upper() == 'TRIANGULAR':
                glop_parameters.initial_basis = parameters_pb2.GlopParameters.TRIANGULAR
            elif glop_initial_basis.upper() == 'MAROS':
                glop_parameters.initial_basis = parameters_pb2.GlopParameters.MAROS
        glop_parameters.use_scaling = glop_use_scaling
        
        if glop_use_scaling:
            glop_parameters.scaling_method = parameters_pb2.GlopParameters.EQUILIBRATION
            glop_parameters.cost_scaling = parameters_pb2.GlopParameters.CONTAIN_ONE_COST_SCALING
        if glop_use_dual_simplex is not None:
            glop_parameters.use_dual_simplex = glop_use_dual_simplex
        glop_parameters.allow_simplex_algorithm_change = glop_allow_simplex_algorithm_change

    proto = text_format.MessageToString(glop_parameters)
    solver.SetSolverSpecificParametersAsString(proto)

    if verbose:
        logger.info(f'Solver parameters:\n{pformat(glop_parameters)}\n')

    n = len(scores)
    logger.info(f'Building problem with {n} variables')
    ell = [solver.NumVar(0, 1, f'ell_{i}') for i in range(n)]
    ell_aux = [solver.NumVar(0, 1, f'ell_aux_{j}') for j in range(n-1)]
    logger.info(f'Setting constraints')
    solver.Add(sum(ell) <= np.floor(budget * n)*1.0, name='budget')
    for i in range(n - 1):
        # log every 10% of the way
        if i % (n//10) == 0:
            logger.info(f'Constraint ({i}/{n})')
        solver.Add(ell_aux[i] >= ell[i] - ell[i + 1])
        solver.Add(ell_aux[i] >= -(ell[i] - ell[i + 1]))
    objective = solver.Objective()
    logger.info(f'Setting coefficients in the objective')
    for i in range(n):
        # log every 10% of the way
        if i % (n//10) == 0:
            logger.info(f'Coefficient ({i}/{n})')
        objective.SetCoefficient(ell[i], -1.0*scores[i])
    denom = None
    if scale_gamma:
        if denom_const is None:
            denom_const = 1.0
        if beta is None:
            beta = 1.0/2.0
        logger.info(f'Scaling γ positionally: `γ/(|Score_i - Score_i+1| + ε)^β`, with β={beta} and ε={denom_const}')
        denom = (np.abs(np.diff(scores,1)) + denom_const)**beta

    for j in range(n-1):
        # log every 10% of the way
        if j % (n//10) == 0:
            logger.info(f'Coefficient (aux.) ({j}/{n})')
        if scale_gamma:
            objective.SetCoefficient(ell_aux[j], (1.0*gamma)/denom[j])
        elif denom is None:
            objective.SetCoefficient(ell_aux[j], (1.0*gamma))
    objective.SetMinimization()

    if hint is not None:
        try:
            solver.SetHint(ell, hint)
        except Exception as e:
            logger.info(f'Could not apply hint solution:\n{e}')

    # threads = 0: use default number of threads
    if threads > 0:
        solver.SetNumThreads(threads)
    elif threads < 0:
        solver.SetNumThreads(max(multiprocessing.cpu_count()-1,1))

    if verbose:
        solver.EnableOutput()
    logger.info('Solving...')
    solver.Solve()
    ell_arr = np.zeros(n)
    for i in range(n):
        ell_arr[i] = ell[i].solution_value()
    optimal_value = solver.Objective().Value()

    summary = {
    'iterations': solver.iterations(),
    'time': solver.wall_time(),
    }

    if verbose:
        logger.info('Solver summary:')
        logger.info(summary)

    if save_model is not None:
        try:
            with open(save_model, 'w') as f:
                f.write(solver.ExportModelAsMpsFormat(fixed_format= True, obfuscated=False))
        except Exception as e:
            logger.info(f'Could not save model:\n{e}')

    return ell_arr, optimal_value


def get_floor_eps_sol(chrom_lp_sol:np.ndarray, budget:float,
                      int_tol:float=1e-6,
                      eps_mult:float=1.01) -> np.ndarray:
    r"""Compute the `floor_eps` heuristic from the relaxation
    
    Adds a small :math:`\epsilon` to each decision variable before applying the floor function. Addresses solutions in the interior of the feasible region near integer corner point solutions.  The `floor_eps` approach is a crude heuristic to obtain an initial integer-feasible solution from the LP relaxation, but is
    well-suited as an initial reference point for the ROCCO-RR heuristic. In many cases, the relaxed LP solution is
    nearly integral and the `floor_eps` heuristic and randomized `ROCCO-RR` procedure has little effect.
    
    :param chrom_lp_sol: Solution vector from the LP relaxation
    :type chrom_lp_sol: np.ndarray
    :param budget: :math:`b` upper bounds the proportion of the chromosome that can be selected as open/accessible
    :type budget: float
    :param int_tol: If a decision variable is within `int_tol` of 0 or 1, it is considered integral
    :type int_tol: float
    :param eps_mult: Value by which to divide the initial epsilon value
    :type eps_mult: float
    :return: Initial integer-feasible solution
    :rtype: np.ndarray
    
    """
    
    # check if LP solution is already integral
    if _check_integrality(chrom_lp_sol, int_tol):
        return np.round(chrom_lp_sol)

    # Round the others where possible (within int_tol)
    for i in range(len(chrom_lp_sol)):
        if chrom_lp_sol[i] < int_tol:
            chrom_lp_sol[i] = 0
        elif chrom_lp_sol[i] > 1 - int_tol:
            chrom_lp_sol[i] = 1

    if eps_mult <= 1:
        raise ValueError('`eps_mult` must be greater than 1')

    n = len(chrom_lp_sol)
    plus_half = np.array([x for x in chrom_lp_sol if x > 0 and x < 1])
    if plus_half is None or len(plus_half) == 0:
        eps_cpy = 0
    else:
        eps_cpy = (1/2) + 1.0e-4
    init_sol = np.floor(chrom_lp_sol + eps_cpy)
    floor_eps_iter_ct = 0
    while np.sum(init_sol) > np.floor(n*budget):
        eps_cpy = eps_cpy/eps_mult
        init_sol = np.array([np.floor(chrom_lp_sol[i] + eps_cpy) for i in range(len(chrom_lp_sol))])
        floor_eps_iter_ct += 1
    return init_sol


def get_rround_sol(chrom_lp_sol, scores, budget, gamma,
                   rand_iter=1000, int_tol=1.0e-6,
                   eps_mult:float=1.01) -> Tuple[np.ndarray, float]:
    r"""Get the ROCCO-RR solution from the solution to the relaxation
    
    :param chrom_lp_sol: Solution vector from the LP relaxation
    :type chrom_lp_sol: np.ndarray
    :param scores: Scores for each genomic position within a given chromosome
    :type scores: np.ndarray
    :param budget: :math:`b` upper bounds the proportion of the chromosome that can be selected as open/accessible
    :type budget: float
    :param gamma: :math:`\gamma` is the coefficient for the fragmentation penalty used to promote spatial consistency in distinct open genomic regions and sparsity elsewhere.
    :type  gamma: float
    :param rand_iter: Number of randomizations to obtain the ROCCO-RR solution
    :type rand_iter: int
    :param int_tol: If a decision variable is within `int_tol` of 0 or 1, it is considered integral and ignored in the randomization step
    :type int_tol: float
    :param eps_mult: Value by which to divide the initial :func:`get_floor_eps_sol` epsilon value
    :type eps_mult: float
    :return: ROCCO-RR solution and optimal value    

    """
    n = len(scores)
    floor_sol = np.floor(chrom_lp_sol)
    floor_score = _objective_function(sol=floor_sol, scores=scores, gamma=gamma)
    
    floor_eps_sol = get_floor_eps_sol(chrom_lp_sol, budget, int_tol=int_tol, eps_mult=eps_mult)
    floor_eps_score = _objective_function(sol=floor_eps_sol, scores=scores, gamma=gamma)
    init_sol = None
    init_score = None
    if floor_eps_score < floor_score and np.sum(floor_eps_sol) <= np.floor(n*budget):
        init_sol = floor_eps_sol
        init_score = floor_eps_score
    else:
        init_sol = floor_sol
        init_score = floor_score

    if rand_iter <= 0:
        logger.info('`rand_iter < 0`, returning floor-based solution.')
        return init_sol, init_score

    nonint_loci = np.array([i for i in range(len(chrom_lp_sol)) if chrom_lp_sol[i] > int_tol and chrom_lp_sol[i] < 1-int_tol], dtype=int)
    logger.info(f'Number of fractional decision variables in LP (tol={int_tol}): {len(nonint_loci)}')
    if len(nonint_loci) > 0:
        logger.info(f"{nonint_loci}, {chrom_lp_sol[nonint_loci]}\n")
    if len(nonint_loci) == 0:
        logger.info('LP solution is integral. Returning.\n')
        return init_sol, init_score
 
    rround_sol = copy.copy(init_sol)
    rround_sum = np.sum(init_sol)
    best_score = init_score
    for j in range(rand_iter):
        ell_rand_n = copy.copy(chrom_lp_sol)
        # for efficiency, only round `nonint_loci`
        for idx in nonint_loci:
            if random.random() <= chrom_lp_sol[idx]:
                ell_rand_n[idx] = 1
            else:
                ell_rand_n[idx] = 0

        score = _objective_function(sol=ell_rand_n, scores=scores, gamma=gamma)
        sol_sum = np.sum(ell_rand_n)
        is_feas = (sol_sum <= np.floor(n*budget))
        # if two solutions are equal in the objective, choose the sparser of the two
        if (is_feas and score < best_score) or (is_feas and score == best_score and sol_sum < rround_sum):
            rround_sol = ell_rand_n
            rround_sum = sol_sum
            best_score = score

    return rround_sol, best_score


def chrom_solution_to_bed(chromosome, intervals, solution, ID=None,
                          check_gaps_intervals=True, min_length_bp=None)-> str:
    r"""Convert the ROCCO-generated vector of decision variables for a given chromosome to a BED file
    
    :param chromosome: Chromosome name
    :type chromosome: str
    :param intervals: Intervals for the chromosome
    :type intervals: np.ndarray
    :param solution: Solution vector for the chromosome
    :type solution: np.ndarray
    :param ID: Unique identifier for the solution
    :type ID: str
    :param check_gaps_intervals: Check if intervals are contiguous and fixed width
    :type check_gaps_intervals: bool
    :param min_length_bp: Minimum length of a region to be included in the output BED file
    :type min_length_bp: int
    :return: Output BED file
    :rtype: str
    
    """
    if len(intervals) != len(solution):
        raise ValueError(f'Intervals and solution must have the same length at the pre-merge stage: {len(intervals)} != {len(solution)}')

    if check_gaps_intervals:
        if len(set(np.diff(intervals))) > 1:
            raise ValueError(f'Intervals must be contiguous: {set(np.diff(intervals))}')
        
    if ID is None:
        output_file = f'rocco_{chromosome}.bed'
    else:
        output_file = f'rocco_{ID}_{chromosome}.bed'

    with open(output_file, 'w') as f:
        for i in range(len(intervals)-1):
            # At this point, solutions in the default implementation
            # should be binary, but for potential future use in other
            # just use a threshold of 0.50
            if solution[i] > 0.50:
                f.write(f'{chromosome}\t{intervals[i]}\t{intervals[i+1]}\n')
    chrom_pbt = pybedtools.BedTool(output_file).sort().merge()
    # filter out regions less than min_length_bp if specified
    if min_length_bp is not None:
        chrom_pbt = chrom_pbt.filter(lambda x: int(x[2]) - int(x[1]) >= min_length_bp)
    chrom_pbt.saveas(output_file)
    if os.path.exists(output_file):
        return output_file


def combine_chrom_results(chrom_bed_files:list, output_file:str, name_features:bool=False) -> str:
    r"""Combine the results from individual chromosome solutions into a single BED file after running
    ROCCO on each chromosome
    
    :param chrom_bed_files: List of BED files for each chromosome
    :type chrom_bed_files: list
    :param output_file: Output BED file
    :type output_file: str
    :param name_features: Name the features in the output BED file
    :type name_features: bool
    :return: Output BED file
    :rtype: str
    
    """
    printed_colct_msg = False
    
    if os.path.exists(output_file):
        logger.info(f'Removing existing output file: {output_file}')
        try:
            os.remove(output_file)
        except:
            logger.info(f'Could not remove existing output file: {output_file}.')
    with open (output_file, 'w') as f:
        for chrom_bed_file in chrom_bed_files:
            if not os.path.exists(chrom_bed_file):
                raise FileNotFoundError(f'File does not exist: {chrom_bed_file}')
            try:
                chrom_pbt = pybedtools.BedTool(chrom_bed_file).sort().merge()
            except Exception as e:
                logger.info(f'Could not read or merge BED file: {chrom_bed_file}\n{e}\n')
                raise
            if chrom_pbt.field_count() > 3:
                if not printed_colct_msg:
                    logger.info('More than 3 columns detected in the input BED files. Extra columns will be ignored.')
                    printed_colct_msg = True
            for feature_ in chrom_pbt:
                if name_features:
                    feature_name = f'{feature_.chrom}_{feature_.start}_{feature_.stop}'

                    f.write(f'{feature_.chrom}\t{feature_.start}\t{feature_.stop}\t{feature_name}\n')
                else:
                    f.write(f'{feature_.chrom}\t{feature_.start}\t{feature_.stop}\n')
    return output_file


def main():
    ID = str(int(uuid.uuid4().hex[:5], base=16))
    logger.info(f'\nID: {ID}')
    epilog_cli_help = (
        "\nGitHub (Homepage): <https://github.com/nolan-h-hamilton/ROCCO/>\n"
        "Documentation: <https://nolan-h-hamilton.github.io/ROCCO/>\n"
        "Paper: <https://doi.org/10.1093/bioinformatics/btad725>\n"
    )
    parser = argparse.ArgumentParser(description='ROCCO Consensus Peak Detection Algorithm for Multisample HTS Datasets', add_help=True, formatter_class=argparse.RawTextHelpFormatter, epilog=epilog_cli_help)
    parser.add_argument('--input_files', '-i', nargs='+', help='BAM alignment files or BigWig files corresponding to samples')
    parser.add_argument('--output', '--outfile', '-o', type=str, default=f"rocco_peaks_output_{ID}.bed")
    parser.add_argument('--genome', '-g', default=None, help='Genome assembly. Invoking this argument with a supported assembly (hg38, hg19, mm10, mm39, dm6) will use default resources (`--chrom_sizes_file`), parameters (`--params`) and EGS (`--effective_genome_size`) specific to that assembly that come with the ROCCO package by default. If this argument is not invoked, you can just supply the required arguments manually with the command-line arguments.')
    parser.add_argument('--chrom_sizes_file', '-s', default=None, help='Chromosome sizes file. Required if genome is not specified')
    parser.add_argument('--effective_genome_size', type=int, default=None, help='Effective genome size. Required if genome is not specified and using RPGC normalization')
    parser.add_argument('--chroms', nargs='+', type=str, default=[], help='Chromosomes to process. If not specified, all chromosomes will be processed')
    parser.add_argument('--skip_chroms', nargs='+', type=str, default=[], help='Chromosomes to skip')
    parser.add_argument('--verbose', action='store_true', help='Invoke for verbose output')


    # optimization-related arguments
    parser.add_argument('--solver', default='pdlp', choices=['glop', 'pdlp', 'GLOP', 'PDLP'], help='Solver to use for optimization. Default is pdlp (first-order pdhg method) but glop (simplex-based method) is also supported. See module documentation for more information.')
    parser.add_argument('--int_tol', type=float, default=1.0e-6, help='If a decision variable is within `int_tol` of 0 or 1, it is considered integral.')
    parser.add_argument('--eps_mult', type=float, default=1.01, help='Successively divides the floor_eps solution epsilon value until the floor_eps solution, floor(lp_solution + ε) is feasible.')
    parser.add_argument('--rand_iter', type=int, default=1000, help='Number of random iterations for the ROCCO-RR randomization procedure. If less than or equal to zero, ROCCO-RR is not applied and one of the floor-based heuristics is used to obtain an integer feasible solution.')
    parser.add_argument('--budget', type=float, default=None, help='Upper bounds the proportion of the genome that can be selected as open chromatin. NOTE: if invoked, this argument value will override the chromosome-specific budget values in (`--params`), assuming they were provided.')
    parser.add_argument('--gamma', type=float, default=None, help='Gamma penalty in the optimization problem. Controls weight of the "fragmentation penalty" to promote spatial consistency over enriched regions and sparsity elsewhere.  NOTE: if invoked, this value will override the chromosome-specific gamma values (`--params`), assuming they were provided.')
    parser.add_argument('--params', type=str, default=None, help='CSV file containing chromosome-specific optimization parameters. Each supported genome has a custom `--params` file packaged with ROCCO, but a custom file can be easily be created manually and passed to ROCCO with this argument. Consider using this argument to set chromosome-specific budget and gamma values custom to your data and preferences.')
    parser.add_argument('--threads', type=int, default=-1, help='Number of threads to use for optimization. Default is -1 (use all available threads)')
    parser.add_argument('--eps_optimal_absolute', '--atol', type=float, default=1.0e-8, help="pdlp solver only. One component of the bound on the duality gap used to check for convergence. If computational resources are limited, consider using `1.0e-4` per the `ortools` documentation")
    parser.add_argument('--eps_optimal_relative', '--rtol', type=float, default=1.0e-8, help="pdlp solver only. One component of the bound on the duality gap used to check for convergence. If computational resources are limited, consider using `1.0e-4` per the `ortools` documentation")
    parser.add_argument('--primal_feasibility_tolerance', type=float, default=1.0e-8, help="glop solver only.")
    parser.add_argument('--dual_feasibility_tolerance', type=float, default=1.0e-8, help="glop solver only.")
    parser.add_argument('--pdlp_presolve_use_glop', action='store_true', help="pdlp solver only. Use glop's presolve routines but solve with pdlp. Recommended for most cases.")
    parser.add_argument('--loose_solve', action='store_true', help="This will run pdlp (not glop) with weakened termination criteria. Consider using if computational resources are limited.")
    parser.add_argument('--save_model', type=str, default=None, help='Save the optimization model as an MPS file. If specified, the model will be saved to the provided path.')
    parser.add_argument('--scale_gamma', action='store_true', help='If True, the fragmentation penalty (TV) is scaled by the difference in scores at the same position so that the penalty is inversely related to the score difference. Avoids penalizing fragmentation over intervals where the score change is large, and a shift in state is expected. Default is False.')
    parser.add_argument('--scale_gamma_eps','--denom_const', type=float, default=None, help='If `scale_gamma` is True, this value is added to the denominator in the fragmentation penalty to avoid division by zero. Only relevant if `--scale_gamma` is invoked.')
    parser.add_argument('--scale_gamma_beta', '--beta', type=float, default=None, help='If `scale_gamma` is True, this value is used to exponentiate the denominator of the (dynamically-scaled) fragmentation penalty. Only relevant if `--scale_gamma` is also invoked.')


    # Scoring-related arguments
    ## central tendency
    parser.add_argument('--c_1', type=float, default=1.0, help='Score parameter: coefficient for central tendency measure. Assumed positive in the default implementation')
    parser.add_argument('--method_central_tendency', default='quantile', choices=['quantile', 'tmean', 'mean'], help='Central tendency measure. Default is `quantile` with `quantile_value` set to 0.50 (median)')
    parser.add_argument('--quantile_value', type=float, default=0.50, help='Quantile value for central tendency measure--Only applies if `method_central_tendency` is set to `quantile`')
    parser.add_argument('--tprop', type=float, default=0.05, help='Trim proportion for  (`tmean`)--Only applies for trimmed mean central tendency measure.')

    ## dispersion
    parser.add_argument('--c_2', type=float, default=-1.0, help='Score parameter: coefficient for dispersion measure. Assumed negative in the default implementation')
    parser.add_argument('--method_dispersion', default='mad', choices=['mad', 'std', 'iqr', 'tstd'], help='Dispersion measure')

    ## boundary
    parser.add_argument('--c_3', type=float, default=1.0, help='Score parameter: coefficient for boundary measure. Assumed positive in the default implementation')

    ## parametric-sigmoid
    parser.add_argument('--use_parsig', action='store_true', help='Apply `parsig` function to scores. Consider invoking this argument to promote already-integral solutions in the LP relaxed version of the problem, thereby reducing unnecessary dependence on the stochastic ROCCO-RR procedure.')    
    parser.add_argument('--parsig_B', type=float, default=None, help='parsig function `B` parameter')
    parser.add_argument('--parsig_M', type=float, default=None, help='parsig function `M` parameter')
    parser.add_argument('--parsig_R', type=float, default=None, help='parsig function `R` parameter')
    parser.add_argument('--eps_neg', type=float, default=-1.0e-3, help='Negative constant added to scores prior to optimization avoid selection of low-scoring regions only to exhaust the budget. Default is -1.0e-3. Not relevant if scores are already distributed over positive and negative values.')


    # Count track/matrix generation and processing arguments
    parser.add_argument('--step', '-w', type=int, default=50,
                        help='Size of contiguous, fixed-width genomic intervals over which reads are counted (referred to as "loci" in the paper). Larger `--step` values will capture less substructure within peak regions, but this may be desired for some analyses and reduces runtime/memory use, as well. Default is 50 bp. If using BigWig files as input, this value is inferred from data and ignored. For BigWig input, the step must be consistent across all input files.')
    parser.add_argument('--norm_method', default='RPGC',
                        choices=['RPGC', 'CPM', 'RPKM', 'BPM', 'rpgc', 'cpm', 'rpkm', 'bpm'], help='Normalization method. Default is RPGC (Reads Per Genomic Content), for which the `--effective_genome_size` argument is required (default EGS values supplied automatically for supported genomes). If using BigWig files as input, this value has no effect.')
    parser.add_argument('--min_mapping_score', type=int, default=-1, help='Equivalent to samtools view -q. If using BigWig files as input, this value has no effect.')
    parser.add_argument('--flag_include', type=int, default=66, help='Equivalent to samtools view -f.  If using BigWig files as input, this value has no effect.')
    parser.add_argument('--flag_exclude', type=int, default=1284, help='Equivalent to samtools view -F.  If using BigWig files as input, this value has no effect.')
    parser.add_argument('--extend_reads', type=int, default=-1, help='See `deeptools bamCoverage --extendReads`. If using BigWig files as input, this value has no effect.')
    parser.add_argument('--center_reads', action='store_true', help='See `deeptools bamCoverage --centerReads`.  If using BigWig files as input, this value has no effect.')
    parser.add_argument('--ignore_for_norm', nargs='+', default=[], help='Chromosomes to ignore for normalization.  If using BigWig files as input, this value has no effect.')
    parser.add_argument('--scale_factor', type=float, default=1.0,
                        help='bamCoverage scale factor argument. See `deeptools bamCoverage --scaleFactor`.  If using BigWig files as input, this value has no effect.')
    parser.add_argument('--use_existing_bigwigs', action='store_true',
                        help='If True, use existing BigWig files that were generated previously for the same BAM files instead of generating new ones.')
    parser.add_argument('--round_digits', type=int, default=5,
                        help='Number of digits to round values to where applicable')
    parser.add_argument('--use_savgol_filter', action='store_true',
                        help='Use Savitzky-Golay filter (local least squares) on count tracks after normalization.')
    parser.add_argument('--savgol_window_bp', type=int, default=None, help='Window size for Savitzky-Golay filter in base pairs.')
    parser.add_argument('--savgol_order', type=int, default=None, help='Polynomial degree for the least-squares approximation at each step. If None, the degree is set to roughly window_size-3.')
    parser.add_argument('--use_median_filter', action='store_true',
                        help='Apply median filter to count tracks after normalization.')
    parser.add_argument('--median_filter_kernel', type=int, default=None, help='Kernel (window) size for median filter in units of base pairs. If None, the window size is set in `readtracks.py`')
    parser.add_argument('--transform_log_pc', action='store_true',
                        help='If invoked, count matrices will have their elements transformed as `log2(x + c)` where `c` is a constant (see `--log_const`)')
    parser.add_argument('--log_const', type=float, default=None,
                        help='Constant to add before log transformation.')
    parser.add_argument('--transform_local_ratio', action='store_true')
    parser.add_argument('--local_ratio_window_bp', type=int, default=None, help='Window size for local ratio transformation in base pairs.')
    parser.add_argument('--local_ratio_window_steps', type=int, default=None, help='Window size for local ratio transformation in steps.')
    parser.add_argument('--local_ratio_pc', type=float, default=None, help='Add to local_ref when computing local ratio.')


    # post-processing-related arguments
    parser.add_argument('--min_length_bp', type=int, default=None,
                        help='Minimum length of regions to output in the final BED file')
    args = vars(parser.parse_args())
    
    if len(sys.argv)==1 or args['input_files'] is None or len(args['input_files']) == 0:
        parser.print_help(sys.stdout)
        sys.exit(0)
    
    if any([_get_input_type(args['input_files'][i]) == 'bw' for i in range(len(args['input_files']))]):
        bigwig_notice = (
        """\nNote, some read counting/filtering/etc. options cannot be applied to
        BigWig files, even if corresponding CLI arguments have been invoked. Ensure
        that the read counts in the BigWig files are sufficient to represent the
        experiments' data/sample alignments. Alternatively, you can directly supply
        BAM files as input to ROCCO and use the appropriate command-line arguments
        (e.g., `--step`, `--norm_method`, `--flag_include`, etc.) to ensure that the
        read counts are computed according to your preferences.\n"""
        )
        logger.info(bigwig_notice)

    args['solver'] = clean_string(args['solver'])
    args['norm_method'] = clean_string(args['norm_method']).upper()
    if args['genome'] is not None:
        args['genome'] = clean_string(args['genome'])
        if args['genome'] in GENOME_DICT:
            args['effective_genome_size'] = GENOME_DICT[args['genome']]['effective_genome_size']
            args['chrom_sizes_file'] = GENOME_DICT[args['genome']]['sizes_file']
            if args['params'] is None:
                args['params'] = GENOME_DICT[args['genome']]['params']
        else:
            raise ValueError(f'Genome not found: {args["genome"]}.\nAvailable genomes: {list(GENOME_DICT.keys())}.\nYou can also provide genome resources manually (e.g., chromosome sizes file `-s/--chrom_sizes_file)`, effective genome size `--effective_genome_size`.')
    if args['chrom_sizes_file'] is None:
        raise ValueError('A genome with default resources available was not specified with `-g/--genome`, and so a chromosome sizes file must be supplied with `-s/--chrom_sizes_file`')
    if args['effective_genome_size'] is None and args['norm_method'] == 'RPGC':
        raise ValueError('A genome with default resources available was not specified with `-g/--genome`, and so `--norm_method RPGC` normalization, which requires an effective genome size, must be specified along with use of the `--effective_genome_size` argument.')
    bam_files = []
    bigwig_files = []
    for file_ in args['input_files']:
        if not os.path.exists(file_):
            raise FileNotFoundError(f'File not found: {file_}')
        if _get_input_type(file_) == 'bam':
            logger.info(f'Input file: {file_}')
            bw_file_ = generate_bigwig(file_,step=args['step'],
                    effective_genome_size=args['effective_genome_size'],
                    norm_method=args['norm_method'],
                    min_mapping_score=args['min_mapping_score'],
                    flag_include=args['flag_include'],
                    flag_exclude=args['flag_exclude'],
                    extend_reads=args['extend_reads'],
                    center_reads=args['center_reads'],
                    ignore_for_norm=args['ignore_for_norm'],
                    scale_factor=args['scale_factor'],
                    num_processors=args['threads'],
                    overwrite=args['use_existing_bigwigs'])
            bigwig_files.append(bw_file_)
            bam_files.append(file_)
        elif _get_input_type(file_) == 'bw':
            bigwig_files.append(file_)
        else:
            raise ValueError('Input file must be a BAM alignment file or BigWig file')
    logger.info(f'BigWig files: {bigwig_files}')

    chroms_to_process = get_chroms_and_sizes(args['chrom_sizes_file']).keys()
    if args['chroms']:
        chroms_to_process = [
            chrom for chrom in chroms_to_process if chrom in args['chroms']]
    if args['skip_chroms']:
        chroms_to_process = [
            chrom for chrom in chroms_to_process if chrom not in args['skip_chroms']]

    logger.info(f'Chromosomes: {chroms_to_process}')
    tmp_chrom_bed_files = []
    for chrom_ in chroms_to_process:
        chrom_budget = None
        chrom_gamma = None
        #check if chrom in params_df
        if args['params'] is not None:
            logger.info(f'Attempting to use custom parameters file: {args["params"]}')
            params_df = pd.read_csv(args['params'], sep=',')
            if chrom_ in params_df['chrom'].values:
                # assign float values
                chrom_budget = params_df.loc[params_df['chrom'] == chrom_]['budget'].values[0]
                chrom_gamma = params_df.loc[params_df['chrom'] == chrom_]['gamma'].values[0]

        if chrom_budget is None or args['budget'] is not None:
            chrom_budget = args['budget']

        if chrom_gamma is None or args['gamma'] is not None:
            chrom_gamma = args['gamma']

        # Computing chromosome-specific matrix of read counts/densities/enrichments...
        logger.info(f'Generating chromosome matrix: {chrom_}')
        chrom_intervals = None
        chrom_matrix = None
        if args['use_savgol_filter']:
            chrom_intervals, chrom_matrix = generate_chrom_matrix(chrom_, bigwig_files, args['chrom_sizes_file'], args['step'], round_digits=args['round_digits'], filter_type='savgol',savgol_window_bp=args['savgol_window_bp'], savgol_order=args['savgol_order'], transform_log_pc=args['transform_log_pc'], log_const=args['log_const'], transform_local_ratio=args['transform_local_ratio'], local_ratio_window_bp=args['local_ratio_window_bp'], local_ratio_window_steps=args['local_ratio_window_steps'], local_ratio_pc=args['local_ratio_pc'])
        elif args['use_median_filter']:
            chrom_intervals, chrom_matrix = generate_chrom_matrix(chrom_, bigwig_files, args['chrom_sizes_file'], args['step'], round_digits=args['round_digits'], filter_type='median', medfilt_kernel_bp=args['median_filter_kernel'], transform_log_pc=args['transform_log_pc'], log_const=args['log_const'], transform_local_ratio=args['transform_local_ratio'], local_ratio_window_bp=args['local_ratio_window_bp'], local_ratio_window_steps=args['local_ratio_window_steps'], local_ratio_pc=args['local_ratio_pc'])
        else:
            chrom_intervals, chrom_matrix = generate_chrom_matrix(chrom_, bigwig_files, args['chrom_sizes_file'], args['step'], round_digits=args['round_digits'], transform_log_pc=args['transform_log_pc'], log_const=args['log_const'], transform_local_ratio=args['transform_local_ratio'], local_ratio_window_bp=args['local_ratio_window_bp'], local_ratio_window_steps=args['local_ratio_window_steps'], local_ratio_pc=args['local_ratio_pc'])
        if chrom_intervals is None or chrom_matrix is None:
            logger.warning(f'Skipping chromosome {chrom_}... no data found.')
            continue

        logger.info(f'Chromosome {chrom_} Matrix: {chrom_matrix.shape}')

        try:
            step = float(chrom_intervals[1] - chrom_intervals[0])
            chrom_gamma = (chrom_gamma / step)*50.0
        except Exception as e:
            logger.info(f'Could not scale gamma based on step size: {e}')

        
        # Scoring phase
        logger.info(f'Scoring regions: {chrom_}')
        if clean_string(args['method_central_tendency']) == 'quantile':
            ct_scores = score_central_tendency_chrom(chrom_matrix, quantile=args['quantile_value'])
        elif clean_string(args['method_central_tendency']) == 'tmean':
            ct_scores = score_central_tendency_chrom(chrom_matrix, tprop=args['tprop'], method='tmean')
        elif clean_string(args['method_central_tendency']) == 'mean':
            ct_scores = score_central_tendency_chrom(chrom_matrix, method='mean')
        if clean_string(args['method_dispersion']) == 'mad':
            disp_scores = score_dispersion_chrom(chrom_matrix, method='mad')
        elif clean_string(args['method_dispersion']) == 'iqr':
            disp_scores = score_dispersion_chrom(chrom_matrix, method='iqr')
        elif clean_string(args['method_dispersion']) == 'std':
            disp_scores = score_dispersion_chrom(chrom_matrix, method='std')
        elif clean_string(args['method_dispersion']) == 'tstd':
            disp_scores = score_dispersion_chrom(chrom_matrix, tprop=args['tprop'], method='tstd')
        boundary_scores = score_boundary_chrom(ct_scores)
        
        if (args['parsig_B'] is not None or args['parsig_R'] is not None or args['parsig_M'] is not None) and args['use_parsig'] is False:
            args['use_parsig'] = True
        if args['use_parsig'] and args['parsig_B'] is None:
            args['parsig_B'] = (1 - chrom_budget) + 1e-3
        if args['use_parsig'] and args['parsig_R'] is None:
            args['parsig_R'] = 2

        chrom_scores = score_chrom_linear(ct_scores, disp_scores, boundary_scores, gamma=chrom_gamma, c_1=args['c_1'], c_2=args['c_2'], c_3=args['c_3'], eps_neg=args['eps_neg'], parsig_B=args['parsig_B'], parsig_R=args['parsig_R'], parsig_M=args['parsig_M'])
        
        score_output = pformat({
            'Quantile=0.01': round(np.quantile(chrom_scores, q=0.01, method='higher'), 5),
            'Quantile=0.10': round(np.quantile(chrom_scores, q=0.10, method='higher'), 5),
            'Quantile=0.25': round(np.quantile(chrom_scores, q=0.25, method='higher'), 5),
            'Quantile=0.50': round(np.median(chrom_scores), 5),
            'Quantile=0.75': round(np.quantile(chrom_scores, q=0.75, method='higher'), 5),
            'Quantile=0.90': round(np.quantile(chrom_scores, q=0.90, method='higher'), 5),
            'Quantile=0.99': round(np.quantile(chrom_scores, q=0.99, method='higher'), 5)})

        logger.info(f"\nChromosome {chrom_} scores:\n{score_output}\n")
        # Optimization phase
        logger.info(f'Solving LP relaxation using {args["solver"]}')
        logger.info(f'{chrom_}: budget: {chrom_budget}\tgamma: {chrom_gamma}')
        if args['solver'] == 'glop':
            chrom_lp_sol, chrom_lp_score = solve_relaxation_chrom_glop(chrom_scores, budget=chrom_budget, gamma=chrom_gamma, threads=args['threads'], verbose=args['verbose'], glop_dual_feasibility_tolerance=args['dual_feasibility_tolerance'], glop_primal_feasibility_tolerance=args['primal_feasibility_tolerance'], save_model=args['save_model'], scale_gamma=args['scale_gamma'], denom_const=args['scale_gamma_eps'], beta=args['scale_gamma_beta'])
        elif args['solver'] == 'pdlp':
            chrom_lp_sol, chrom_lp_score = solve_relaxation_chrom_pdlp(chrom_scores, budget=chrom_budget, gamma=chrom_gamma, threads=args['threads'], verbose=args['verbose'], pdlp_presolve_options_use_glop=args['pdlp_presolve_use_glop'],
            pdlp_termination_criteria_eps_optimal_absolute=args['eps_optimal_absolute'],
            pdlp_termination_criteria_eps_optimal_relative=args['eps_optimal_relative'],
            pdlp_use_low_precision=args['loose_solve'], save_model=args['save_model'], scale_gamma=args['scale_gamma'], denom_const=args['scale_gamma_eps'], beta=args['scale_gamma_beta'])

        logger.info(f'Refining relaxed solution for integrality')
        chrom_rround_sol, chrom_rround_score = get_rround_sol(chrom_lp_sol, chrom_scores, chrom_budget, chrom_gamma, rand_iter=args['rand_iter'], int_tol=args['int_tol'], eps_mult=args['eps_mult'])
        logger.info(f'\n\nChromosome {chrom_} Optimization Results (Desired: Ratio ~ 1)\
            \nLP-Bound Ideal: {round(chrom_lp_score,4)}\
            \nAttained ROCCO-RR Performance: {round(chrom_rround_score,4)}\
            \nROCCO-RR:LP-Ideal Ratio: {round(chrom_rround_score/chrom_lp_score,8)}\n\n')

        logger.info(f'Chromosome {chrom_}: Writing solution to BED')
        chrom_outfile = chrom_solution_to_bed(chrom_, chrom_intervals, chrom_rround_sol, ID, check_gaps_intervals=True)
        tmp_chrom_bed_files.append(chrom_outfile)

    logger.info('Combining chromosome solutions')
    final_output = combine_chrom_results(tmp_chrom_bed_files, args['output'])
    if os.path.exists(final_output):
        logger.info(f'Final output: {final_output}')

    logger.info('Cleaning up temporary files')
    for tmp_file in tmp_chrom_bed_files:
        try:
            os.remove(tmp_file)
        except Exception as e:
            logger.info(f'Could not remove chromosome-specific temp. file: {tmp_file}\n{e}')
    
if __name__ == '__main__':
    main()
