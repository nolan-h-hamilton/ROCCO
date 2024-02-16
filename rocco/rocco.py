r"""
##############################################################################
ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization
##############################################################################

.. image:: ../docs/logo.png
  :width: 400
  :align: center
  :alt: logo

Underlying ROCCO is a :ref:`constrained optimization problem <problem>` with solutions dictating accessible chromatin regions across multiple samples.


GitHub (Homepage)
==================

https://github.com/nolan-h-hamilton/ROCCO/

Paper
========

If using ROCCO in your research, please cite the `original paper <https://doi.org/10.1093/bioinformatics/btad725>`_ in *Bioinformatics*

.. code-block:: text

    Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
    Bioinformatics, Volume 39, Issue 12, December 2023

**DOI**: `10.1093/bioinformatics/btad725 <https://doi.org/10.1093/bioinformatics/btad725>`_


Installation
===============

``pip install rocco``


Example Use
==============

Both an API and command-line interface (CLI) are available to run ROCCO. Toy examples are provided below using both the CLI and API.

API Use
----------------------

Example One
^^^^^^^^^^^^^^^^

Run ROCCO with BAM input files for each sample using the default chromosome-specific budget, gamma, etc. parameters for hg38 assembly (See ``Rocco.HG38_PARAMS``). Compute coverage tracks at intervals of 100 base pairs (Default is 50).

.. doctest::

    >>> import rocco
    >>> bamfiles = ['sample1.bam', 'sample2.bam', 'sample3.bam']
    >>> rocco_obj = rocco.Rocco(
    ...     input_files=bamfiles,
    ...     genome_file='genome.sizes',
    ...     chrom_param_file='hg38',
    ...     step=100)
    >>> rocco_obj.run() # output peak annotation saved as a BED file `rocco_obj.outfile`

Example Two
^^^^^^^^^^^^^^^^

Run ROCCO with BigWig input files for each sample using default chromosome-specific budget, gamma, etc. parameters for the ``hg38`` assembly in ``Rocco.HG38_PARAMS``

Accepting BigWig input directly allows users the option to generate the coverage tracks independently with, e.g., `deepTools bamCoverage <https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html>`_ or
another utility that produces BigWig signal files with additional features for normalization, smoothing, read extension, etc. that are useful in many experimental settings.

.. doctest::

    >>> import rocco
    >>> bw_files = ['sample1.bw', 'sample2.bw', 'sample3.bw']
    >>> rocco_obj = rocco.Rocco(
    ...     input_files=bw_files,
    ...     genome_file='genome.sizes',
    ...     chrom_param_file='hg38')
    >>> rocco_obj.run()

Example Three
^^^^^^^^^^^^^^^^

Scale coverage tracks of samples before calling peaks

.. doctest::

    >>> import rocco
    >>> bamfiles = ['sample1.bam', 'sample2.bam', 'sample3.bam']
    >>> rocco_obj = rocco.Rocco(
    ...     input_files=bamfiles,
    ...     genome_file='genome.sizes',
    ...     chrom_param_file='hg38',
    ...     sample_weights=[0.50, 1.0, 1.0]
    ... )
    >>> rocco_obj.run()

Example Four
^^^^^^^^^^^^^^^^

Use a custom parameter file to set chromosome-specific budgets, gamma, tau, etc.

With the following saved to **``custom_params.csv``**

.. code-block:: text

    chrom,budget,gamma,tau,c_1,c_2,c_3
    chr19,0.05,1.0,0,1.0,1.0,1.0
    chr21,0.03,1.0,0,1.0,1.0,1.0

.. doctest::

    >>> import rocco
    >>> bamfiles = ['sample1.bam', 'sample2.bam', 'sample3.bam']
    >>> rocco_obj = rocco.Rocco(
    ...     input_files=bamfiles,
    ...     genome_file='genome.sizes',
    ...     chrom_param_file='custom_params.csv')
    >>> rocco_obj.run()

For a template chromosome parameter file see ``docs``.


Example Five
^^^^^^^^^^^^^^^^

Use constant, default genome-wide parameters for all chromosomes. Users can modify these default genome-wide parameters via the ``constant_`` arguments at the command line or as below

.. doctest::

    >>> import rocco
    >>> bamfiles = ['sample1.bam', 'sample2.bam', 'sample3.bam']
    >>> rocco_obj = rocco.Rocco(
    ...     input_files=bamfiles,
    ...     genome_file='genome.sizes',
    ...     chrom_param_file=None,
    ...     constant_budget=0.035,
    ...     constant_gamma=1.0,
    ...     constant_tau=0.0)
    >>> rocco_obj.run()


CLI Use
-------------------------------

See :meth:`rocco.main` or run ``rocco --help`` at the command line for a full list of available arguments.

Example One
^^^^^^^^^^^^^^^^

.. code-block:: text

    rocco -i sample1.bam sample2.bam sample3.bam --genome_file genome.sizes --chrom_param_file hg38 --step 100

ROCCO will also accept wildcard/regex filename patterns:

.. code-block:: text

    rocco -i *.bam --genome_file genome.sizes --chrom_param_file hg38 --step 100

Example Two
^^^^^^^^^^^^^^^^

.. code-block:: text

    rocco -i *.bw --genome_file genome.sizes --chrom_param_file hg38


Example Three
^^^^^^^^^^^^^^^^^

List the input files and weights explicitly to ensure the correct order

.. code-block:: text

    rocco -i sample1.bam sample2.bam sample3.bam --genome_file genome.sizes --chrom_param_file hg38 --sample_weights 0.50 1.0 1.0

Example Four
^^^^^^^^^^^^^^^^^

.. code-block:: text

    rocco -i sample1.bam sample2.bam sample3.bam --genome_file genome.sizes --chrom_param_file custom_params.csv

Example Five
^^^^^^^^^^^^^^^^^^

.. code-block:: text

    rocco -i *.bam --genome_file genome.sizes --constant_budget 0.035 --constant_gamma 1.0 --constant_tau 0.0


Testing ROCCO
===============

Run PyTest unit tests

.. code-block::

    cd tests
    pytest -v -rPA -l -k "regular" test_rocco.py


Notes/Miscellaneous
======================

* Default parameters (constant or chromosome-specific) generally provide strong results, but users may consider tweaking the default parameters using a custom ``--chrom_param_file`` or by modifying the default 'constant' genome-wide parameters with the ``--constant_budget``, ``--constant_gamma``, etc. arguments. See `paper and/or supplement <https://doi.org/10.1093/bioinformatics/btad725>`_ for further guidance on optional parameter tuning.

* Run with ``--plot_hist`` to generate a histogram of peak scores that may be useful if tuning the ``--peak_score_filter`` argument.

* If RAM is a special consideration, you can try increasing `--step` from its default of `50` to, e.g., `100` and/or using a lightweight solver for the optimization, e.g., `pip` install `ortools` and run ROCCO with `--solver PDLP`

* Ensure [samtools](https://samtools.github.io) and [bedtools](https://bedtools.readthedocs.io/en/latest/) are installed and in your PATH. These tools are utilized for several auxiliary features.


To-Do Items
======================

.. todolist::


"""
#!/usr/bin/env python
import argparse
import copy
import logging
import math
import multiprocessing
import os
from datetime import datetime
from io import StringIO
from pprint import pformat
import random
import subprocess
import sys
import types
import warnings

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import pybedtools
import pyBigWig as pbw
import scipy

def file_basename(filename):
    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]
    return basename

def nearest_step_idx(val, step):
    return round((val // step)*step)

def get_chroms_and_sizes(genome_file):
    r"""
    Parse chromosomes and their respective sizes in `genome_file`

    :raises FileNotFoundError: If genome file cannot be found
    :return: chromosome names and sizes in `genome_file`
    :rtype: tuple(list,list)
    """
    if not os.path.exists(genome_file) or genome_file is None:
        raise FileNotFoundError(f'Genome file, {genome_file}, not found or is `None`')
    chrom_names = pd.read_csv(genome_file,sep='\t',header=None)[0]
    chrom_sizes = pd.read_csv(genome_file,sep='\t',header=None)[1]
    return dict(zip(chrom_names, chrom_sizes))


class Sample:
    r"""
    Sample

    Used to generate/parse coverage tracks of samples' BAM or BigWig files

    :param input_file: a BAM (.bam), or bigwig (.bw) file for a given sample
    :type input_file: str
    :param genome_file: Path to the genome sizes file containing chromosome sizes
    :type genome_file: str
    :param skip_chroms: List of chromosomes to exclude when parsing/generating coverage tracks
    :type skip_chroms: list, optional
    :param proc_num: Number of processes to use when computing chromosome-specific coverage tracks
    :type proc_num: int, optional
    :param step: Step size for coverage tracks. This is overwritten and inferred from the data if a BigWig file is used as input
    :type step: int, optional
    :param weight: Weight to scale coverage values by. Can be used to apply scaling factor for normalization, etc. Defaults to 1.0
    :type weight: float, optional
    """
    def __init__(self, input_file, genome_file, **kwargs):
        logging.basicConfig(level=logging.INFO, format='LOG: %(asctime)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        logging.info(f"Current working directory: {os.getcwd()}")

        self.input_file = os.path.relpath(input_file)
        self.genome_file = os.path.relpath(genome_file)

        self.proc_num = kwargs.get('proc_num', max(multiprocessing.cpu_count()-1,1))
        self.step = kwargs.get('step', 50)
        self.bamcov_extra_args = kwargs.get('bamcov_extra_args', None)
        self.bamcov_cmd = kwargs.get('bamcov_cmd', None)
        self.weight = kwargs.get('weight', 1.0)

        self.skip_chroms = kwargs.get('skip_chroms', [])
        genome_chrom_sizes_dict = get_chroms_and_sizes(self.genome_file)
        self.chroms = kwargs.get('chroms', [x for x in genome_chrom_sizes_dict.keys() if x not in self.skip_chroms])
        self.chrom_sizes_dict = {}
        for chrom in genome_chrom_sizes_dict:
            if chrom in self.chroms:
                self.chrom_sizes_dict.update({chrom: genome_chrom_sizes_dict[chrom]})

        if self.get_input_type() == 'bam':
            logging.info(f"Calling bam_to_bigwig(): {self.input_file}")
            self.bigwig = self.bam_to_bigwig(bamcov_cmd=self.bamcov_cmd, additional_args=self.bamcov_extra_args)
            
        elif self.get_input_type() == 'bw':
            self.bigwig = self.input_file

    def __str__(self):
        attributes = {}
        for attr, value in vars(self).items():
            attributes[attr] = value
        return pformat(attributes)


    def get_input_type(self):
        r"""
        Determine if self.input_file is a BAM, or BigWig file

        The file type is determined by the file extension: '.bam', '.bw', etc. and is not robust
        to incorrectly labelled files.

        :raises ValueError: If file extension is not supported
        :return: a string (extension) representing the file type
        :rtype: str
        """
        file_type = None
        file_ext = os.path.splitext(self.input_file.lower())[1][1:]

        if file_ext in ['bam']:
            file_type = 'bam'
        elif file_ext in ['bw','bigwig']:
            file_type = 'bw'

        if file_type is None:
            raise ValueError('Input file must be a BAM alignment file or bigwig file')

        return file_type


    def bam_to_bigwig(self, bamcov_cmd=None, additional_args=None):
        r"""
        Wraps deeptools' bamCoverage to generate a BigWig coverage track from the input BAM file
        """
        
        if bamcov_cmd is None:
            bamcov_cmd = ['bamCoverage', '--bam', self.input_file,
                    '--binSize', str(self.step),
                    '-o', f"{self.input_file + '.bw'}",
                    '-p', str(self.proc_num), '--samFlagInclude', '64']
            if additional_args:
                bamcov_cmd.extend(additional_args)
        try:
            subprocess.run(bamcov_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except:
            logging.info(f"{bamcov_cmd} failed. Ensure input is a sorted/indexed BAM file and that deepTools is installed.")
            raise
        return f"{self.input_file + '.bw'}"


    def get_chrom_data(self, chromosome):
        r"""
        Parse a bigwig file and return chromosome-specific coverage loci, coverage values
        """
        input_ = self.bigwig 
        try:
            input_bw = pbw.open(input_)
        except:
            logging.info(f"Could not read {input_} as a BigWig via pyBigWig")
            raise

        loci = []
        vals = []
        idx = 0
        first_nonzero = -1
        for interval in input_bw.intervals(chromosome,0, self.chrom_sizes_dict[chromosome]):
            loci.append(interval[0])
            vals.append(interval[2])
            if interval[2] > 0 and first_nonzero < 0:
                first_nonzero = idx
            idx += 1
        loci = np.array(loci[first_nonzero:])
        vals = self.weight * np.array(vals[first_nonzero:])
        # Step sizes larger than the 'correct' step size are common in bigwig files
        # generated with popular tools from BAM files as a result of compression.
        # Step sizes smaller than the correct step size should be much less common,
        # so, `min()`. Ideally, the BigWig file will explicitly list each interval
        # so that this issue is avoided altogether. 
        step = min([x for x in np.diff(loci[np.nonzero(loci)]) if x > 0])
        if step != self.step:
            warnings.warn(f"Step size inferred from BigWig file ({step}) doesn't match `self.step`. Resetting `self.step=step`.")
            self.step = step
        gap_indices = np.where(np.diff(loci) > step)[0]
        new_loci = []
        new_vals = []

        for gap_index in gap_indices:
            gap_start = loci[gap_index] + step
            gap_end = loci[gap_index + 1]
            fill_val = vals[gap_index]
            gap_loci = np.arange(gap_start, gap_end + step, step)
            new_loci.extend(gap_loci)
            new_vals.extend(fill_val * np.ones_like(gap_loci))
        combined_loci = np.concatenate([loci, new_loci])
        combined_vals = np.concatenate([vals, new_vals])
        sort_index = np.argsort(combined_loci)
        loci_ = combined_loci[sort_index]
        vals_ = combined_vals[sort_index]
        unique_indices = np.unique(loci_, return_index=True)[1]
        loci = loci_[unique_indices]
        vals = vals_[unique_indices]
        step = min([x for x in np.diff(loci[np.nonzero(loci)]) if x > 0])
        return (loci,vals)


class Rocco:

    r"""
    :param input_files: a list of samples' BAM files, or a list of samples' BigWig
    :type input_files: list
    :param genome_file: Path to the genome sizes file containing chromosome sizes
    :type genome_file: str
    :param chrom_param_file: Path to the chromosome parameter file
    :type chrom_param_file: str
    :param skip_chroms: List of chromosomes to skip
    :type skip_chroms: list, optional
    :param proc_num: Number of processes to use when computing chromosome-specific coverage tracks
    :type proc_num: int, optional
    :param step: Step size for coverage tracks. Inferred from the intervals if bigwig input is supplied.
    :type step: int, optional
    :param sample_weights: List/array of weights used to scale the coverage tracks for each sample. Defaults to np.ones()
    :type sample_weights: numpy array, optional
    :param filler_params: Filler parameters used if missing from `chrom_param_file`
    :type filler_params: dict, optional
    :param solver: Solver to use for optimization (default: 'CLARABEL').
    :type solver: str, optional
    :param solver_maxiter: Maximum number of solving iterations
    :type solver_maxiter: int, optional
    :param solver_reltol: Relative optimality gap tolerance when solving the relaxed optimization problem
    :type solver_reltol: float, optional
    :param solver_abstol: Absolute optimality gap tolerance when solving the relaxed optimization problem
    :type solver_abstol: float, optional
    :param solver_feastol: Feasibility tolerance when solving the relaxed optimization problem
    :type solver_feastol: float, optional
    :param rand_iter: Number of RR iterations.
    :type rand_iter: int, optional
    :param verbose_solving: Whether to print solver logging data
    :type verbose_solving: bool, optional

"""

    def __init__(self, input_files, genome_file, chrom_param_file=None, **kwargs):
        logging.basicConfig(level=logging.INFO, format='LOG: %(asctime)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        self.curr_time = datetime.now().strftime('%m%d%Y_%H%M%S')
        self.HG38_PARAMS =\
"""chrom,budget,gamma,tau,c_1,c_2,c_3
chr1,0.03,1.0,0,1.0,1.0,1.0
chr2,0.02,1.0,0,1.0,1.0,1.0
chr3,0.02,1.0,0,1.0,1.0,1.0
chr4,0.02,1.0,0,1.0,1.0,1.0
chr5,0.02,1.0,0,1.0,1.0,1.0
chr6,0.02,1.0,0,1.0,1.0,1.0
chr7,0.025,1.0,0,1.0,1.0,1.0
chr8,0.025,1.0,0,1.0,1.0,1.0
chr9,0.025,1.0,0,1.0,1.0,1.0
chr10,0.02,1.0,0,1.0,1.0,1.0
chr11,0.035,1.0,0,1.0,1.0,1.0
chr12,0.035,1.0,0,1.0,1.0,1.0
chr13,0.02,1.0,0,1.0,1.0,1.0
chr14,0.025,1.0,0,1.0,1.0,1.0
chr15,0.03,1.0,0,1.0,1.0,1.0
chr16,0.03,1.0,0,1.0,1.0,1.0
chr17,0.04,1.0,0,1.0,1.0,1.0
chr18,0.02,1.0,0,1.0,1.0,1.0
chr19,0.045,1.0,0,1.0,1.0,1.0
chr20,0.03,1.0,0,1.0,1.0,1.0
chr21,0.02,1.0,0,1.0,1.0,1.0
chr22,0.03,1.0,0,1.0,1.0,1.0
chrX,0.015,1.0,0,1.0,1.0,1.0
chrY,0.005,1.0,0,1.0,1.0,1.0
"""

        self.MM10_PARAMS=\
"""chrom,budget,gamma,tau,c_1,c_2,c_3
chr1,0.03,1.0,0,1.0,1.0,1.0
chr2,0.03,1.0,0,1.0,1.0,1.0
chr3,0.03,1.0,0,1.0,1.0,1.0
chr4,0.03,1.0,0,1.0,1.0,1.0
chr5,0.035,1.0,0,1.0,1.0,1.0
chr6,0.025,1.0,0,1.0,1.0,1.0
chr7,0.035,1.0,0,1.0,1.0,1.0
chr8,0.03,1.0,0,1.0,1.0,1.0
chr9,0.035,1.0,0,1.0,1.0,1.0
chr10,0.025,1.0,0,1.0,1.0,1.0
chr11,0.04,1.0,0,1.0,1.0,1.0
chr12,0.03,1.0,0,1.0,1.0,1.0
chr13,0.03,1.0,0,1.0,1.0,1.0
chr14,0.025,1.0,0,1.0,1.0,1.0
chr15,0.03,1.0,0,1.0,1.0,1.0
chr16,0.025,1.0,0,1.0,1.0,1.0
chr17,0.03,1.0,0,1.0,1.0,1.0
chr18,0.03,1.0,0,1.0,1.0,1.0
chr19,0.03,1.0,0,1.0,1.0,1.0
chrX,0.02,1.0,0,1.0,1.0,1.0
chrY,0.01,1.0,0,1.0,1.0,1.0
"""
        self.input_files = input_files
        self.genome_file = genome_file
        self.skip_chroms = kwargs.get('skip_chroms', [])
        self.chroms = kwargs.get('chroms', [x for x in get_chroms_and_sizes(self.genome_file).keys() if x not in self.skip_chroms])
        self.peak_score_filter = kwargs.get('peak_score_filter', 0.0)

        self.proc_num = kwargs.get('proc_num', max(multiprocessing.cpu_count()-1,1))
        self.step = kwargs.get('step', 50)
        self.sample_weights = kwargs.get('sample_weights')
        if self.sample_weights is None or self.sample_weights == []:
            self.sample_weights = np.ones(len(self.input_files))
        self.chrom_param_file = chrom_param_file
        self.filler_params = kwargs.get('filler_params',
                                        {'budget':0.035, 'gamma':1.0, 'tau':0.0, 'c_1':1.0, 'c_2':1.0, 'c_3':1.0})

        self.chrom_param_file = chrom_param_file
        self.param_df = None
        expected_columns = ['chrom','budget','gamma','tau','c_1','c_2','c_3']

        if self.chrom_param_file is None:
            csv_rows = ['chrom,budget,gamma,tau,c_1,c_2,c_3']
            for chrom in self.chroms:
                csv_rows.append(','.join([chrom] + [str(self.filler_params[param]) for param in ['budget', 'gamma', 'tau', 'c_1', 'c_2', 'c_3']]))
            gen_csv = '\n'.join(csv_rows)
            csvStringIO = StringIO(gen_csv)

        elif self.chrom_param_file.lower() in ['hg', 'hg38', 'grch38']:
            csvStringIO = StringIO(self.HG38_PARAMS)

        elif self.chrom_param_file.lower() in ['mm','mm10','grcm38']:
            csvStringIO = StringIO(self.MM10_PARAMS)

        elif os.path.exists(self.chrom_param_file):
            param_csv = None
            with open(self.chrom_param_file, "r", encoding="utf-8") as params:
                param_csv = params.read()
            csvStringIO = StringIO(param_csv)

        self.param_df = pd.read_csv(csvStringIO, sep=",")
        selected_columns = list(set(expected_columns) & set(self.param_df.columns))
        for column, filler_value in self.filler_params.items():
            if column not in self.param_df.columns:
                self.param_df[column] = filler_value
            else:
                self.param_df[column].fillna(filler_value, inplace=True)
        self.param_df = self.param_df[~self.param_df['chrom'].isin(self.skip_chroms)]
        self.param_df = self.param_df[selected_columns]


        logging.info(f"Chromosome Parameter Dataframe:\n{self.param_df}")

        self.solver = kwargs.get('solver')
        if self.solver is None:
            self.solver = 'CLARABEL'

        self.solver_maxiter = kwargs.get('solver_maxiter', 10000)
        self.solver_reltol = kwargs.get('solver_reltol', 1.0e-8)
        self.solver_abstol = kwargs.get('solver_abstol', 1.0e-8)
        self.solver_feastol = kwargs.get('solver_feastol', 1.0e-8)
        self.rand_iter = kwargs.get('rand_iter', 100)
        self.verbose_solving = kwargs.get('verbose_solving', True)
        self.keep_chrom_bedfiles = kwargs.get('keep_chrom_bedfiles', False)
        self.eps_l = kwargs.get('eps_l', 1.0e-4)

        samples = []
        for j,input_file in enumerate(self.input_files):
            samples.append(Sample(input_file, self.genome_file, weight=self.sample_weights[j], step=self.step, proc_num=self.proc_num))
        self.samples = samples
        self.step = self.samples[0].step
        self.outfile = kwargs.get('outfile', f"rocco_peaks_{self.curr_time}.bed")
        self.tempfiles = []
        self.pr_bed = kwargs.get('pr_bed', '')
        self.plot_hist = kwargs.get('plot_hist', False)

    def __str__(self):
        attributes = {}
        for attr, value in vars(self).items():
            attributes[attr] = value
        return pformat(attributes)


    def get_Smat(self, chromosome, samples=None):
        r"""
        Generates a :math:`K \times n` coverage signal matrix :math:`\mathbf{S}_{chr}` from samples' data for input to ROCCO

        .. math::

            \mathbf{S}_{chr} = \begin{pmatrix}
            s_{1_1} & s_{2_1}  & \ldots & s_{n_1}\\
            s_{1_2} & s_{2_2}  & \ldots & s_{n_2}\\
            \ldots &\ldots&\ldots&\ldots\\
            s_{1_k} & s_{2_k}  & \ldots & s_{n_K}
            \end{pmatrix}


        :param chromosome: Chromosome over which to compute :math:`\mathbf{S}_{chr}`
        :type chromosome: str
        :param samples: list of Sample objects
        :type samples: list
        :return: a matrix representing the coverage signals of each sample over common loci
        :rtype: numpy.ndarray

        """
        samples_loci = []
        samples_vals = []
        if samples is None:
            samples = self.samples

        for j,samp in enumerate(samples):
            loci,vals = samp.get_chrom_data(chromosome)
            samples_loci.append(list(loci))
            samples_vals.append(list(vals))
        common_loci = sorted([x for x in set.intersection(*map(set,samples_loci))])
        Smat = np.zeros(shape=(len(samples),len(common_loci)))
        for j,samp in enumerate(samples):
            logging.info(f"Constructing Smat: ({j+1}/{len(samples)})")
            samp_chrom_dict = dict(zip(samples_loci[j],samples_vals[j]))
            for i,loc in enumerate(common_loci):
                Smat[j][i] = samp_chrom_dict[loc]
        return common_loci, Smat


    def get_scores(self, chromosome, Smat_chr, c_1=None, c_2=None, c_3=None, tau=None):
        r"""
        Compute locus scores :math:`\mathcal{S}(i),~ \forall i=1\ldots n`

        :param chromosome: Name of chromosome, e.g., 'chr19'
        :type chromosome: str
        :param Smat_chr: Matrix of coverage signals
        :type Smat_chr: np.ndarray

        :return: An :math:`n` -length array of locus scores, each of which quantifies the appeal of selecting the respective locus as open (accessible)
        :rtype: np.ndarray

        """
        if c_1 is None:
            c_1 = self.param_df.loc[self.param_df['chrom'] == chromosome, 'c_1'].values[0]
        if c_2 is None:
            c_2 = self.param_df.loc[self.param_df['chrom'] == chromosome, 'c_2'].values[0]
        if c_3 is None:
            c_3 = self.param_df.loc[self.param_df['chrom'] == chromosome, 'c_3'].values[0]
        if tau is None:
            tau = self.param_df.loc[self.param_df['chrom'] == chromosome, 'tau'].values[0]

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

        med_vec = np.median(Smat_chr,axis=0) # g_1
        mad_vec = scipy.stats.median_abs_deviation(Smat_chr,axis=0) # g_2
        g3_vec_ = g3_vec(med_vec) # g_3
        s_vec = c_1*med_vec - c_2*mad_vec + c_3*g3_vec_
        for i, val in enumerate(s_vec):
            if med_vec[i] <= tau:
                s_vec[i] = 0
        return s_vec - self.eps_l


    def delete_tempfiles(self):
        for fname in self.tempfiles:
            if os.path.exists(fname):
                os.remove(fname)


    def run_rr_proc(self, lp_sol, scores, budget, gamma, rand_iter=None, eps_rand = 1e-4) -> np.ndarray:
        r"""
        Execute the randomization procedure to obtain an integral solution from the
        solution to the relaxed problem,  ``lp_sol``

        :param lp_sol: A solution to the continuous, relaxed optimization problem.
        :type lp_sol: np.nadarray

        :param scores: The score for a given locus quantifies the benefit in selecting it as accessible.

            .. math::

                \mathcal{S}: \mathbb{R}^{K\times 1}_{\geq 0} \rightarrow \mathbb{R}

        :type scores: np.ndarray

        :param budget: Defines the budget constraint to upper bound the proportion of loci selected as open/accessible.

            .. math::
                \sum^{n}_{i=1} \ell_i \leq \lfloor nb \rfloor

            where :math:`\textsf{budget} := b` and :math:`n` is the number of loci (fixed-step contiguous genomic intervals).  In general, a lower budget will result in more conservative results.
        :type budget: float, optional

        :param gamma: Weight for the 'fragmentation penalty' in decision space.

            .. math::
                \gamma \sum^{n-1}_{i=1} |\ell_i - \ell_{i+1}|
            
            In general, higher values of gamma yield fewer but broader peaks. Gamma should be set relative to the step size in the coverage signals, with gamma increasing as the step size decreases, as discussed in the Supplement of the paper.
        :type gamma: type, optional

        :param rand_iter: Number of randomized integral solutions to draw as

            .. math::
                \mathbf{\ell}^{\mathsf{~rand}} \sim \text{Bernoulli}(\mathbf{\ell^{\mathsf{~LP}}})

        :type rand_iter: int, optional

        :param eps_rand: Defines the initial 'floor-epsilon' reference solution with the following procedure

            .. code-block:: python

                eps_cpy = eps_rand
                # initialize as floor_eps solution
                init_sol = np.floor(lp_sol + eps_rand)
                while np.sum(init_sol) > np.floor(n*budget):
                    # loop guarantees the feasibility of solutions
                    eps_cpy = eps_cpy/2
                    init_sol = np.floor(lp_sol + eps_cpy)
                # now start randomization procedure...

        :type eps_rand: float, optional

        """
        def obj(sol, scores, gamma):
            return (-scores@sol
                       + gamma*np.sum(np.abs(np.diff(sol,1))))

        if rand_iter is None:
            rand_iter = self.rand_iter

        n = len(scores)
        eps_cpy = eps_rand
        # initialize as floor_eps solution
        init_sol = np.floor(lp_sol + eps_rand)
        while np.sum(init_sol) > np.floor(n*budget):
            # loop guarantees the feasibility of solutions
            eps_cpy = eps_cpy/2
            init_sol = np.floor(lp_sol + eps_cpy)
        if self.rand_iter <= 0:
            return init_sol

        init_score = obj(sol=init_sol, scores=scores, gamma=gamma)

        nonint_loci = [i for i in range(len(lp_sol)) if lp_sol[i] > 0 and lp_sol[i] < 1]
        rr_sol = init_sol
        best_score = init_score
        for j in range(self.rand_iter):
            ell_rand_n = copy.copy(lp_sol)
            # for efficiency, only round `nonint_loci`
            for idx in nonint_loci:
                if random.random() <= lp_sol[idx]:
                    ell_rand_n[idx] = 1
                else:
                    ell_rand_n[idx] = 0

            score = obj(sol=ell_rand_n, scores=scores, gamma=gamma)

            is_feas = (np.sum(ell_rand_n) <= math.floor(n*budget))
            if is_feas and score < best_score:
                rr_sol = ell_rand_n
                best_score = score
        return rr_sol


    def solve_chrom(self, chromosome,
              common_loci=None,
              Smat_chr=None,
              scores=None,
              budget=None,
              gamma=None,
              tau=None,
              solver=None,
              solver_reltol=None,
              solver_abstol=None,
              solver_feastol=None,
              solver_maxiter=None,
              verbose_solving=None,
              step=None,
              outfile=None,
              pr_bed=None):
        r"""
        Executes ROCCO on a given chromosome.

        In the case users wish to construct the signal matrix :math:`\mathbf{S}_{chr}` (``Smat_chr``), locus scores :math:`\mathcal{S}`, etc.
        with their own custom methods, this function has been written to accept these as parameters, only calling the default methods if left
        as `None`.

        :param chromosome: Name of chromosome on which to execute Rocco
        :type chromosome: str

        :param common_loci: A list/array of loci indices in ``chromosome`` over which all samples' coverage signals are defined. Inferred from data if ``None``.
        :type common_loci: np.ndarray, optional

        :param Smat_chr: A :math:`K \times n` np.ndarray (samples by loci) coverage signal matrix. Inferred from data if ``None``.
        :type Smat_chr: np.ndarray, optional

        :param scores: A list/array of scores for each locus quantifying their appeal for selection as open/accessible.
        :type scores: np.ndarray, optional

        :param budget: Defines the budget constraint to upper bound the proportion of loci selected as open/accessible.

            .. math::

                \sum^{n}_{i=1} \ell_i \leq \lfloor nb \rfloor

            where :math:`\textsf{budget} := b` and :math:`n` is the number of loci (fixed-step contiguous genomic intervals). In general, a lower budget will result in more conservative results.
        :type budget: float, optional

        :param gamma: Weight for the 'fragmentation penalty' in decision space.

            .. math::
                \gamma \sum^{n-1}_{i=1} |\ell_i - \ell_{i+1}|
            
            In general, higher values of gamma yield fewer but broader peaks. Gamma should be set relative to the step size in the coverage signals, with gamma increasing as the step size decreases, as discussed in the Supplement of the paper.
        :type gamma: type, optional

        :param outfile: Name of the chromosome-specific BED output file.
        :type outfile: str, optional

        :return: a tuple (chromosome, selected_loci)
        :rtype: tuple(str, np.ndarray)

        .. _problem:

        :notes:

            The initial, unrelaxed optimization problem ROCCO addresses is the following integer program

            .. math::

                \begin{aligned}
                & \underset{\boldsymbol{\ell}}{\text{ Minimize:}}
                & & f_{\mathsf{IP}}(\boldsymbol{\ell}) = \sum_{i=1}^{n}-\left(\mathcal{S}(i)\cdot\ell_i\right) + \gamma\sum_{i=1}^{n-1} |\ell_i - \ell_{i+1}| \\
                & \text{Subject To:} & &  \text{(i)}~~\sum_{i=1}^{n}\ell_i \leq \lfloor nb\rfloor \\
                & & &  \text{(ii)}~~\ell_i \in \{0,1\}, ~\forall i=1 \ldots n.
                \end{aligned}

            where :math:`\ell_i` denotes the binary decision variable for the :math:`i` -th locus, :math:`\mathcal{S}(i)` denotes the 'score' for the :math:`i` -th locus, etc.
            See paper for more details.
        """
        if Smat_chr or common_loci is None:
            common_loci, Smat_chr = self.get_Smat(chromosome)
        if scores is None:
            scores = self.get_scores(chromosome, Smat_chr)
        if budget is None:
            budget = self.param_df.loc[self.param_df['chrom'] == chromosome, 'budget'].values[0]
        if gamma is None:
            gamma = self.param_df.loc[self.param_df['chrom'] == chromosome, 'gamma'].values[0]
        if tau is None:
            tau = self.param_df.loc[self.param_df['chrom'] == chromosome, 'tau'].values[0]
        if solver is None:
            solver = self.solver
        if solver_reltol is None:
            solver_reltol = float(self.solver_reltol)
        if solver_abstol is None:
            solver_abstol = float(self.solver_abstol)
        if solver_feastol is None:
            solver_feastol = float(self.solver_feastol)
        if solver_maxiter is None:
            solver_maxiter = int(self.solver_maxiter)
        if verbose_solving is None:
            verbose_solving = self.verbose_solving
        if step is None and common_loci is None:
            step = self.step
        elif len(common_loci) > 1:
            step = common_loci[1]-common_loci[0]
        if outfile is None:
            outfile = self.outfile
        if pr_bed is None:
            pr_bed = self.pr_bed

        
        
        n = len(scores)
        logging.info(f"solve_chrom() {chromosome}\nloci: {n}\nbudget: {budget}\ngamma: {gamma}")
        # define problem in CVXPY
        ell = cp.Variable(n)
        ell_aux = cp.Variable(n-1)
        constraints = [ell <= 1,
                          ell >= 0,
                          cp.sum(ell) <= math.floor(budget*n),
                          ell_aux >= cp.diff(ell,1),
                          ell_aux >= -1*cp.diff(ell,1)]
        problem = cp.Problem(cp.Minimize(-scores@ell + gamma*cp.sum(ell_aux)),
                             constraints)


        if solver.lower() == 'clarabel':
            try:
                problem.solve(solver=cp.CLARABEL,
                        tol_gap_rel=solver_reltol,
                        tol_feas=solver_feastol,
                        max_iter=solver_maxiter,
                        verbose=verbose_solving)
            except cp.error.SolverError as ex:
                logging.info(f"Solver {solver} failed")
                raise

        elif solver.lower() == 'pdlp':
            try:
                problem.solve(solver=cp.PDLP, verbose=verbose_solving)
            except cp.error.SolverError as ex:
                logging.info("Ensure a supported version of ortools/PDLP is installed for cvxpy")
                raise
        elif solver.lower() == "ecos":
            problem.solve(solver=cp.ECOS,
                          reltol=solver_reltol,
                          max_iters=solver_maxiter,
                          feastol=solver_feastol,
                          verbose=verbose_solving)
        p_stat = problem.status
        if p_stat is None or p_stat in ['infeasible','unbounded'] or problem.variables()[0].value is None:
            raise cp.error.SolverError(f'\nFailed to obtain optimal solution.\
                \nProblem status: {problem.status}.\n')

        # NOTE:`lp_sol` is the truncated solution to the LP
        lp_sol = problem.variables()[0].value
        rr_sol = self.run_rr_proc(lp_sol, scores, budget, gamma)
        # NOTE: write output as BED file
        tmpfile = f"{outfile}.{chromosome}.tmp"

        selected_loci = []
        with open(tmpfile , 'w') as outfile:
            iter_idx = 0
            for loc,dec_var in zip(common_loci, rr_sol):
                if dec_var > 0:
                    selected_loci.append(loc)
                    outfile.write(f"{chromosome}\t{int(loc)}\t{int(loc+step)}\t{chromosome + '_' + str(int(loc)) + '_' + str(int(loc+step))}\t{np.mean(Smat_chr[:,iter_idx])}\t{'.'}\n")
                iter_idx += 1
        tmpfile_cpy = copy.deepcopy(tmpfile)
        try:
            if pr_bed is not None and len(pr_bed) > 0 and os.path.exists(pr_bed):
                cleaned_bed = pybedtools.BedTool(tmpfile).subtract(pr_bed)
                cleaned_bed.saveas(tmpfile)
        except Exception as ex:
            logging.info(f"The following exception was ignored:\n{ex}\nSkipping removal of problematic regions...")
            tmpfile = tmpfile_cpy
        self.tempfiles.append(tmpfile)
        return (chromosome, np.array(selected_loci))


    def run(self, peak_score_scale = 1000.0):
        r"""
        Execute ROCCO over each given chromosome, merge, score results, and create an output BED file.

        .. todo:: :meth:`Rocco.run`: Update peak score histogram with subplots for each chromosome
        """
        for chrom in self.chroms:
            logging.info(f"Evaluating {chrom}")
            self.solve_chrom(chrom)
        bedtools_list = [pybedtools.BedTool(file_path) for file_path in self.tempfiles]
        first = pybedtools.BedTool(bedtools_list[0])

        pb = first.cat(*bedtools_list[0:],postmerge=False)
        merged_bed = pb.sort().merge(c=5, o='sum')
        peak_scores = []
        with open(self.outfile, 'w') as outfile:
            for feature in merged_bed:
                peak_score = peak_score_scale*(float(feature[3])/(float(feature[2]) - float(feature[1])))
                peak_scores.append(peak_score)
                if peak_score >= self.peak_score_filter:
                    outfile.write(f"{feature[0]}\t{feature[1]}\t{feature[2]}\t{feature[0] + '_' + str(feature[1]) + '_' + str(feature[2])}\t{peak_score}\t.\n")

        if self.plot_hist:
            plt.hist(peak_scores, bins=50, color='blue', edgecolor='black', alpha=0.75, label='Peak score')
            plt.title(f"Peak Score Histogram")
            plt.legend()
            plt.savefig(f"peak_score_hist_{self.curr_time}.pdf")
            plt.close()

        if not self.keep_chrom_bedfiles:
            self.delete_tempfiles()

def main():
    parser = argparse.ArgumentParser(description="ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization", add_help=True)
    parser.add_argument('--input_files', '-i', nargs='+', type=str, help="Samples' corresponding BAM or BigWig files. Accepts wildcard values, e.g., '*.bam', '*.bw'")
    parser.add_argument('--chrom_param_file', type=str, default=None, help="(Optional) Path to CSV param_file OR use `hg38`/`mm10` for human/mouse default parameters. If left unspecified, the 'constant_'/filler parameters are used for all chromosomes.")
    parser.add_argument('--skip_chroms', nargs='+', type=str, default=[], help="Skip these chromosomes")
    parser.add_argument('--genome_file', type=str, help="Genome sizes file. A tab-separated file of the genome's chromosomes and their respective sizes measured in base pairs, e.g., `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes`")
    parser.add_argument('--sample_weights', nargs='+', type=float, default=None)
    parser.add_argument('--pr_bed', type=str, help="BED file of blacklisted/problematic regions to exclude from peak annotation", default=None)
    parser.add_argument('--proc_num', '-p', default=max(multiprocessing.cpu_count()-1,1), type=int,
                        help='Number of processes to run simultaneously when generating coverage signals from BAM files')
    parser.add_argument('--constant_budget', default=0.035, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--constant_gamma', default=1.0, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--constant_tau', default=0.0, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--constant_c_1', default=1.0, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--constant_c_2', default=1.0, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--constant_c_3', default=1.0, type=float, help="'constant' parameters are used to fill in missing or NaN entries in the chromosome-specific parameter files/tables")
    parser.add_argument('--step', default=50, type=int, help='step size in coverage signal tracks. This argument is overwritten and inferred from the intervals if BigWig input is supplied')
    parser.add_argument('--rand_iter', '-N', type=int, default=100, help = 'Number of RR iterations')
    parser.add_argument('--solver', default='CLARABEL', type=str, help='Optimization software used to solve the relaxation')
    parser.add_argument('--solver_maxiter', default=10000, type=int, help='Maximum number of solver iterations')
    parser.add_argument('--solver_reltol', default=1e-8, type=float, help='Maximum allowed relative optimality gap when solving the relaxation')
    parser.add_argument('--solver_feastol', default=1e-8, type=float, help='Maximum allowed feasibility gap when solving the relaxation')
    parser.add_argument('--solver_abstol', default=1e-8, type=float, help='Maximum allowed absolute optimality gap when solving the relaxation')
    parser.add_argument('--peak_score_filter', default = 0.0, type=float, help='Only include peaks in the final annotation with peak scores above `--peak_score_filter`')
    parser.add_argument('--plot_hist', action='store_true', help="If `True` save a plotted histogram of peak scores")
    parser.add_argument('--outfile', '-o',
                        default=f"rocco_peaks_{datetime.now().strftime('%m%d%Y_%H%M%S')}.bed",
                        help='Name of output peak/BED file')
    parser.add_argument('--verbose_solving', action='store_true', default=False)
    args = vars(parser.parse_args())
    
    if args['input_files'] is None:
        print('No input files `--input_files` were supplied. See `rocco --help`')
        sys.exit(1)
    
    if args['genome_file'] is None:
        print('No genome sizes file `--genome_file` was supplied. See `rocco --help`')
        sys.exit(1)
        
    filler_params = {'budget':args['constant_budget'],
                     'gamma':args['constant_gamma'],
                     'tau':args['constant_tau'],
                     'c_1':args['constant_c_1'],
                     'c_2':args['constant_c_2'],
                     'c_3':args['constant_c_3']}

    rocco_obj = Rocco(input_files=args['input_files'],
          genome_file=args['genome_file'],
          chrom_param_file=args['chrom_param_file'],
          skip_chroms=args['skip_chroms'],
          proc_num=args['proc_num'],
          step=args['step'],
          filler_params=filler_params,
          solver=args['solver'],
          sample_weights=args['sample_weights'],
          solver_reltol=args['solver_reltol'],
          solver_maxiter=args['solver_maxiter'],
          solver_feastol=args['solver_feastol'],
          solver_abstol=args['solver_abstol'],
          rand_iter=args['rand_iter'],
          verbose_solving=args['verbose_solving'],
          outfile=args['outfile'],
          pr_bed=args['pr_bed'],
          plot_hist=args['plot_hist'])
    logging.info(rocco_obj)
    rocco_obj.run()

if __name__ == "__main__":
    main()
