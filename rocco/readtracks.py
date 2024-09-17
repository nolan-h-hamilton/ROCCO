r"""
ROCCO: Readtracks
==================================================================================

Basic suite of functions for calculating read
density tracks from BAM files and collating data
for downstream analysis.

"""
import logging
import multiprocessing
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pyBigWig as pbw
import scipy.signal as signal

from typing import Tuple


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                     format='%(filename)s: %(asctime)s - %(levelname)s - %(message)s')


def get_shape(matrix: np.ndarray) -> Tuple:
    r"""Helper function to get the shape of a 1D/2D numpy array"""
    if len(matrix.shape) == 1:
        return 1, len(matrix)
    return matrix.shape


def _next_odd_number(number):
    return number + 1 if number % 2 == 0 else number


def _check_deeptools_installed():
    try:
        subprocess.run(['bamCoverage', '--version'], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        logger.info("Ensure that deepTools is installed and available in your PATH: `pip install deeptools`. `deeptools==3.5.5` is current as of this release. Upgrading matplotlib for conformability with the latest version of deepTools may also resolve issues.")
        raise


def _run_cmd(cmd):
    logger.info(f'Running command: {" ".join(cmd)}\n')
    try:
        proc = subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        raise

def clean_string(string_input):
    if string_input is None:
        return ''
    return string_input.lower().replace(' ', '')


def get_chroms_and_sizes(chrom_sizes_file):
    r"""Parse a chromosome sizes file and return a dictionary of chromosome names and sizes.
    :param chrom_sizes_file: Path to the chromosome sizes file.
    :type chrom_sizes_file: str
    
    :return: A dictionary of chromosome names and sizes.
    :rtype: dict
    
    :raises FileNotFoundError: If the chromosome sizes file is not found.
    """
    if not os.path.exists(chrom_sizes_file):
        raise FileNotFoundError(
            f'Sizes file, {chrom_sizes_file}, not found or is `None`')
    try:
        chrom_names = pd.read_csv(chrom_sizes_file, sep='\t', header=None)[0]
        chrom_sizes = pd.read_csv(chrom_sizes_file, sep='\t', header=None)[1]
    except Exception as e:
        logging.info(f'Error reading chromosome sizes file: {chrom_sizes_file}.\
            \nExpected format is tab-separated with two columns, e.g.,\
            \nchr1\t248956422\nchr2\t242193529\n...')
        raise
    return dict(zip(chrom_names, chrom_sizes))


def _get_count_track_name(bam_file: str, effective_genome_size: int = -1, step: int = 50, norm_method: str = 'RPGC', min_mapping_score: int = -1, flag_include: int = 66, flag_exclude: int = 1284, extend_reads: int = -1, center_reads: bool = False, scale_factor: float = 1.0) -> str:
    r"""Generate a BigWig file name based on the parameters used in `generate_bigwig`.

    .. note:
        The file name will not reflect the `ignore_for_norm` parameter used in `generate_bigwig`.

    """
    basename = os.path.basename(bam_file).replace('.bam', '')
    norm_suffix = f'step_{step}_norm_{norm_method}'
    if effective_genome_size > 0:
        norm_suffix += f'_EGS_{int(effective_genome_size)}'
    if min_mapping_score > 0:
        norm_suffix += f'_minMappingScore_{min_mapping_score}'
    if flag_include > 0:
        norm_suffix += f'_flagInclude_{flag_include}'
    if flag_exclude > 0:
        norm_suffix += f'_flagExclude_{flag_exclude}'
    if extend_reads >= 0:
        if extend_reads > 0:
            norm_suffix += f'_extendReads_{extend_reads}'
        else:
            norm_suffix += f'_extendReads'
    if center_reads:
        norm_suffix += f'_centerReads'
    if scale_factor != 1.0:
        norm_suffix += f'_scaleFactor_{scale_factor}'
    return f"{basename}_{norm_suffix}.bw"


def generate_bigwig(bam_file: str, step: int = 50,
                    effective_genome_size: float = -1,
                    norm_method: str = 'RPGC',
                    min_mapping_score: int = -1,
                    flag_include: int = 66,
                    flag_exclude: int = 1284,
                    extend_reads: int = -1,
                    center_reads: bool = False,
                    ignore_for_norm: list = ['chrX', 'chrY', 'chrM'],
                    scale_factor: float = 1.0,
                    ouput_file: str = None, num_processors: int = -1,
                    overwrite: bool = False) -> str:
    r"""Generate a BigWig file from a BAM file.
    
    :param bam_file: Path to the BAM file.  
    :type bam_file: str
    :param step: Step size for the intervals. Used synonymously with 'bin size' throughout this documentation.
    :type step: int
    :param effective_genome_size: Effective genome size for normalization. Required for RPGC normalization.
    :type effective_genome_size: float
    :param norm_method: Normalization method. Must be one of: RPGC, RPKM, CPM, BPM. See `deeptools` documentation for more information.
    :type norm_method: str
    :param min_mapping_score: Minimum mapping score for reads to be included in the count track.
    :type min_mapping_score: int
    :param flag_include: SAM flag to include in the count track. Equivalent to `-f` in `samtools view`.
    :type flag_include: int
    :param flag_exclude: SAM flag to exclude from the count track. Equivalent to `-F` in `samtools view`.
    :type flag_exclude: int
    :param extend_reads: Extend reads by a fixed number of base pairs. If < 0, ignore. If 0, extend reads to 'estimated fragment length'. If > 0, extend reads by that number of base pairs.
    :type extend_reads: int
    :param center_reads: Center reads at the midpoint of the fragment length.
    :type center_reads: bool
    :param ignore_for_norm: List of chromosome names to ignore for normalization.
    :type ignore_for_norm: list
    :param scale_factor: Factor to scale the values by.
    :type scale_factor: float
    :param ouput_file: Path to the output BigWig file. If `None`, the file name will be generated based on the parameters used.
    :type ouput_file: str
    :param num_processors: Number of processes to use for the calculation. If < 0, use all but one of the available processors.
    :type num_processors: int
    :param overwrite: If `True`, overwrite the output file if it already exists. If `False`, skip the calculation if the output file already exists.
    :type overwrite: bool
    :return: Path to the output BigWig file.
    :rtype: str

    :raises FileNotFoundError: If the BAM file is not found.
    :raises ValueError: If the normalization method is not one of: RPGC, RPKM, CPM, BPM.
    :raises ValueError: If effective genome size is not specified and the normalization method is RPGC.

    """

    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if norm_method not in ['RPGC', 'RPKM', 'CPM', 'BPM']:
        if norm_method.upper() not in ['RPGC', 'RPKM', 'CPM', 'BPM']:
            raise ValueError(
                f'Normalization method must be one of: RPGC, RPKM, CPM')
        else:
            norm_method = norm_method.upper()

    if norm_method == 'RPGC' and effective_genome_size < 0:
        raise ValueError(
            'Effective genome size must be provided for RPGC normalization.')

    if num_processors < 1:
        num_processors = max(multiprocessing.cpu_count() - 1, 1)

    bw_outfile = ouput_file if ouput_file else _get_count_track_name(bam_file, effective_genome_size,
                                                                     step=step, norm_method=norm_method,
                                                                     min_mapping_score=min_mapping_score,
                                                                     flag_include=flag_include, flag_exclude=flag_exclude, extend_reads=extend_reads,
                                                                     center_reads=center_reads, scale_factor=scale_factor)

    if os.path.exists(bw_outfile) and not overwrite:
        last_modified_time = time.strftime("%a, %d %b %Y %H:%M:%S",
                                           time.localtime(os.path.getmtime(bw_outfile)))
        logger.info(f"Count track {bw_outfile} for {bam_file} already exists.\
            \nLast modified: {last_modified_time}.")
        return bw_outfile

    _check_deeptools_installed()

    cmd = ['bamCoverage', '-b', bam_file, '-o', bw_outfile,
           '--binSize', str(step), '--normalizeUsing', norm_method]
    if effective_genome_size > 0:
        cmd += ['--effectiveGenomeSize', str(int(effective_genome_size))]
    if min_mapping_score > 0:
        cmd += ['--minMappingQuality', str(min_mapping_score)]
    if flag_include > 0:
        cmd += ['--samFlagInclude', str(flag_include)]
    if flag_exclude > 0:
        cmd += ['--samFlagExclude', str(flag_exclude)]
    if extend_reads == 0:
        cmd += ['--extendReads']
    if extend_reads > 0:
        cmd += ['--extendReads', str(extend_reads)]
    if center_reads:
        cmd += ['--centerReads']
    if ignore_for_norm:
        cmd += ['--ignoreForNormalization', ' '.join(ignore_for_norm)]
    if scale_factor != 1.0:
        cmd += ['--scaleFactor', str(scale_factor)]
    cmd += ['--numberOfProcessors', str(num_processors)]
    _run_cmd(cmd)

    if not os.path.exists(bw_outfile):
        logger.warning(f"Cannot find bw file: {bw_outfile}.")
    return bw_outfile


def decompress_features_vals(intervals: np.ndarray, vals: np.ndarray, step: int):
    r"""Reformat data from BedGraph-like intervals to have fixed-width, *contiguous* intervals. To avoid redundancy, methods generating BedGraph-like ouput will often create a single entry for repeated
    values over multiple adjacent intervals. For the use-cases in scope of this module, the savings achieved
    by this compressed representation are of little practical value and can be annoying to work with
    downstream. So this function decompresses the intervals and values.
    
    :param intervals: list-like of integer values corresponding to genomic positions
    :type intervals: np.ndarray
    :param vals: list-like of floats representing the measured signal value at each interval
    :type vals: np.ndarray
    :param step: integer specifying the step size in base pairs. 'step' is used
        synonymously with 'bin size'.
    :type step: int
    :return: A tuple of np.ndarray arrays: the first is the gap-filled intervals,
        the second is the signal values over the new intervals
    :rtype: tuple

    :seealso: :func:`get_chrom_reads()`
    :seealso: :func:`generate_chrom_matrix()`

    """

    gap_indices = np.where(np.diff(intervals) > step)[0]
    new_intervals = []
    new_vals = []
    # loop through gaps, fill in with the last non-zero value
    for gap_index in gap_indices:
        gap_start = intervals[gap_index] + step
        gap_end = intervals[gap_index + 1]
        fill_val = vals[gap_index]
        gap_intervals = np.arange(gap_start, gap_end + step, step)
        new_intervals.extend(gap_intervals)
        new_vals.extend(fill_val * np.ones_like(gap_intervals))

    # combine the original intervals and values with the filled gaps
    combined_intervals = np.concatenate([intervals, new_intervals])
    # ensure intervals are of type int
    combined_intervals = combined_intervals.astype(int)
    combined_vals = np.concatenate([vals, new_vals])

    sort_index = np.argsort(combined_intervals)
    intervals_ = combined_intervals[sort_index]
    vals_ = combined_vals[sort_index]
    unique_indices = np.unique(intervals_, return_index=True)[1]
    return intervals_[unique_indices].astype(int), vals_[unique_indices]


def get_chrom_reads(bigwig_file: str, chromosome: str, chrom_sizes_file: str,
                    step: int, const_scale: float = 1.0, round_digits: int = 5,
                    scale_by_step: bool = False):
    r"""Extract, decompress, and scale signal values (e.g., read counts) from a BigWig file for a specific chromosome.

    :param bigwig_file: Path to the BigWig file.
    :type bigwig_file: str
    :param chromosome: Chromosome name to extract data for.
    :type chromosome: str
    :param chrom_sizes_file: Path to the chromosome sizes file.
    :type chrom_sizes_file: str
    :param step: Step size for the intervals.
    :type step: int
    :param const_scale: Factor to scale the signal values by after bigwig generation. Note, this is not the same
        parameter as the scale factor used in `generate_bigwig`.
    :type const_scale: float
    :param round_digits: Number of decimal places to round the wig values to (not intervals).
    :type round_digits: int
    :param scale_by_step: If True, scale the values by the step size.
    :type scale_by_step: bool

    :return: A tuple containing two np arrays. One for the genomic intervals, one for the values.
    :rtype: tuple

    :raises FileNotFoundError: If the BigWig file or chromosome sizes file is not found.
    :raises ValueError: If no non-zero values are found in the BigWig file.
    
    :seealso: :func:`decompress_features_vals()`
    :seealso: :func:`generate_chrom_matrix()`

    """

    if not os.path.exists(bigwig_file):
        raise FileNotFoundError(f"BigWig file not found: {bigwig_file}")
    if not os.path.exists(chrom_sizes_file):
        raise FileNotFoundError(
            f"Chromosome sizes file not found: {chrom_sizes_file}")

    chrom_sizes_dict = get_chroms_and_sizes(chrom_sizes_file)

    if chromosome not in chrom_sizes_dict:
        raise ValueError(
            f"Chromosome {chromosome} not found in chromosome sizes file: {chrom_sizes_file}")
    try:
        with pbw.open(bigwig_file) as input_bw:
            if chromosome not in input_bw.chroms():
                logger.warning(
                    f"Chromosome {chromosome} not found in BigWig file: {bigwig_file}. Returning (None,None).")
                return None, None
            intervals = []
            vals = []
            idx = 0
            first_nonzero = -1
            for interval in input_bw.intervals(
                    chromosome, 0, chrom_sizes_dict[chromosome]):
                intervals.append(interval[0])
                vals.append(interval[2])
                if interval[2] > 0 and first_nonzero < 0:
                    first_nonzero = idx
                idx += 1
    except Exception as e:
        logger.info(f"Error reading BigWig file with pyBigWig: {bigwig_file}\n{e}")
        raise

    if first_nonzero < 0:
        logger.warning(f"No non-zero values found in BigWig file: {bigwig_file} for chromosome: {chromosome}. Returning (None,None).")
        return None, None

    intervals = np.array(intervals[first_nonzero:])
    vals = np.array(vals[first_nonzero:])

    observed_step = int(np.min(np.abs(np.diff(intervals))))
    if observed_step != step:
        logger.warning(f"Observed step size in BigWigs: {observed_step}. Expected step size: {step}.")
    step = observed_step

    # ensure contiguous intervals
    intervals, vals = decompress_features_vals(intervals, vals, step)

    if scale_by_step:
        vals = vals / step
        logger.info(f"Scaling values by step size: {step}")

    if const_scale >= 0:
        if const_scale == 0:
            logger.warning("You are scaling the values by 0.")
        vals = vals * const_scale
    vals = np.round(vals, round_digits)
    return intervals, vals


def generate_bigwigs(bam_files: list, step: int = 50,
                     effective_genome_size: float = -1,
                     norm_method: str = 'RPGC',
                     min_mapping_score: int = -1,
                     flag_include: int = 66,
                     flag_exclude: int = 1284,
                     extend_reads: int = -1,
                     center_reads: bool = False,
                     ignore_for_norm: list = ['chrX', 'chrY', 'chrM'],
                     scale_factor: float = 1.0, num_processors: int = -1,
                     overwrite: bool = True) -> list:
    r"""Generate BigWig files from a list of BAM files. This function is just a wrapper for :func:`generate_bigwig()` 
    
    :param bam_files: List of paths to the BAM files.
    :type bam_files: list
    :param step: Step size for the intervals.
    :type step: int
    :param effective_genome_size: Effective genome size for normalization. Required for RPGC normalization.
    :type effective_genome_size: float
    :param norm_method: Normalization method. Must be one of: RPGC, RPKM, CPM, BPM. See `deeptools` documentation for more information.
    :type norm_method: str
    :param min_mapping_score: Minimum mapping score for reads to be included in the count track.
    :type min_mapping_score: int
    :param flag_include: SAM flag to include in the count track. Equivalent to `-f` in `samtools view`.
    :type flag_include: int
    :param flag_exclude: SAM flag to exclude from the count track. Equivalent to `-F` in `samtools view`.
    :type flag_exclude: int
    :param extend_reads: Extend reads by a fixed number of base pairs. If < 0, ignore. If 0, extend reads to 'estimated fragment length'. If > 0, extend reads by that number of base pairs.
    :type extend_reads: int
    :param center_reads: Center reads at the midpoint of the fragment length.
    :type center_reads: bool
    :param ignore_for_norm: List of chromosome names to ignore for normalization.
    :type ignore_for_norm: list
    :param scale_factor: Factor to scale the values by.
    :type scale_factor: float
    :param num_processors: Number of processes to use for the calculation. If < 0, use all but one of the available processors.
    :type num_processors: int
    :param overwrite: If `True`, overwrite the output file if it already exists. If `False`, skip the calculation if the output file already exists.
    :type overwrite: bool
    :return: List of paths to the output BigWig files.
    :rtype: list
    
    :raises FileNotFoundError: If a BAM file is not found.
    :raises ValueError: If the normalization method is not one of: RPGC, RPKM, CPM, BPM.
    :raises ValueError: If effective genome size is not specified and the normalization method is RPGC.
    
    :seealso: :func:`generate_bigwig()`
        
    """
    bigwig_files = []
    for bam_file in bam_files:
        bigwig_files.append(generate_bigwig(bam_file, step=step,
                                            effective_genome_size=effective_genome_size,
                                            norm_method=norm_method,
                                            min_mapping_score=min_mapping_score,
                                            flag_include=flag_include,
                                            flag_exclude=flag_exclude,
                                            extend_reads=extend_reads,
                                            center_reads=center_reads,
                                            ignore_for_norm=ignore_for_norm,
                                            scale_factor=scale_factor,
                                            ouput_file=None, num_processors=num_processors,
                                            overwrite=overwrite))
    return bigwig_files


def generate_chrom_matrix(chromosome: str, bigwig_files: list, chrom_sizes_file: str, step: int,
                          const_scale: float=1.0, round_digits: int = 5, scale_by_step: bool = False, 
                          filter_type: str=None, savgol_window_bp: int=None, savgol_window_steps: int=5, 
                          savgol_order: int=None, medfilt_kernel_bp: int=None, medfilt_kernel_steps: int=5, 
                          transform_log_pc: bool=False, log_const: float=None, transform_local_ratio: bool=False, 
                          local_ratio_window_bp: int=None, local_ratio_window_steps: int=None, 
                          local_ratio_pc: float=None):
    r"""
    Create a matrix of read counts for a given chromosome from a list of BigWig files with potentially 
    varying start and end positions. Provides options for transforming and filtering the tracks comprising the matrix.

    :param chromosome: The chromosome to generate the matrix for.
    :type chromosome: str
    :param bigwig_files: List of paths to the BigWig files.
    :type bigwig_files: list
    :param chrom_sizes_file: Path to the chromosome sizes file.
    :type chrom_sizes_file: str
    :param step: Step size for the intervals.
    :type step: int
    :param const_scale: Factor to scale the signal values by after bigwig generation.
    :type const_scale: float
    :param round_digits: Number of decimal places to round the wig values to (not intervals).
    :type round_digits: int
    :param scale_by_step: If True, scale the values by the step size.
    :type scale_by_step: bool
    :param filter_type: Type of filter to apply. Options are 'savitzky-golay' or 'median'.
    :type filter_type: str
    :param savgol_window_bp: Window size for the Savitzky-Golay filter in base pairs.
    :type savgol_window_bp: int
    :param savgol_window_steps: Window size for the Savitzky-Golay filter in steps.
    :type savgol_window_steps: int
    :param savgol_order: Polynomial order for the Savitzky-Golay filter.
    :type savgol_order: int
    :param medfilt_kernel_bp: Kernel size for the median filter in base pairs.
    :type medfilt_kernel_bp: int
    :param medfilt_kernel_steps: Kernel size for the median filter in steps.
    :type medfilt_kernel_steps: int
    :param transform_log_pc: If True, apply a log2(x + c) transformation to data.
    :type transform_log_pc: bool
    :param log_const: Constant to add before taking the log.
    :type log_const: float
    :param transform_local_ratio: If True, apply local ratio transformation to the values. See :func:`apply_transformation()`.
    :type transform_local_ratio: bool
    :param local_ratio_window_bp: Window size for the local ratio transformation in base pairs.
    :type local_ratio_window_bp: int
    :param local_ratio_window_steps: Window size for the local ratio transformation in steps.
    :type local_ratio_window_steps: int
    :param local_ratio_pc: Constant to add to the local reference before taking the ratio.
    :type local_ratio_pc: float
    
    :return: A tuple containing two numpy arrays. The first array contains genomic intervals, and the second contains the corresponding values.
    :rtype: tuple(numpy.ndarray, numpy.ndarray)

    :raises FileNotFoundError: If a BigWig file or chromosome sizes file is not found.

    :seealso: :func:`get_chrom_reads()`
    :seealso: :func:`apply_transformation`
    :seealso: :func:`apply_filter`

    """

    # get sync'd intervals for all bigwigs
    interval_matrix = []
    vals_matrix = []
    for bigwig_file in bigwig_files:
        intervals_, vals_ = get_chrom_reads(bigwig_file,
                                            chromosome,
                                            chrom_sizes_file,
                                            step,
                                            const_scale,
                                            round_digits,
                                            scale_by_step)
        if intervals_ is None or vals_ is None:
            logger.warning(f"No data found for {bigwig_file} in chromosome {chromosome}. Excluding this track for {chromosome}.")
            continue
        interval_matrix.append(intervals_)
        vals_matrix.append(vals_)
    if len(interval_matrix) == 0:
        logger.warning(f"No data found in the files {str(bigwig_files)} for chromosome {chromosome}. Returning (None,None).")
        return None, None
    common_intervals = np.sort(
        np.unique(np.concatenate(interval_matrix, axis=0)))
    # initialize with zeroes
    count_matrix = np.zeros((len(interval_matrix), len(common_intervals)))
    for i, (intervals_, vals_) in enumerate(zip(interval_matrix, vals_matrix)):
        # all rows must have same length
        idx = np.searchsorted(common_intervals, intervals_)
        count_matrix[i, idx] = vals_

    if get_shape(count_matrix)[0] == 1:
        count_matrix = count_matrix.reshape(1, -1)

    # apply transformations first
    _, count_matrix = apply_transformation(common_intervals, count_matrix, transform_log_pc=transform_log_pc, log_const=log_const,
                                            transform_local_ratio=transform_local_ratio, local_ratio_window_bp=local_ratio_window_bp,
                                            local_ratio_window_steps=local_ratio_window_steps, local_ratio_pc=local_ratio_pc)
    
    # apply filters
    _, count_matrix = apply_filter(common_intervals, count_matrix, filter_type=filter_type, savgol_window_bp=savgol_window_bp, savgol_window_steps=savgol_window_steps, savgol_order=savgol_order, medfilt_kernel_bp=medfilt_kernel_bp, medfilt_kernel_steps=medfilt_kernel_steps)

    return np.array(common_intervals).astype(int), count_matrix


def apply_filter(intervals: np.ndarray, count_matrix: np.ndarray, filter_type: str = None, savgol_window_bp: int = None, savgol_window_steps: int = 5, savgol_order: int = None, medfilt_kernel_bp: int = None, medfilt_kernel_steps: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    r"""Apply filtering to the count matrix.
    
    :param intervals: The genomic intervals.
    :type intervals: np.ndarray
    :param count_matrix: The count matrix.
    :type count_matrix: np.ndarray
    :param filter_type: The type of filter to apply. Options are 'savitzky-golay' or 'median'.
    :type filter_type: str
    :param savgol_window_bp: Window size for the Savitzky-Golay filter in base pairs.
    :type savgol_window_bp: int
    :param savgol_window_steps: Window size for the Savitzky-Golay filter in steps.
    :type savgol_window_steps: int
    :param savgol_order: Polynomial order for the Savitzky-Golay filter.
    :type savgol_order: int
    :param medfilt_kernel_bp: Kernel size for the median filter in base pairs.
    :type medfilt_kernel_bp: int
    :param medfilt_kernel_steps: Kernel size for the median filter in steps.
    :type medfilt_kernel_steps: int
    
    """

    if filter_type is None:
        if savgol_window_bp is not None or savgol_order is not None:
            filter_type = 'savitzky-golay'
        elif medfilt_kernel_bp is not None:
            filter_type = 'median'
        else:
            return intervals, count_matrix
    
    step_sizes = np.unique(np.abs(np.diff(intervals)))
    if len(step_sizes) > 1:
        logger.warning(f"Warning: non-uniform step sizes detected: {step_sizes}. This may lead to unexpected results.")
    step_size = step_sizes[0]

    if clean_string(filter_type) in ['savitzky-golay', 'savitzky', 'golay', 'savgol']:
        if savgol_window_bp is None or savgol_window_bp < 0: 
            savgol_window_bp = step_size * savgol_window_steps
        filter_window = _next_odd_number(savgol_window_bp//step_size)
        
        if savgol_order is None or savgol_order < 0:
            savgol_order = filter_window//2
        filter_order = max(0, savgol_order)
        logger.info(f"Applying Savitzky-Golay filter with window size: {filter_window*step_size}bp and polynomial degree: {filter_order}")
        for i in range(count_matrix.shape[0]):
            count_matrix[i] = np.array([max(0,x) for x in signal.savgol_filter(count_matrix[i], filter_window, filter_order)])

    elif clean_string(filter_type) in ['median', 'med', 'medfilt', 'median_filter', 'med_filt']:
        if medfilt_kernel_bp is None or medfilt_kernel_bp < 0:
            medfilt_kernel_bp = step_size * medfilt_kernel_steps
        filter_kernel = _next_odd_number(medfilt_kernel_bp//step_size)
        logger.info(f"Applying median filter with kernel size: {filter_kernel*step_size}bp")
        for i in range(count_matrix.shape[0]):
            count_matrix[i] = signal.medfilt(count_matrix[i], filter_kernel)

    return intervals, count_matrix


def apply_transformation(intervals: np.ndarray, count_matrix: np.ndarray,
                        transform_log_pc: bool = False, log_const: float = None,
                        transform_local_ratio: bool = False, local_ratio_window_bp: int = None,
                        local_ratio_window_steps: int = None, local_ratio_pc: float = None) -> Tuple[np.ndarray, np.ndarray]:
    r"""Transform the count matrix after computing read densities.

    :param intervals: The genomic intervals.
    :type intervals: np.ndarray
    :param count_matrix: The count matrix.
    :type count_matrix: np.ndarray
    :param transform_log_pc: If True, apply a `log2(x + c)` transformation to data.
    :type transform_log_pc: bool
    :param log_const: Constant to add before taking the log.
    :type log_const: float
    :param transform_local_ratio: If True, apply local ratio transformation to the values: `x/(median_filter(x,row,kernel=local_ratio_window_bp) + ε)`.
    :type transform_local_ratio: bool
    :param local_ratio_window_bp: Window size for the local ratio transformation in base pairs.
    :type local_ratio_window_bp: int
    :param local_ratio_window_steps: Window size for the local ratio transformation in steps.
    :type local_ratio_window_steps: int
    :param local_ratio_pc: Constant to add to the local reference before taking the ratio.
    :type local_ratio_pc: float

    :return: A tuple containing two numpy arrays. The first array contains genomic intervals, and the second is a matrix (2d array) containing the corresponding values.
    :rtype: tuple(numpy.ndarray, numpy.ndarray)
    
    """
    if transform_log_pc or log_const is not None:
        if log_const is None:
            log_const = 0.25
        logger.info(f"Transforming data as log2(x + {log_const})")
        count_matrix = np.log2(count_matrix + log_const)
    
    if transform_local_ratio or local_ratio_window_bp is not None:
        if local_ratio_window_bp is None or local_ratio_window_bp < 0:
            local_ratio_window_bp = 5000
        step_sizes = np.unique(np.abs(np.diff(intervals)))
        if len(step_sizes) > 1:
            logger.warning(f"Warning: non-uniform step sizes detected: {step_sizes}. This may lead to unexpected results.")
        step_size = step_sizes[0]
        local_ratio_window_steps = _next_odd_number(local_ratio_window_bp // step_size)
        logger.info(f"Transforming data as ratio x/(median_filter(x,row,kernel={local_ratio_window_bp}bp) + ε)")
        for i in range(count_matrix.shape[0]):

            local_ref = signal.medfilt(count_matrix[i], local_ratio_window_steps)
            if local_ratio_pc is None:
                local_ratio_pc = 1.0
            local_ref = local_ref + local_ratio_pc
            count_matrix[i] = (count_matrix[i]) / (local_ref)

    return intervals, count_matrix