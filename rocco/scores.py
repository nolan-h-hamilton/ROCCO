r"""
ROCCO: Scores
==================================================================================

Compute *post hoc* peak statistics given samples' BAM files and ROCCO peak output.

"""

import logging
import multiprocessing
import os
import subprocess
from typing import Tuple
from collections import OrderedDict

import numpy as np
import pandas as pd
import pysam
from scipy import ndimage, signal, stats

from rocco.readtracks import get_chroms_and_sizes, check_type_bam_files 


logging.basicConfig(level=logging.INFO,
                     format='%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING,
                    format='%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def _check_read(read: pysam.AlignedSegment, min_mapping_quality: int=30):
    r"""Called by `pysam.AlignmentFile.count()` to determine if a read is suitable for counting
    :param read: A `pysam` read from a BAM file.
    """
    return not read.is_unmapped and read.mapping_quality >= min_mapping_quality


def _gmean(vals, trim_left: float=None, trim_right: float=None):
    r"""Helper function to compute geometric mean for a list of (optionally) truncated values.
    """
    truncated_vals = None
    if trim_left is not None and trim_right is None:
        truncated_vals = vals[vals > np.quantile(vals, trim_left)]
    elif trim_left is None and trim_right is not None:
        truncated_vals = vals[vals < np.quantile(vals, 1-trim_right)]
    elif trim_left is not None and trim_right is not None:
        truncated_vals = vals[(vals > np.quantile(vals, trim_left)) & (vals < np.quantile(vals, 1-trim_right))]
    if truncated_vals is not None:
        return np.power(np.prod(truncated_vals), 1/len(truncated_vals))
    return np.power(np.prod(vals), 1/len(vals))


def raw_count_matrix(bam_files: str, peak_file: str, output_file: str, bed_columns: int=3, overwrite=True):
    r"""Generate a 'raw' count matrix from BAM files and a ROCCO output peak file.

    :param bam_files: List of paths to the BAM files OR a single filepath to a text file containing a list of BAM files (one filepath per line).
    :type bam_files: list or str
    :param bam_list_file: list of filepaths to BAM files or a single filepath to a textfile containing a list of BAM files (one per line, order will be preserved in the output matrix).
    :type bam_list_file: str
    :param peak_file: Path to the BED file containing peaks.
    :type bed_file: str
    :param output_file: name of file to write the count matrix in.
    :type output_file: str
    :return: Name/path of the output file if successful, otherwise None.
    :rtype: str

    .. note::
        Output Header. 
          The header of the output peak-count files (TSV) will preserve the order that the BAM files were supplied in.
    .. note::
        BED3 format.
          Default expects three-column BED peak files. If the peak file contains additional columns, set `bed_columns` accordingly.

    """
    bam_files = check_type_bam_files(bam_files)
    samples = []
    
    for bam in bam_files:
        name = os.path.basename(bam)
        if name.endswith(".bam"):
            name = name[:-4]
        samples.append(name)
    header = "peak_name\t" + "\t".join(samples)
    logger.info(f"\n\nCount matrix header: {header}\n\n")

    cmd = f"bedtools multicov -bams {' '.join(bam_files)} -bed {peak_file}"
    linecount = sum(1 for _ in open(peak_file))
    linecount_tenpct = int(linecount * 0.10) - 1
    logger.info(f"Running `bedtools multicov` on {len(bam_files)} alignments and {linecount} peak regions.")
    if output_file is not None and os.path.exists(output_file):
        if overwrite:
            logger.warning(f"{output_file} already exists...overwriting.")
            os.remove(output_file)
        else:
            logger.warning(f"{output_file} already exists...returning None.")
            return None
    elif output_file is None:
        output_file = peak_file + ".raw_counts.tsv"
        logger.info(f"Output file not specified...writing to {output_file}.")
    with open(output_file, 'a', encoding='utf-8') as f:
        f.write(header + "\n")
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', bufsize=1) as proc:
            i = 0
            for line in proc.stdout:
                line = line.strip()
                fields = line.split("\t")
                peak_name = '_'.join(fields[0:3])
                counts = fields[bed_columns:]
                f.write(f"{peak_name}\t" + "\t".join(counts) + "\n")
                if i % linecount_tenpct == 0:
                    print(f"Processed {i} of {linecount} peaks.")
                i += 1
    logger.info(f"Count matrix written to {output_file}")
    return output_file


def get_read_length(bam_file: str, num_reads: int = 1000, min_mapping_quality: int = 10):
    r"""Get the *mapped* read length from a BAM file's first `num_reads`.
    
    :param bam_file: Path to the BAM file. 
    :type bam_file: str
    :param num_reads: Number of reads to sample for the approximation.
    :type num_reads: int
    :return: estimated mapped read length
    :rtype: int
    """

    with pysam.AlignmentFile(bam_file, "rb", threads=max(multiprocessing.cpu_count()//2 - 1,1)) as bam:
        read_lengths = []
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            # check mapping quality 
            if read.mapping_quality < min_mapping_quality:
                continue
            read_lengths.append(read.infer_query_length())
            if len(read_lengths) >= num_reads:
                break
    return int(np.percentile(read_lengths, 75))


def score_peaks(bam_files, chrom_sizes_file: str=None,
                peak_file: str=None,
                count_matrix_file: str=None,
                effective_genome_size: float=None,
                skip_for_norm: list = ['chrX', 'chrY', 'chrM'],
                row_scale=1000, ucsc_base=500,
                threads: int=None, pc=1, ecdf_nsamples=500,
                output_file = 'scored_peaks.bed', seed: int=None):
    r"""Compute peak scores based on geometric mean of 1x-normalized, length-scaled counts. p-values based on ECDFs of read counts in each peak length.
    
    :param bam_files: List of paths to the BAM files OR a single filepath to a text file containing a list of BAM files (one filepath per line).
    :type bam_files: list or str
    :param chrom_sizes_file: Path to the chromosome sizes file.
    :type chrom_sizes_file: str
    :param peak_file: Path to (BED) file containing peak regions. 
    :type peak_file: str
    
    .. note::
        Selection of length-specific 'null' regions.
          Should be revisited in future verisons. Though the current approach appears to yield approximately uniformly-distributed p-values
          under the null, there's no explicit check to verify the sampled regions are 'background'. At the cost of some execution
          speed, we might consider checking :math:`A_{xx}(\tau=0) > a\left(\max_{\tau > 0} A_{xx}(0)\right)`,
          :math:`\max - \min < \textsf{threshold}`, distance to nearest peak, etc. The `nullranges` R package provides a rich set of features for sampling null regions--consider integrating in future efforts to migrate to R.

    """
    # check type of bam_files
    bam_files_ = check_type_bam_files(bam_files)

    matrix_df = None
    matrix_ = None
    try:
        matrix_df = pd.read_csv(count_matrix_file, sep="\t", header=0, index_col=0)
        matrix_ = matrix_df.values
    except Exception as e:
        if peak_file is None:
            raise e
        logger.info(f"Generating count matrix from {len(bam_files_)} BAM files and {peak_file} --> {count_matrix_file}")
        count_matrix_file = raw_count_matrix(bam_files_, peak_file, count_matrix_file, bed_columns=3)
        matrix_df = pd.read_csv(count_matrix_file, sep="\t", header=0, index_col=0)
        matrix_ = matrix_df.values
    if matrix_ is None or matrix_df is None:
        raise ValueError(f"Failed to generate/read count matrix from {len(bam_files_)} BAM files and {peak_file}.")


    # Get an array of peak lengths (in bp) from the bed/peak file
    # ...if we only receive a count_matrix file as input, we can try
    # ...to extract lengths from the index
    lengths = np.zeros(matrix_.shape[0])
    bed_strings = []
    names = []
    try:
        with open(peak_file) as f:
            for i, line in enumerate(f):
                fields = line.strip().split("\t")
                lengths[i] = int(fields[2]) - int(fields[1])
                bed_string = '\t'.join(fields[0:3])
                name = '_'.join(fields[0:3])
                bed_strings.append(bed_string)
                names.append(name)
    except Exception as e:
        if matrix_df is None:
            raise e

        # if no peak file is provided, try to extract lengths from the count matrix file header
        # this assumes if the index column of the count matrix has peaks named as: `chr_start_end`
        # ...which will be the case if the count matrix was generated by `rocco.readtracks.raw_count_matrix()`
        lengths = np.array([int(x.split('_')[2]) - int(x.split('_')[1]) for x in matrix_df.index])
        bed_strings = ['\t'.join(x.split('_')[0:3]) for x in matrix_df.index]
        names = [str(x) for x in matrix_df.index]
        logger.info(f"Extracted peak lengths from count matrix file: {count_matrix_file}")

    # Normalize each sample's counts so their total coverage is comparable to the 'effective genome size'
    if effective_genome_size is None:
        effective_genome_size = np.sum([x[1] for x in get_chroms_and_sizes(chrom_sizes_file).items() if x[0] not in skip_for_norm])
    threads = threads if threads is not None else max(multiprocessing.cpu_count()//2 - 1,1)
    mapped_counts = np.zeros(len(bam_files_), dtype=int)
    mapped_rlens = np.zeros(len(bam_files_), dtype=int)
    # get number of mapped counts in chromosomes that are not in `skip_for_norm`
    for i,sample in enumerate(bam_files_):
        aln_sample = pysam.AlignmentFile(sample, "rb", threads=threads)
        aln_all_mapped = aln_sample.mapped
        for chrom in skip_for_norm:
            aln_all_mapped -= aln_sample.count(chrom)
        mapped_counts[i] = aln_all_mapped
        mapped_rlens[i] = get_read_length(sample)
    mapped_sizes = mapped_counts * mapped_rlens

    # compute sample scaling constants:
    #   Case 1: (mapped_rlen*mapped_counts) >= effective_genome_size, i.e., 'covers' genome --> scale DOWN to 1x genome
    #   Case 2: (mapped_rlen*mapped_counts < effective_genome_size) < effective_genome_size --> scale UP to 1x genome
    sample_scaling_constants = (effective_genome_size)/mapped_sizes
    for sample_idx in range(len(bam_files_)):
        # scale by (mapped_aln_count*mapped_read_length)
        matrix_[:, sample_idx] = matrix_[:, sample_idx] * sample_scaling_constants[sample_idx]


    # Get empirical cdf for each unique peak length
    # Generate a random seed if None
    if seed is None:
        seed = np.random.randint(1, 10000)
        logger.info(f"Using random seed: {seed} to sample length-matched regions for ECDFs.")
    ecdf_dict = multi_ecdf(bam_files, lengths, chrom_sizes_file, nsamples_per_length=ecdf_nsamples, sample_scaling_constants=sample_scaling_constants, seed=seed)

    # (i) Signal value for a peak is the geometric mean of the 1x-normalized, length-scaled counts
    # ...observed in each sample over the peak
    sig_vals = np.zeros(matrix_.shape[0])

    # p-value for a peak is the proportion of randomly selected, length-matched regions
    # ...with greater enrichment than the mean of counts observed in the peak over samples.
    pvals = np.zeros(matrix_.shape[0])

    for feature_idx in range(len(lengths)):
        if feature_idx % 1000 == 0:
            logger.info(f"Processing peak {feature_idx} of {len(lengths)}")
        # (i)
        sig_vals[feature_idx] = _gmean((matrix_[feature_idx, :]*(row_scale/lengths[feature_idx])) + pc)
        # (ii)
        pvals[feature_idx] = 1 - ecdf_dict[lengths[feature_idx]].evaluate(np.mean(matrix_[feature_idx,:]))
    scores = sig_vals
    
    # We will use BH FDR correction to compute q-values
    # ...may later consider an alternative that does not
    # ...assume independence or at least gives users options
    # ...to apply, e.g., yekutieli, bonferroni(FWER), etc.
    qvals = stats.false_discovery_control(pvals, method='bh')

    # Scale signal values, p-values, q-values according to narrowPeak convention
    bed6_scores = np.minimum(np.array(ucsc_base + sig_vals/np.quantile(sig_vals,q=0.975)*(1000 - ucsc_base), dtype=int),1000)
    pvals_out = np.round(-np.log10(pvals + 1e-10),4)
    qvals_out = np.round(-np.log10(qvals + 1e-10),4)
    sig_vals = np.round(sig_vals,4)


    with open(output_file, 'w') as f:
        for i, peak in enumerate(bed_strings):
            # for now, call summit as midpoint of peak
            mp_ = int((int(peak.split('\t')[1]) + int(peak.split('\t')[2]))/2)
            f.write(f"{peak}\t{names[i]}\t{bed6_scores[i]}\t.\t{sig_vals[i]}\t{pvals_out[i]}\t{qvals_out[i]}\t{mp_}\n")
        logger.info(f"Scored output: {output_file}")
    return scores, bed6_scores, pvals


def get_ecdf(bam_files, length: int, chrom_sizes_file: str,
             nsamples=500, sample_scaling_constants: list = None, seed: int=None):
    r"""Compute an empirical cdf of counts from randomly selected genomic regions of `length` bases
    
    In the current implementation, this function is called for each unique peak length by `score_peaks()`
    so that the empirical cdf used to approximate p-values is based on regions of the same length as the peaks.
    
    :param bam_files: List of paths to the BAM files OR a single filepath to a text file containing a list of BAM files (one filepath per line).
    :type bam_files: list or str
    :param length: Defines length of regions that will be randomly sampled to estimate the background distribution
    :type length: int
    :param chrom_sizes_file: Filepath to chromosome sizes file.
    :type chrom_sizes_file: str
    :param nsamples: Number of random 'background' regions to sample.
    :type nsamples: int
    :param sample_scaling_constants: List of scaling constants (floats) for each sample (provided by `score_peaks()` in current implementation).
    :type sample_scaling_constants: list
    :param seed: Random seed supplied to `bedtools random`
    :type seed: int
    :return: ECDF object with an `evaluate()` method.
    :rtype: scipy.stats._distn_infrastructure.rv_frozen

    :seealso: `score_peaks()`, `multi_ecdf()`
    """

    bam_files_ = check_type_bam_files(bam_files)
    len_avgs = []
    cmd = f"bedtools random  -l {length} -n {nsamples} -g {chrom_sizes_file}"
    if seed is not None:
        cmd = cmd + f" -seed {seed}"
    logger.info(f"Computing ECDF for length: {length} with {nsamples} samples.")
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', bufsize=1) as proc:
        sample_vals = []
        for line in proc.stdout:
            line = line.strip()
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            cperlen = np.zeros(len(bam_files_))
            for j,bam_file in enumerate(bam_files_):
                with pysam.AlignmentFile(bam_file, "rb") as aln_sample:
                    cperlen[j] = aln_sample.count(chrom, start, end, read_callback=_check_read)*sample_scaling_constants[bam_files_.index(bam_file)]
            sample_vals.extend(cperlen)
    len_avgs = np.array(sample_vals)
    return stats.ecdf(len_avgs).cdf


def wrap_run_ecdf(*args):
    r"""Wrapper for `get_ecdf()` for compliance with `multiprocessing.Pool`"""
    return get_ecdf(*args)


def multi_ecdf(bam_files, lengths, chrom_sizes_file: str,
                nsamples_per_length, sample_scaling_constants=None, seed=None,
                proc: int=None):
    r"""Compute ECDFs in parallel for each unique peak length and a correspondding dictionary of `scipy` ECDF objects.
    
    See `get_ecdf()` for details.
    """
    bam_files_ = check_type_bam_files(bam_files)
    if proc is None:
        proc = min(max(multiprocessing.cpu_count() - 2,1), 8)

    uniq_lengths = np.unique(lengths)
    ecdf_len_dict = OrderedDict.fromkeys(uniq_lengths, None)
    ctx = multiprocessing.get_context('fork')
    with ctx.Pool(processes=proc) as pool:
        args = [(bam_files_, len_, chrom_sizes_file, nsamples_per_length, sample_scaling_constants, seed) for len_ in ecdf_len_dict.keys()]
        results = pool.starmap(wrap_run_ecdf, args)
    
    for idx, result in enumerate(results):
        ecdf_len_dict[uniq_lengths[idx]] = result
    return ecdf_len_dict