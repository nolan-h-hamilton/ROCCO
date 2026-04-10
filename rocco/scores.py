r"""
ROCCO: Scores
==================================================================================

(Experimental) Compute *post hoc* peak statistics given samples' BAM files and ROCCO peak output.
"""

import logging
import multiprocessing
import os
from typing import Tuple, Optional
from collections import OrderedDict
from collections.abc import Callable

import numpy as np
import pandas as pd
import pysam
from scipy import ndimage, signal, stats

try:
    from . import _hts_counts
except ImportError:  # pragma: no cover - exercised in source-only environments
    _hts_counts = None

from rocco.readtracks import (
    get_chroms_and_sizes,
    check_type_bam_files,
)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def _random_intervals(
    chrom_sizes_file: str,
    length: int,
    nsamples: int,
    seed: int | None = None,
) -> list[tuple[str, int, int]]:
    chrom_sizes = get_chroms_and_sizes(chrom_sizes_file)
    length_ = int(max(1, length))
    chroms = []
    max_starts = []
    for chrom, chrom_size in chrom_sizes.items():
        max_start = int(chrom_size) - length_ + 1
        if max_start <= 0:
            continue
        chroms.append(str(chrom))
        max_starts.append(int(max_start))
    if len(chroms) == 0:
        raise ValueError(
            f"No chromosome in {chrom_sizes_file} is long enough for intervals of length {length_}."
        )
    weights = np.asarray(max_starts, dtype=np.float64)
    weight_sum = float(np.sum(weights))
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise ValueError("Could not construct a valid random-interval sampler.")
    weights = weights / weight_sum
    rng = np.random.default_rng(seed)
    chrom_indices = rng.choice(
        len(chroms),
        size=int(max(1, nsamples)),
        replace=True,
        p=weights,
    )
    starts = [
        int(rng.integers(0, max_starts[int(chrom_idx)]))
        for chrom_idx in chrom_indices
    ]
    return [
        (chroms[int(chrom_idx)], int(start), int(start + length_))
        for chrom_idx, start in zip(chrom_indices, starts)
    ]


def _require_native_counter():
    if _hts_counts is None:
        raise ImportError(
            "The HTSlib counting code is unavailable. "
            "Reinstall ROCCO so the `rocco._hts_counts` extension is built."
        )
    return _hts_counts


def _read_peak_intervals(
    peak_file: str,
    min_columns: int = 3,
) -> tuple[list[str], list[str], list[int], list[int], list[str]]:
    chroms: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    bed_strings: list[str] = []
    names: list[str] = []

    with open(peak_file, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, start=1):
            line_ = line.strip()
            if line_ == "":
                continue
            fields = line_.split("\t")
            if len(fields) < int(max(3, min_columns)):
                raise ValueError(
                    f"Peak file row {line_num} has fewer than {max(3, min_columns)} columns."
                )
            chrom = str(fields[0])
            start = int(fields[1])
            end = int(fields[2])
            chroms.append(chrom)
            starts.append(start)
            ends.append(end)
            bed_strings.append("\t".join(fields[0:3]))
            names.append("_".join(fields[0:3]))
    return chroms, starts, ends, bed_strings, names


class EmpiricalNull:
    r"""Tiny helper for a finite-sample empirical null.

    We keep the implementation simple on purpose. The right-tail survival gets a
    plus-one correction so the largest observed statistic never receives a zero
    p-value just because we drew a finite number of background regions.
    """

    def __init__(self, values: np.ndarray):
        values_ = np.sort(np.asarray(values, dtype=np.float64))
        if values_.ndim != 1 or values_.size == 0:
            raise ValueError("`values` must be a non-empty one-dimensional array.")
        self.values = values_
        self.size = int(values_.size)

    def survival(self, x):
        x_ = np.asarray(x, dtype=np.float64)
        idx = np.searchsorted(self.values, x_, side="left")
        survival = (self.size - idx + 1.0) / (self.size + 1.0)
        if x_.ndim == 0:
            return float(survival)
        return survival

    def evaluate(self, x):
        x_ = np.asarray(x, dtype=np.float64)
        idx = np.searchsorted(self.values, x_, side="right")
        cdf = idx / float(self.size)
        if x_.ndim == 0:
            return float(cdf)
        return cdf


def _check_read(
    read: pysam.AlignedSegment, min_mapping_quality: int = 10
):
    r"""Called by `pysam.AlignmentFile.count()` to determine if a read is suitable for counting
    :param read: A `pysam` read from a BAM file.
    """
    return (
        not read.is_unmapped
        and read.mapping_quality >= min_mapping_quality
    )


def _null_stat(vals: np.ndarray, percentile: float = 75.0):
    r"""Default statistic for `get_ecdf()` and `score_peaks()`.

    Computes the specified percentile of the input values.
    :param vals: Array of values to compute the null statistic from.
    :type vals: np.ndarray
    :return: The 75th percentile of the input values.
    :rtype: float
    """
    return np.percentile(vals, percentile)


def _peak_signal_stat(
    vals: np.ndarray,
    length: int,
    row_scale: float = 1000.0,
    pc: float = 1.0,
    percentile: float = 75.0,
):
    r"""Use the same transformed statistic for signal and empirical-null scoring."""
    length_ = max(int(length), 1)
    vals_ = np.asarray(vals, dtype=np.float64)
    transformed = np.log2(
        np.maximum(
            vals_ * (float(row_scale) / float(length_)) + float(pc),
            float(pc),
        )
    )
    return float(np.percentile(transformed, percentile))


def _assign_length_bins(
    lengths: np.ndarray,
    max_bins: int = 24,
    min_bin_width_bp: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    lengths_ = np.maximum(np.asarray(lengths, dtype=np.int64), 1)
    if lengths_.ndim != 1 or lengths_.size == 0:
        raise ValueError("`lengths` must be a non-empty one-dimensional array.")

    uniq_lengths = np.unique(lengths_)
    span_bp = int(uniq_lengths[-1] - uniq_lengths[0])
    width_limited_max_bins = 1
    if span_bp >= int(min_bin_width_bp):
        width_limited_max_bins = max(1, span_bp // int(min_bin_width_bp))
    effective_max_bins = max(
        1,
        min(int(max_bins), int(width_limited_max_bins)),
    )
    if uniq_lengths.size <= effective_max_bins:
        return lengths_.astype(np.int64, copy=False), uniq_lengths.astype(
            np.int64,
            copy=False,
        )

    log_edges = np.linspace(
        np.log(float(uniq_lengths[0])),
        np.log(float(uniq_lengths[-1])),
        num=int(effective_max_bins) + 1,
    )
    bin_ids = np.digitize(
        np.log(uniq_lengths.astype(np.float64)),
        log_edges[1:-1],
        right=False,
    )

    length_to_bin: dict[int, int] = {}
    bin_representatives: list[int] = []
    for bin_id in np.unique(bin_ids):
        members = uniq_lengths[bin_ids == bin_id]
        representative = int(np.median(members))
        representative = max(representative, 1)
        bin_representatives.append(representative)
        for length in members:
            length_to_bin[int(length)] = representative

    binned_lengths = np.asarray(
        [length_to_bin[int(length)] for length in lengths_],
        dtype=np.int64,
    )
    return binned_lengths, np.asarray(
        sorted(set(bin_representatives)),
        dtype=np.int64,
    )


def raw_count_matrix(
    bam_files: str,
    peak_file: str,
    output_file: str,
    bed_columns: int = 3,
    overwrite=True,
):
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

    .. todo::
        Parallelize wrt `bam_files` to speed up the counting process.

    """
    bam_files = check_type_bam_files(bam_files)
    native = _require_native_counter()
    chroms, starts, ends, _, peak_names = _read_peak_intervals(
        peak_file,
        min_columns=bed_columns,
    )
    num_peaks = len(peak_names)
    if num_peaks == 0:
        raise ValueError("Peak file does not contain any intervals.")

    samples = []
    for bam in bam_files:
        name = os.path.basename(bam)
        if name.endswith(".bam"):
            name = name[:-4]
        samples.append(name)
    header = "peak_name\t" + "\t".join(samples)
    logger.info(f"\n\nCount matrix header: {header}\n\n")
    logger.info(
        "Counting %s peak regions across %s alignments.",
        int(num_peaks),
        int(len(bam_files)),
    )

    count_matrix = np.zeros((num_peaks, len(bam_files)), dtype=np.int64)
    thread_count = max(multiprocessing.cpu_count() // 2 - 1, 1)
    for sample_idx, bam_file in enumerate(bam_files):
        logger.info(
            "Counting sample %s of %s: %s",
            sample_idx + 1,
            len(bam_files),
            bam_file,
        )
        counts = native.count_alignment_intervals(
            bam_file,
            chroms,
            starts,
            ends,
            one_read_per_bin=1,
            thread_count=thread_count,
            flag_exclude=0,
            min_mapping_quality=10,
            count_mode="coverage",
        )
        count_matrix[:, sample_idx] = np.rint(
            np.asarray(counts, dtype=np.float64)
        ).astype(np.int64, copy=False)

    if output_file is not None and os.path.exists(output_file):
        logger.warning(f"{output_file} already exists...overwriting.")
        os.remove(output_file)

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for peak_idx, peak_name in enumerate(peak_names):
            counts = "\t".join(
                str(int(value)) for value in count_matrix[peak_idx]
            )
            handle.write(f"{peak_name}\t{counts}\n")
    logger.info("Count matrix written to %s", output_file)
    return output_file


def get_read_length(
    bam_file: str,
    num_reads: int = 1000,
    min_mapping_quality: int = 10,
):
    r"""Get the *mapped* read length from a BAM file's first `num_reads`.

    :param bam_file: Path to the BAM file.
    :type bam_file: str
    :param num_reads: Number of reads to sample for the approximation.
    :type num_reads: int
    :return: estimated mapped read length
    :rtype: int
    """

    with pysam.AlignmentFile(
        bam_file,
        "rb",
        threads=max(multiprocessing.cpu_count() // 2 - 1, 1),
    ) as bam:
        read_lengths = []
        for read in bam.fetch():
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
            ):
                continue
            # check mapping quality
            if read.mapping_quality < min_mapping_quality:
                continue
            read_lengths.append(read.infer_query_length())
            if len(read_lengths) >= num_reads:
                break
    return int(np.percentile(read_lengths, 75))


def score_peaks(
    bam_files,
    chrom_sizes_file: str = None,
    peak_file: str = None,
    count_matrix_file: str = None,
    effective_genome_size: float = None,
    skip_for_norm: list = ["chrX", "chrY", "chrM"],
    row_scale=1000,
    ucsc_base=250,
    threads: int = None,
    pc=1,
    ecdf_nsamples=500,
    ecdf_max_length_bins: int = 24,
    output_file="scored_peaks.bed",
    seed: int = None,
    proc: int = None,
    null_stat: Callable[[np.ndarray], float] = _null_stat,
    summit_offsets_file: str | None = None,
):
    r"""Score ROCCO peaks and write narrowPeak-like output.

    The signal column is a robust high quantile of log2-transformed,
    length-scaled counts. The p-values come from an empirical null built on
    that exact same statistic. Lengths are grouped into nearby bins first, so
    we do not waste time fitting a separate ECDF for every exact width.

    :param bam_files: List of paths to the BAM files OR a single filepath to a text file containing a list of BAM files (one filepath per line).
    :type bam_files: list or str
    :param chrom_sizes_file: Path to the chromosome sizes file.
    :type chrom_sizes_file: str
    :param peak_file: Path to (BED) file containing peak regions.
    :type peak_file: str
    :param uscc_base: Base value for UCSC score column. Default is 250 such that no peaks are indiscernible.
    :type uscc_base: int
    :param pc: Pseudocount used inside the log2 transform applied to length-scaled counts.
    :type pc: float
    :param ecdf_max_length_bins: Maximum number of length bins used for empirical-null fitting.
    :type ecdf_max_length_bins: int
    :param summit_offsets_file: Optional TSV mapping `peak_name` to the narrowPeak `peak`
        offset from `chromStart`. Missing peaks default to `-1`.
    :type summit_offsets_file: str or None
    """

    threads_ = 1
    threads_ = (
        threads
        if threads is not None
        else max(multiprocessing.cpu_count() // 2 - 1, 1)
    )

    proc_ = 1
    if proc is None or (isinstance(proc, int) and proc < 1):
        proc_ = min(max(multiprocessing.cpu_count() // 2 - 1, 1), 8)
    else:
        proc_ = proc

    bam_files_ = check_type_bam_files(bam_files)
    matrix_df = None
    matrix_ = None
    try:
        matrix_df = pd.read_csv(
            count_matrix_file, sep="\t", header=0, index_col=0
        )
        matrix_ = matrix_df.values
    except Exception as e:
        if peak_file is None:
            raise e
        logger.info(
            f"Generating count matrix from {len(bam_files_)} BAM files and {peak_file} --> {count_matrix_file}"
        )
        count_matrix_file = raw_count_matrix(
            bam_files_, peak_file, count_matrix_file, bed_columns=3
        )
        matrix_df = pd.read_csv(
            count_matrix_file, sep="\t", header=0, index_col=0
        )
        matrix_ = matrix_df.values
    if matrix_ is None or matrix_df is None:
        raise ValueError(
            f"Failed to generate/read count matrix from {len(bam_files_)} BAM files and {peak_file}."
        )

    # Get an array of peak lengths (in bp) from the bed/peak file
    # ...if we only receive a count_matrix file as input, we can try
    # ...to extract lengths from the index
    lengths = np.zeros(matrix_.shape[0])
    bed_strings = []
    names = []
    try:
        chroms, starts, ends, bed_strings, names = _read_peak_intervals(
            peak_file,
            min_columns=3,
        )
        lengths = np.asarray(
            [end - start for start, end in zip(starts, ends)],
            dtype=np.float64,
        )
    except Exception as e:
        if matrix_df is None:
            raise e

        # if no peak file is provided, try to extract lengths from the count matrix file header
        # this assumes if the index column of the count matrix has peaks named as: `chr_start_end`
        # ...which will be the case if the count matrix was generated by `rocco.readtracks.raw_count_matrix()`
        lengths = np.array(
            [
                int(x.split("_")[2]) - int(x.split("_")[1])
                for x in matrix_df.index
            ]
        )
        bed_strings = [
            "\t".join(x.split("_")[0:3]) for x in matrix_df.index
        ]
        names = [str(x) for x in matrix_df.index]
        logger.info(
            f"Extracted peak lengths from count matrix file: {count_matrix_file}"
        )

    # Normalize each sample's counts so their total coverage is comparable to the 'effective genome size'
    if effective_genome_size is None:
        effective_genome_size = np.sum(
            [
                x[1]
                for x in get_chroms_and_sizes(
                    chrom_sizes_file
                ).items()
                if x[0] not in skip_for_norm
            ]
        )
    mapped_counts = np.zeros(len(bam_files_), dtype=int)
    mapped_rlens = np.zeros(len(bam_files_), dtype=int)
    # get number of mapped counts in chromosomes that are not in `skip_for_norm`
    for i, sample in enumerate(bam_files_):
        aln_sample = pysam.AlignmentFile(
            sample, "rb", threads=threads_
        )
        aln_all_mapped = aln_sample.mapped
        for chrom in skip_for_norm:
            aln_all_mapped -= aln_sample.count(chrom)
        mapped_counts[i] = aln_all_mapped
        mapped_rlens[i] = get_read_length(sample)
    mapped_sizes = mapped_counts * mapped_rlens

    # compute sample scaling constants:
    #   Case 1: (mapped_rlen*mapped_counts) >= effective_genome_size, i.e., 'covers' genome --> scale DOWN to 1x genome
    #   Case 2: (mapped_rlen*mapped_counts < effective_genome_size) < effective_genome_size --> scale UP to 1x genome
    sample_scaling_constants = (effective_genome_size) / mapped_sizes
    for sample_idx in range(len(bam_files_)):
        # scale by (mapped_aln_count*mapped_read_length)
        matrix_[:, sample_idx] = (
            matrix_[:, sample_idx]
            * sample_scaling_constants[sample_idx]
        )

    binned_lengths, ecdf_lengths = _assign_length_bins(
        lengths,
        max_bins=ecdf_max_length_bins,
    )
    logger.info(
        "Using %s ECDF length bins for %s unique peak lengths.",
        int(ecdf_lengths.size),
        int(np.unique(lengths).size),
    )

    # Get empirical cdf for each length bin
    # Generate a random seed if None
    if seed is None:
        seed = np.random.randint(1, 10000)
        logger.info(
            f"Using random seed: {seed} to sample length-binned regions for ECDFs."
        )
    ecdf_dict = multi_ecdf(
        bam_files,
        ecdf_lengths,
        chrom_sizes_file,
        nsamples_per_length=ecdf_nsamples,
        sample_scaling_constants=sample_scaling_constants,
        seed=seed,
        proc=proc_,
        row_scale=row_scale,
        pc=pc,
    )

    # (i) Signal value for a peak is the ~75th percentile~ of 1x-normalized, log2-transformed, length-scaled counts
    # ...observed in each sample over the peak
    sig_vals = np.zeros(matrix_.shape[0])

    # p-value for a peak is the proportion of randomly selected, length-binned regions
    # ...with greater values than the null statistic from random regions
    pvals = np.zeros(matrix_.shape[0])

    for feature_idx in range(len(lengths)):
        if feature_idx % 1000 == 0:
            logger.info(
                f"Processing peak {feature_idx} of {len(lengths)}"
            )
        # (i)
        sig_vals[feature_idx] = _peak_signal_stat(
            matrix_[feature_idx, :],
            lengths[feature_idx],
            row_scale=row_scale,
            pc=pc,
        )
        pvals[feature_idx] = ecdf_dict[binned_lengths[feature_idx]].survival(
            sig_vals[feature_idx]
        )
    scores = sig_vals

    # We will use BH FDR correction to compute q-values
    # ...may later consider an alternative that does not
    # ...assume independence or at least gives users options
    # ...to apply, e.g., yekutieli, bonferroni(FWER), etc.
    qvals = stats.false_discovery_control(pvals, method="bh")

    summit_offsets = {}
    if summit_offsets_file is not None:
        with open(summit_offsets_file, encoding="utf-8") as handle:
            for line_num, line in enumerate(handle, start=1):
                line_ = line.strip()
                if line_ == "":
                    continue
                fields = line_.split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"Summit offset row {line_num} in {summit_offsets_file} has fewer than 2 columns."
                    )
                summit_offsets[str(fields[0])] = int(fields[1])

    # Scale signal values, p-values, q-values according to narrowPeak convention
    bed6_scores = np.minimum(
        np.array(
            ucsc_base
            + sig_vals
            / np.quantile(sig_vals, q=0.99)
            * (1000 - ucsc_base),
            dtype=int,
        ),
        1000,
    )
    pvals_out = np.round(-np.log10(pvals + 1e-10), 4)
    qvals_out = np.round(-np.log10(qvals + 1e-10), 4)
    sig_vals = np.round(sig_vals, 4)

    with open(output_file, "w") as f:
        for i, peak in enumerate(bed_strings):
            summit_offset = int(summit_offsets.get(names[i], -1))
            if summit_offset >= 0:
                summit_offset = int(
                    np.clip(
                        summit_offset,
                        0,
                        max(int(lengths[i]) - 1, 0),
                    )
                )
            f.write(
                f"{peak}\t{names[i]}\t{bed6_scores[i]}\t.\t{sig_vals[i]}\t{pvals_out[i]}\t{qvals_out[i]}\t{summit_offset}\n"
            )
        logger.info(f"Scored output: {output_file}")
    return scores, bed6_scores, pvals


def get_ecdf(
    bam_files,
    length: int,
    chrom_sizes_file: str,
    nsamples=500,
    sample_scaling_constants: list = None,
    seed: int = None,
    null_stat: Callable[[np.ndarray], np.float64] = _null_stat,
    trim_proportion: float = 0.0,
    row_scale: float = 1000.0,
    pc: float = 1.0,
):
    r"""Approximate a null distribution for one representative length bin.

    This samples random genomic intervals of a representative bin length.
    Future versions could condition on more than just width.

    :param bam_files: List of paths to the BAM files OR a single filepath to a text file containing a list of BAM files (one filepath per line).
    :type bam_files: list or str
    :param length: Representative length for the bin that will be sampled.
    :type length: int
    :param chrom_sizes_file: Filepath to chromosome sizes file.
    :type chrom_sizes_file: str
    :param nsamples: Number of random 'background' regions to sample.
    :type nsamples: int
    :param sample_scaling_constants: List of scaling constants (floats) for each sample, typically provided by `score_peaks()`.
    :type sample_scaling_constants: list
    :param seed: Random seed for interval sampling
    :type seed: int
    :return: Empirical null object with `survival()` and `evaluate()` methods.
    :rtype: EmpiricalNull

    :seealso: `score_peaks()`, `multi_ecdf()`
    """

    bam_files_ = check_type_bam_files(bam_files)
    len_avgs = []
    sample_scaling_constants_ = (
        np.ones(len(bam_files_), dtype=np.float64)
        if sample_scaling_constants is None
        else np.asarray(sample_scaling_constants, dtype=np.float64)
    )
    if sample_scaling_constants_.shape[0] != len(bam_files_):
        raise ValueError(
            "`sample_scaling_constants` must match the number of BAM files."
        )
    logger.info(
        f"Computing ECDF for representative length bin: {length} with {nsamples} samples."
    )
    random_intervals = _random_intervals(
        chrom_sizes_file,
        length=int(length),
        nsamples=int(nsamples),
        seed=seed,
    )
    aln_handles = [
        pysam.AlignmentFile(bam_file, "rb")
        for bam_file in bam_files_
    ]
    try:
        for chrom, start, end in random_intervals:
            cperlen = np.zeros(len(bam_files_), dtype=np.float64)
            for j, aln_sample in enumerate(aln_handles):
                cperlen[j] = (
                    aln_sample.count(
                        chrom,
                        start,
                        end,
                        read_callback=_check_read,
                    )
                    * sample_scaling_constants_[j]
                )
            transformed = np.log2(
                np.maximum(
                    np.asarray(cperlen, dtype=np.float64)
                    * (float(row_scale) / float(max(int(length), 1)))
                    + float(pc),
                    float(pc),
                )
            )
            len_avgs.append(null_stat(transformed))
    finally:
        for aln_sample in aln_handles:
            aln_sample.close()
    len_avgs = np.array(len_avgs)
    if trim_proportion > 0:
        len_avgs = stats.trim1(
            len_avgs,
            proportiontocut=trim_proportion,
            tail="right",
        )
    return EmpiricalNull(len_avgs)


def wrap_run_ecdf(*args):
    r"""Wrapper for `get_ecdf()` for compliance with `multiprocessing.Pool`"""
    return get_ecdf(*args)


def multi_ecdf(
    bam_files,
    lengths,
    chrom_sizes_file: str,
    nsamples_per_length,
    sample_scaling_constants = None,
    seed=None,
    proc: int = None,
    null_stat: Callable[[np.ndarray], float] = _null_stat,
    row_scale: float = 1000.0,
    pc: float = 1.0,
):
    r"""Compute ECDFs in parallel for each unique representative length bin.

    See `get_ecdf()` for details.
    """
    np.random.seed(seed)
    bam_files_ = check_type_bam_files(bam_files)
    if proc is None:
        proc = min(max(multiprocessing.cpu_count() // 2 - 1, 1), 8)

    uniq_lengths = np.unique(lengths)
    ecdf_len_dict = OrderedDict.fromkeys(uniq_lengths, None)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(processes=proc) as pool:
        args = [
            (
                bam_files_,
                len_,
                chrom_sizes_file,
                nsamples_per_length,
                sample_scaling_constants,
                seed,
                null_stat,
                0.0,
                row_scale,
                pc,
            )
            for len_ in ecdf_len_dict
        ]
        results = pool.starmap(wrap_run_ecdf, args)

    for idx, result in enumerate(results):
        ecdf_len_dict[uniq_lengths[idx]] = result
    return ecdf_len_dict
