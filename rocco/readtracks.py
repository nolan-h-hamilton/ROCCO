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
from typing import Dict, Tuple

import numpy as np
import pandas as pd

try:
    from . import _hts_counts
except ImportError:  # pragma: no cover - exercised in build environments
    _hts_counts = None


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

_BAM_COUNT_METADATA_CACHE: dict[tuple, Dict[str, float | int | bool]] = {}


def get_shape(matrix: np.ndarray) -> Tuple:
    r"""Helper function to get the shape of a 1D/2D numpy array"""
    if len(matrix.shape) == 1:
        return 1, len(matrix)
    return matrix.shape


def _resolve_num_processors(num_processors: int) -> int:
    if num_processors is None or int(num_processors) < 1:
        return max(multiprocessing.cpu_count() - 1, 1)
    return int(num_processors)


def _resolve_parallel_bam_counting(
    bam_count: int,
    num_processors: int,
) -> Tuple[int, int]:
    threads = _resolve_num_processors(num_processors)
    if int(bam_count) <= 1 or int(threads) <= 1:
        return 1, int(threads)
    if "fork" not in multiprocessing.get_all_start_methods():
        return 1, int(threads)
    # Keep this small. Each process is holding one chromosome slice.
    process_count = min(int(bam_count), int(threads), 4)
    # Split the core budget across BAMs so we do not launch N processes that each try to use every thread.
    threads_per_bam = max(1, int(threads) // int(process_count))
    return int(process_count), int(threads_per_bam)


def _require_native_counter():
    if _hts_counts is None:
        raise ImportError(
            "Can't find ROCCO-native htslib-based counting extension...try building from source `python -m build; python -m pip install .`"
        )
    return _hts_counts


def _get_track_type(track_file: str) -> str:
    ext = os.path.splitext(track_file)[1].lower().lstrip(".")
    if ext == "bam":
        return "bam"
    raise ValueError(f"Unsupported input file type for `{track_file}`. Expected BAM.")


def _estimate_fragment_length(
    bam_file: str,
    flag_exclude: int = 0,
    max_samples: int = 4096,
    num_processors: int = 1,
) -> int | None:
    native = _require_native_counter()
    fragment_length = int(
        native.get_alignment_fragment_length(
            bam_file,
            thread_count=max(1, int(num_processors)),
            flag_exclude=max(0, int(flag_exclude)),
            max_iterations=max(1, int(max_samples)),
            fallback=0,
        )
    )
    if fragment_length <= 0:
        return None
    return fragment_length


def _compute_native_scale_factor(
    norm_method: str,
    effective_genome_size: float,
    step: int,
    mapped_reads: int,
    norm_read_length: int,
    scale_factor: float = 1.0,
) -> float:
    norm_method_ = clean_string(norm_method).upper()
    mapped_reads_ = max(int(mapped_reads), 1)
    tile_len_kb = float(step) / 1000.0
    scale = float(scale_factor)
    if norm_method_ == "RPGC":
        if effective_genome_size is None or float(effective_genome_size) <= 0:
            raise ValueError(
                "Effective genome size must be positive for RPGC normalization."
            )
        current_coverage = (
            float(mapped_reads_) * float(max(int(norm_read_length), 1))
        ) / float(effective_genome_size)
        return float(scale * (1.0 / max(current_coverage, 1.0e-12)))
    if norm_method_ == "RPKM":
        million_reads_mapped = float(mapped_reads_) / 1.0e6
        return float(scale * (1.0 / max(million_reads_mapped * tile_len_kb, 1.0e-12)))
    if norm_method_ in {"CPM", "BPM"}:
        million_reads_mapped = float(mapped_reads_) / 1.0e6
        return float(scale * (1.0 / max(million_reads_mapped, 1.0e-12)))
    raise ValueError(
        f"Normalization method must be one of `RPGC`, `RPKM`, `CPM`, or `BPM`, not `{norm_method}`."
    )


def _get_bam_count_metadata(
    bam_file: str,
    step: int,
    norm_method: str,
    effective_genome_size: float,
    ignore_for_norm: list | None,
    flag_exclude: int = 0,
    extend_reads: int = -1,
    num_processors: int = 1,
    scale_factor: float = 1.0,
) -> Dict[str, float | int | bool]:
    ignore_for_norm_ = tuple(ignore_for_norm or [])
    threads = _resolve_num_processors(num_processors)
    cache_key = (
        bam_file,
        int(step),
        clean_string(norm_method).upper(),
        float(effective_genome_size if effective_genome_size is not None else -1.0),
        ignore_for_norm_,
        int(flag_exclude),
        int(extend_reads),
        int(threads),
        float(scale_factor),
    )
    if cache_key in _BAM_COUNT_METADATA_CACHE:
        return _BAM_COUNT_METADATA_CACHE[cache_key]

    native = _require_native_counter()
    paired_end = bool(
        native.is_alignment_paired_end(
            bam_file,
            max_reads=1024,
            thread_count=threads,
        )
    )
    read_length = int(
        native.get_alignment_read_length(
            bam_file,
            min_reads=32,
            thread_count=threads,
            max_iterations=4096,
            flag_exclude=max(0, int(flag_exclude)),
        )
    )
    mapped_reads, _ = native.get_alignment_mapped_read_count(
        bam_file,
        exclude_chromosomes=list(ignore_for_norm_),
        thread_count=threads,
        count_mode="coverage",
        one_read_per_bin=0,
    )

    norm_read_length = int(read_length)
    resolved_extend_bp = int(extend_reads)
    paired_end_mode = False
    if int(extend_reads) == 0:
        fragment_length = _estimate_fragment_length(
            bam_file,
            flag_exclude=max(0, int(flag_exclude)),
            num_processors=threads,
        )
        if paired_end:
            if fragment_length is not None and fragment_length > 0:
                norm_read_length = int(fragment_length)
                paired_end_mode = True
                resolved_extend_bp = 0
            else:
                logger.warning(
                    "Could not estimate fragment length for %s; falling back to read length %s.",
                    bam_file,
                    read_length,
                )
        else:
            if fragment_length is not None and fragment_length > int(read_length):
                norm_read_length = int(fragment_length)
                resolved_extend_bp = int(fragment_length)
                logger.info(
                    "Using inferred single-end fragment length %s for %s.",
                    fragment_length,
                    bam_file,
                )
            else:
                logger.warning(
                    "`extend_reads=0` requests fragment-length inference, but %s did not yield a larger single-end fragment length; using read length %s.",
                    bam_file,
                    read_length,
                )
                resolved_extend_bp = -1
    elif int(extend_reads) > 0:
        norm_read_length = int(extend_reads)
        resolved_extend_bp = int(extend_reads)

    norm_scale = _compute_native_scale_factor(
        norm_method=norm_method,
        effective_genome_size=effective_genome_size,
        step=step,
        mapped_reads=int(mapped_reads),
        norm_read_length=int(norm_read_length),
        scale_factor=float(scale_factor),
    )
    metadata = {
        "paired_end": paired_end,
        "paired_end_mode": paired_end_mode,
        "read_length": int(read_length),
        "norm_read_length": int(norm_read_length),
        "resolved_extend_bp": int(resolved_extend_bp),
        "mapped_reads": int(mapped_reads),
        "norm_scale": float(norm_scale),
        "threads": int(threads),
    }
    _BAM_COUNT_METADATA_CACHE[cache_key] = metadata
    return metadata


def clean_string(string_input):
    if string_input is None:
        return ""
    return string_input.lower().replace(" ", "")


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
            f"Sizes file, {chrom_sizes_file}, not found or is `None`"
        )
    try:
        chrom_names = pd.read_csv(chrom_sizes_file, sep="\t", header=None)[0]
        chrom_sizes = pd.read_csv(chrom_sizes_file, sep="\t", header=None)[1]
    except Exception as e:
        logger.info(
            f"Error reading chromosome sizes file: {chrom_sizes_file}.\
            \nExpected format is tab-separated with two columns, e.g.,\
            \nchr1\t248956422\nchr2\t242193529\n..."
        )
        raise
    return dict(zip(chrom_names, chrom_sizes))


def get_bam_chrom_reads(
    bam_file: str,
    chromosome: str,
    chrom_sizes_file: str,
    step: int,
    effective_genome_size: float = -1,
    norm_method: str = "RPGC",
    min_mapping_score: int = 10,
    flag_include: int | None = None,
    flag_exclude: int = 3844,
    extend_reads: int = -1,
    center_reads: bool = False,
    ignore_for_norm: list | None = None,
    scale_factor: float = 1.0,
    num_processors: int = -1,
    const_scale: float = 1.0,
    round_digits: int = 5,
    scale_by_step: bool = False,
):
    r"""Count signal directly from a BAM file for one chromosome."""

    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if not os.path.exists(chrom_sizes_file):
        raise FileNotFoundError(f"Chromosome sizes file not found: {chrom_sizes_file}")

    native = _require_native_counter()
    chrom_sizes_dict = get_chroms_and_sizes(chrom_sizes_file)
    if chromosome not in chrom_sizes_dict:
        raise ValueError(
            f"Chromosome {chromosome} not found in chromosome sizes file: {chrom_sizes_file}"
        )

    if ignore_for_norm is None:
        ignore_for_norm = ["chrX", "chrY", "chrM"]

    chrom_size = int(chrom_sizes_dict[chromosome])
    metadata = _get_bam_count_metadata(
        bam_file,
        step=step,
        norm_method=norm_method,
        effective_genome_size=effective_genome_size,
        ignore_for_norm=ignore_for_norm,
        flag_exclude=flag_exclude,
        extend_reads=extend_reads,
        num_processors=num_processors,
        scale_factor=scale_factor,
    )
    threads = int(metadata["threads"])

    try:
        chrom_start, chrom_end = native.get_alignment_chrom_range(
            bam_file,
            chromosome,
            chrom_size,
            thread_count=threads,
            flag_exclude=max(0, int(flag_exclude)),
        )
    except RuntimeError as exc:
        if "chromosome not found" in str(exc).lower():
            logger.warning(
                "Chromosome %s not found in BAM file: %s. Returning (None,None).",
                chromosome,
                bam_file,
            )
            return None, None
        raise

    chrom_start = int(chrom_start)
    chrom_end = int(chrom_end)
    if chrom_end <= chrom_start:
        logger.warning(
            "No mapped reads found in BAM file: %s for chromosome: %s. Returning (None,None).",
            bam_file,
            chromosome,
        )
        return None, None

    count_start = max(0, (chrom_start // step) * step)
    count_end = min(
        chrom_size,
        int(np.ceil(max(chrom_end, count_start + 1) / float(step)) * step),
    )
    if count_end <= count_start:
        count_end = min(chrom_size, count_start + step)

    counts = native.count_alignment_region(
        bam_file,
        chromosome,
        count_start,
        count_end,
        int(step),
        int(metadata["read_length"]),
        one_read_per_bin=1 if center_reads else 0,
        thread_count=threads,
        flag_include=max(0, int(flag_include or 0)),
        flag_exclude=max(0, int(flag_exclude)),
        extend_bp=max(0, int(metadata["resolved_extend_bp"])),
        paired_end_mode=1 if bool(metadata["paired_end_mode"]) else 0,
        min_mapping_quality=max(0, int(min_mapping_score)),
        count_mode="coverage",
    )

    vals = np.asarray(counts, dtype=np.float64)
    intervals = count_start + (np.arange(vals.size, dtype=np.int64) * int(step))

    vals = vals * float(metadata["norm_scale"])
    if scale_by_step:
        vals = vals / float(step)
        logger.info(f"Dividing `vals` by step size (bp): {step}")

    if const_scale >= 0:
        if const_scale == 0:
            logger.warning("You are scaling the values by 0.")
        vals = vals * const_scale

    positive_idx = np.flatnonzero(vals > 0.0)
    if positive_idx.size == 0:
        logger.warning(
            "No non-zero values found in BAM file: %s for chromosome: %s. Returning (None,None).",
            bam_file,
            chromosome,
        )
        return None, None

    first_idx = int(positive_idx[0])
    last_idx = int(positive_idx[-1]) + 1
    intervals = intervals[first_idx:last_idx]
    vals = np.round(vals[first_idx:last_idx], round_digits)
    return intervals.astype(int), vals


def generate_chrom_matrix(
    chromosome: str,
    input_files: list,
    chrom_sizes_file: str,
    step: int,
    const_scale: float = 1.0,
    round_digits: int = 5,
    scale_by_step: bool = False,
    effective_genome_size: float = -1,
    norm_method: str = "RPGC",
    min_mapping_score: int = 10,
    flag_include: int | None = None,
    flag_exclude: int = 3844,
    extend_reads: int = -1,
    center_reads: bool = False,
    ignore_for_norm: list | None = None,
    scale_factor: float = 1.0,
    num_processors: int = -1,
    low_memory: bool = False,
):
    r"""Create a matrix of read counts for a given chromosome from BAM files."""

    interval_matrix = []
    vals_matrix = []
    count_processes, bam_threads = _resolve_parallel_bam_counting(
        len(input_files),
        num_processors,
    )
    count_args = []
    for input_file in input_files:
        _get_track_type(input_file)
        count_args.append(
            (
                input_file,
                chromosome,
                chrom_sizes_file,
                step,
                effective_genome_size,
                norm_method,
                min_mapping_score,
                flag_include,
                flag_exclude,
                extend_reads,
                center_reads,
                ignore_for_norm,
                scale_factor,
                bam_threads,
                const_scale,
                round_digits,
                scale_by_step,
            )
        )

    if count_processes > 1:
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(processes=count_processes) as pool:
            count_results = pool.starmap(get_bam_chrom_reads, count_args)
    else:
        count_results = [get_bam_chrom_reads(*args_) for args_ in count_args]

    for input_file, (intervals_, vals_) in zip(input_files, count_results):
        if intervals_ is None or vals_ is None:
            logger.warning(
                f"No data found for {input_file} in chromosome {chromosome}. Excluding this track for {chromosome}."
            )
            continue
        interval_matrix.append(intervals_)
        vals_matrix.append(vals_)
    if len(interval_matrix) == 0:
        logger.warning(
            f"No data found in the files {str(input_files)} for chromosome {chromosome}. Returning (None,None)."
        )
        return None, None
    common_intervals = np.sort(np.unique(np.concatenate(interval_matrix, axis=0)))
    matrix_dtype = np.float32 if low_memory else np.float64
    count_matrix = np.zeros(
        (len(interval_matrix), len(common_intervals)),
        dtype=matrix_dtype,
    )
    for i, (intervals_, vals_) in enumerate(zip(interval_matrix, vals_matrix)):
        idx = np.searchsorted(common_intervals, intervals_)
        count_matrix[i, idx] = np.asarray(vals_, dtype=matrix_dtype)

    if get_shape(count_matrix)[0] == 1:
        count_matrix = count_matrix.reshape(1, -1)

    return np.array(common_intervals).astype(int), count_matrix


def check_type_bam_files(bam_files):
    if isinstance(bam_files, str):
        with open(bam_files, "r") as f:
            bam_files_ = [line.strip() for line in f if line.strip()]
        for file_ in bam_files_:
            if not os.path.exists(file_):
                raise FileNotFoundError(f"File in `bam_files_` not found: {file_}")
    elif isinstance(bam_files, list):
        bam_files_ = bam_files
        for file_ in bam_files_:
            if not os.path.exists(file_):
                raise FileNotFoundError(f"File in `bam_files` not found: {file_}")
    else:
        raise ValueError(
            "`bam_files` must be either a list or a path to a text file containing a list of BAM file paths."
        )
    return bam_files_
