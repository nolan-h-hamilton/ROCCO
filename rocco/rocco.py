#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
==========================================================================
ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization
==========================================================================

Run ROCCO on BAM or bigWig files.
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
import multiprocessing as mp
import os
import sys
import tempfile
import uuid
from pprint import pformat
from typing import Tuple

import numpy as np
import scipy.stats as stats

from rocco.constants import GENOME_DICT
from rocco.dp import solve_chrom_exact
from rocco.inference import (
    estimate_empirical_bayes_budgets,
    estimate_budget_nonnull_fraction_from_wild_bootstrap_null,
    estimate_budget_nonnull_fraction_from_score_track,
    score_loci_wls,
)
from rocco._version import __version__
from rocco.readtracks import *
import rocco.scores as posthoc_scores


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s",
)
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)
_CHROM_SOLVE_PROCESS_STATE: dict | None = None


def _read_bed_records(
    bed_file: str,
) -> tuple[list[tuple[str, int, int]], bool]:
    records: list[tuple[str, int, int]] = []
    saw_extra_columns = False
    with open(bed_file, "r", encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, start=1):
            line_ = line.strip()
            if line_ == "":
                continue
            fields = line_.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"BED row {line_num} in {bed_file} has fewer than 3 columns."
                )
            if len(fields) > 3:
                saw_extra_columns = True
            records.append((str(fields[0]), int(fields[1]), int(fields[2])))
    return records, saw_extra_columns


def _merge_bed_records(
    records: list[tuple[str, int, int]],
    min_length_bp: int | None = None,
) -> list[tuple[str, int, int]]:
    if len(records) == 0:
        return []
    merged: list[list[str | int]] = []
    for chrom, start, end in sorted(records, key=lambda x: (x[0], x[1], x[2])):
        if len(merged) == 0:
            merged.append([chrom, int(start), int(end)])
            continue
        prev = merged[-1]
        if chrom == prev[0] and int(start) <= int(prev[2]):
            prev[2] = max(int(prev[2]), int(end))
            continue
        merged.append([chrom, int(start), int(end)])
    merged_records = [
        (str(chrom), int(start), int(end))
        for chrom, start, end in merged
        if min_length_bp is None or (int(end) - int(start)) >= int(min_length_bp)
    ]
    return merged_records


def _write_bed_records(
    records: list[tuple[str, int, int]],
    output_file: str,
    name_features: bool = False,
) -> str:
    with open(output_file, "w", encoding="utf-8") as handle:
        for chrom, start, end in records:
            if name_features:
                feature_name = f"{chrom}_{start}_{end}"
                handle.write(f"{chrom}\t{start}\t{end}\t{feature_name}\n")
            else:
                handle.write(f"{chrom}\t{start}\t{end}\n")
    return output_file


def _get_input_type(input_file: str) -> str:
    r"""Determine the supported ROCCO input type from the filename extension.

    The file type is determined by the file extension: '.bam', etc. and is not robust
    to incorrectly labelled files.

    :raises ValueError: If file extension is not supported
    :return: a string representing the file type
    :rtype: str

    """
    file_type = None
    file_ext = str(os.path.splitext(input_file.lower())[1][1:]).lower()

    if file_ext in ["bam"]:
        file_type = "bam"
    elif file_ext in ["bw", "bigwig"]:
        file_type = "bigwig"
    elif file_ext in ["bed", "bedgraph", "bg", "wig", "wiggle"]:
        bedgraph_notice = "\nBedGraph and wiggle-like inputs are not supported. Input files must be BAM alignments or bigWig tracks.\n"
        raise ValueError(bedgraph_notice)
    if file_type is None:
        raise ValueError("Input file must be a BAM alignment file or bigWig track")
    return file_type


def chrom_solution_to_bed(
    chromosome,
    intervals,
    solution,
    ID=None,
    check_gaps_intervals=True,
    min_length_bp=None,
) -> str:
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
        raise ValueError(
            f"Intervals and solution must have the same length at the pre-merge stage: {len(intervals)} != {len(solution)}"
        )

    if check_gaps_intervals:
        if len(set(np.diff(intervals))) > 1:
            raise ValueError(f"Intervals must be contiguous: {set(np.diff(intervals))}")
    step_ = intervals[1] - intervals[0]
    if ID is None:
        output_file = f"rocco_{chromosome}.bed"
    else:
        output_file = f"rocco_{ID}_{chromosome}.bed"

    selected_records: list[tuple[str, int, int]] = []
    for i in range(len(intervals) - 1):
        # At this point, solutions should be binary. Keep a 0.50 cutoff
        # here so tied or float-valued solutions still behave sensibly.
        if solution[i] > 0.50:
            selected_records.append(
                (str(chromosome), int(intervals[i]), int(intervals[i + 1]))
            )
    merged_records = _merge_bed_records(
        selected_records,
        min_length_bp=min_length_bp,
    )
    return _write_bed_records(merged_records, output_file)


def combine_chrom_results(
    chrom_bed_files: list,
    output_file: str,
    name_features: bool = False,
) -> str:
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
        logger.info(f"Removing existing output file: {output_file}")
        try:
            os.remove(output_file)
        except:
            logger.info(f"Could not remove existing output file: {output_file}.")
    combined_records: list[tuple[str, int, int]] = []
    for chrom_bed_file in chrom_bed_files:
        if not os.path.exists(chrom_bed_file):
            raise FileNotFoundError(f"File does not exist: {chrom_bed_file}")
        try:
            chrom_records, saw_extra_columns = _read_bed_records(chrom_bed_file)
        except Exception as e:
            logger.info(f"Could not read BED file: {chrom_bed_file}\n{e}\n")
            raise
        if saw_extra_columns and not printed_colct_msg:
            logger.info(
                "More than 3 columns detected in the input BED files. Extra columns will be ignored."
            )
            printed_colct_msg = True
        combined_records.extend(chrom_records)
    merged_records = _merge_bed_records(combined_records)
    return _write_bed_records(
        merged_records,
        output_file,
        name_features=name_features,
    )


def score_central_tendency_chrom(
    chrom_matrix,
    method="quantile",
    quantile=0.50,
    tprop=0.05,
    power=1.0,
) -> np.ndarray:
    r"""Return a column-wise location summary across samples."""
    chrom_matrix = np.asarray(chrom_matrix, dtype=float)
    if chrom_matrix.ndim != 2:
        raise ValueError("`chrom_matrix` must be a 2D array.")
    if chrom_matrix.shape[0] == 1:
        return np.power(chrom_matrix[0, :], power)

    method_ = str(method).strip().lower().replace("-", "").replace("_", "")
    central_tendency = None

    if method_ == "quantile":
        if not 0.0 <= quantile <= 1.0:
            logger.warning("`quantile` must be in [0, 1]. Using the median instead.")
            quantile = 0.50
        if quantile == 0.50:
            central_tendency = np.median(chrom_matrix, axis=0)
        else:
            central_tendency = np.quantile(
                chrom_matrix,
                quantile,
                axis=0,
                method="nearest",
            )
    elif method_ == "tmean":
        # Trim each column to :math:`[q_{\alpha}, q_{1-\alpha}]` before averaging.
        lower_limit = np.quantile(
            chrom_matrix,
            tprop,
            axis=0,
            method="nearest",
        )
        upper_limit = np.quantile(
            chrom_matrix,
            1.0 - tprop,
            axis=0,
            method="nearest",
        )
        central_tendency = np.array(
            [
                stats.tmean(
                    chrom_matrix[:, idx],
                    limits=(lower_limit[idx], upper_limit[idx]),
                    inclusive=(True, True),
                )
                for idx in range(chrom_matrix.shape[1])
            ],
            dtype=float,
        )
    elif method_ == "mean":
        central_tendency = np.mean(chrom_matrix, axis=0)

    if central_tendency is None:
        raise ValueError(f"Central tendency method not recognized: {method}")

    return np.power(central_tendency, power)


def score_dispersion_chrom(
    chrom_matrix: np.ndarray,
    method: str = "mad",
    rng: Tuple[int, int] = (25, 75),
    tprop: float = 0.05,
    power: float = 1.0,
) -> np.ndarray:
    r"""Return a column-wise dispersion summary across samples."""
    chrom_matrix = np.asarray(chrom_matrix, dtype=float)
    if chrom_matrix.ndim != 2:
        raise ValueError("`chrom_matrix` must be a 2D array.")
    if chrom_matrix.shape[0] == 1:
        return np.power(np.zeros_like(chrom_matrix[0, :]), power)

    method_ = str(method).strip().lower().replace("-", "").replace("_", "")
    dispersion = None

    if method_ == "mad":
        dispersion = stats.median_abs_deviation(chrom_matrix, axis=0)
    elif method_ == "iqr":
        dispersion = stats.iqr(chrom_matrix, rng=rng, axis=0)
    elif method_ == "std":
        dispersion = np.std(chrom_matrix, axis=0)
    elif method_ == "tstd":
        lower_limit = np.quantile(
            chrom_matrix,
            tprop,
            axis=0,
            method="nearest",
        )
        upper_limit = np.quantile(
            chrom_matrix,
            1.0 - tprop,
            axis=0,
            method="nearest",
        )
        dispersion = stats.tstd(
            chrom_matrix,
            limits=(lower_limit, upper_limit),
            inclusive=(True, True),
            axis=0,
        )

    if dispersion is None:
        raise ValueError(
            f"Dispersion method not recognized or could not execute: {method}"
        )

    return np.power(dispersion, power)


def cscores_quantiles(
    chrom_scores: np.ndarray,
    quantiles: np.ndarray = None,
    add_newlines=True,
) -> str:
    """Return a formatted string of quantiles for a locus-score array.

    :param chrom_scores: locus scores (float) for a given chromosome
    :type chrom_scores: np.ndarray
    :param quantiles: array of quantiles in [0.0,1.0] to compute.
    :type quantiles: np.ndarray, optional
    :return: pformatted string of quantiles
    :rtype: str
    """
    if quantiles is None:
        quantiles = np.array(
            [
                0.0,
                0.01,
                0.05,
                0.25,
                0.50,
                0.75,
                0.95,
                0.975,
                0.99,
                1.0,
            ]
        )
    formatted_string = pformat(
        {
            f"Quantile={q}": round(np.quantile(chrom_scores, q=q, method="higher"), 4)
            for q in quantiles
        }
    )
    if add_newlines:
        return f"\n{formatted_string}\n"
    return f"{formatted_string}"


def json_config(config_path):
    with open(config_path, "r") as json_file:
        return json.load(json_file)


def resolve_config(args):
    """Resolve command-line arguments with a JSON configuration file

    :param args: Command-line arguments obtained with `argparse`
    :type args: dict
    :return: Resolved command-line arguments
    :rtype: dict

    .. note::

        * Modifies/overrides command-line arguments specified explicitly in the JSON file
        * For boolean arguments, use `true` or `false` in the JSON file rather than `True` or `False`

    **Example JSON config file**

    .. code-block:: json

        {
            "input_files": ["sample1.bam", "sample2.bam"],
            "output": "rocco_peaks_output.bed",
            "genome": "hg38",
            "chroms": ["chr21", "chr22"],
            "int_tol": 0.01,
            "verbose": true
        }

    Can then run `rocco --config config.json [...]`.
    """

    args_ = copy.deepcopy(args)
    if args_["config"] is None or not os.path.exists(args_["config"]):
        return args_

    json_args = json_config(args_["config"])
    for key, value in json_args.items():
        if key not in args_.keys():
            continue
        args_[key] = value
        logger.info(f"Setting {key}={value} per {args_['config']}")
    return args_


def _build_parser() -> argparse.ArgumentParser:
    epilog_cli_help = (
        "\nGitHub (Homepage): <https://github.com/nolan-h-hamilton/ROCCO/>\n"
        "Paper: <https://doi.org/10.1093/bioinformatics/btad725>\n"
    )
    parser = argparse.ArgumentParser(
        description="ROCCO Consensus Peak Detection Algorithm for Multisample HTS Datasets",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog_cli_help,
    )
    parser.add_argument(
        "--input_files",
        "-i",
        nargs="+",
        help="BAM alignment files or pre-scored bigWig tracks corresponding to samples",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"rocco {__version__}",
    )
    parser.add_argument(
        "--output",
        "--outfile",
        "-o",
        type=str,
        default=f"rocco_peaks_output_{str(int(uuid.uuid4().hex[:5], base=16))}.bed",
    )
    parser.add_argument(
        "--genome",
        "-g",
        default=None,
        help="Genome assembly. Invoking this argument with a supported assembly (hg38, hg19, mm10, mm39, dm6) will use default `--chrom_sizes_file` and `--effective_genome_size` values specific to that assembly. If this argument is not invoked, you can supply those arguments manually.",
    )
    parser.add_argument(
        "--chrom_sizes_file",
        "-s",
        default=None,
        help="Chromosome sizes file. Required if genome is not specified",
    )
    parser.add_argument(
        "--effective_genome_size",
        type=int,
        default=None,
        help="Effective genome size. Required if genome is not specified and using RPGC normalization",
    )
    parser.add_argument(
        "--chroms",
        nargs="+",
        type=str,
        default=[],
        help="Chromosomes to process. If not specified, all chromosomes will be processed",
    )
    parser.add_argument(
        "--skip_chroms",
        nargs="+",
        type=str,
        default=[],
        help="Chromosomes to skip altogether -- peaks in these chromosomes are effectively omitted",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Invoke for verbose output",
    )

    parser.add_argument(
        "--budget",
        type=float,
        default=None,
        help="Upper bounds the proportion of the genome that can be selected as open chromatin.",
    )
    parser.add_argument(
        "--budget_null_draws",
        type=int,
        default=25,
        help="Maximum number of null draws used when initializing chromosome budgets. Default is 25.",
    )
    parser.add_argument(
        "--scale_chrom_budgets",
        type=float,
        default=1.0,
        help="This constant scales each chromosome-specific budget.",
    )
    parser.add_argument(
        "--budget_posterior_quantile",
        type=float,
        default=0.01,
        help="Lower beta-posterior quantile used to summarize EB chromosome budgets. Smaller values are more conservative.",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=None,
        help="Boundary penalty used by the exact DP. If omitted, ROCCO derives a chromosome-level default from score scale and autocorrelation.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=-1,
        help="Number of threads to use. Default is -1 (use available cores)",
    )
    parser.add_argument(
        "--low_memory",
        action="store_true",
        default=False,
        help="Use less memory by keeping large temporary matrices in single precision and using fewer processes automatically.",
    )
    parser.add_argument(
        "--selection_penalty",
        type=float,
        default=None,
        help="Direct penalty on selected loci. If supplied, `--budget` is ignored.",
    )

    parser.add_argument(
        "--score_lower_bound_z",
        type=float,
        default=1.0,
        help="Direct floor subtracted from the moderated standardized score when `--score_min_effect` is not supplied.",
    )
    parser.add_argument(
        "--score_prior_df",
        type=float,
        default=6.0,
        help="Prior degrees of freedom for EB variance shrinkage.",
    )
    parser.add_argument(
        "--score_min_effect",
        type=float,
        default=None,
        help="Optional minimum effect on the log2(1+x) scale used in the optimization score. If omitted, ROCCO optimizes the moderated standardized score shifted by `--score_lower_bound_z`.",
    )
    parser.add_argument(
        "--score_precision_floor_ratio",
        type=float,
        default=0.01,
        help="Lower bound on moderated variance as a fraction of the prior variance before computing standard errors.",
    )

    parser.add_argument(
        "--step",
        "-w",
        type=int,
        default=50,
        help="Bin width used for BAM inputs. Ignored for bigWig inputs, which use their native binning scheme.",
    )
    parser.add_argument(
        "--norm_method",
        default="RPGC",
        choices=["RPGC", "CPM", "RPKM", "BPM", "rpgc", "cpm", "rpkm", "bpm"],
        help="Normalization method for BAM inputs. Default is RPGC (Reads Per Genomic Content), for which the `--effective_genome_size` argument is required (default EGS values supplied automatically for supported genomes).",
    )
    parser.add_argument(
        "--min_mapping_score",
        type=int,
        default=10,
        help="Equivalent to samtools view -q.",
    )
    parser.add_argument(
        "--flag_include",
        type=int,
        default=-1,
        help="Equivalent to samtools view -f.",
    )
    parser.add_argument(
        "--flag_exclude",
        type=int,
        default=3844,
        help="Equivalent to samtools view -F.",
    )
    parser.add_argument(
        "--extend_reads",
        type=int,
        default=-1,
        help="If < 0, reads are not extended. If 0, fragment length is estimated from the alignments. If > 0, reads are extended by that number of base pairs.",
    )
    parser.add_argument(
        "--center_reads",
        action="store_true",
        help="Center reads at the midpoint of the fragment length.",
    )
    parser.add_argument(
        "--ignore_for_norm",
        nargs="+",
        default=[],
        help="Chromosomes to ignore for normalization.",
    )
    parser.add_argument(
        "--scale_factor",
        type=float,
        default=1.0,
        help="Scale factor applied after normalization.",
    )
    parser.add_argument("--round_digits", type=int, default=5)

    parser.add_argument("--min_length_bp", type=int, default=None)
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--narrowPeak", action="store_true", default=False)
    parser.add_argument(
        "--ecdf_samples",
        type=int,
        default=250,
        help="Number of background regions sampled for each representative length bin in the empirical null.",
    )
    parser.add_argument("--ecdf_seed", type=int, default=42)
    parser.add_argument(
        "--ecdf_proc",
        type=int,
        default=None,
        dest="ecdf_proc",
        help="Number of processes used to fit the binned empirical null for `--narrowPeak`.",
    )
    return parser


def _prepare_args(parser: argparse.ArgumentParser) -> dict:
    args = vars(parser.parse_args())
    args = resolve_config(args)
    if (
        len(sys.argv) == 1
        or args["input_files"] is None
        or len(args["input_files"]) == 0
    ):
        parser.print_help(sys.stdout)
        sys.exit(0)

    args["norm_method"] = clean_string(args["norm_method"]).upper()
    if args["low_memory"]:
        if int(args["threads"]) <= 0:
            total_cores = max(1, os.cpu_count() or 1)
            args["threads"] = int(min(4, max(1, total_cores // 4)))
        if "--budget_null_draws" not in sys.argv and int(
            args["budget_null_draws"]
        ) == int(parser.get_default("budget_null_draws")):
            args["budget_null_draws"] = 16

    if args["genome"] is not None:
        args["genome"] = clean_string(args["genome"])
        if args["genome"] not in GENOME_DICT:
            raise ValueError(
                f"Genome not found: {args['genome']}. Available genomes: {list(GENOME_DICT.keys())}"
            )
        if args["effective_genome_size"] is None:
            args["effective_genome_size"] = GENOME_DICT[args["genome"]][
                "effective_genome_size"
            ]
        if args["chrom_sizes_file"] is None:
            args["chrom_sizes_file"] = GENOME_DICT[args["genome"]]["sizes_file"]

    input_types = {_get_input_type(file_) for file_ in args["input_files"]}
    if len(input_types) != 1:
        raise ValueError("All input files must share the same type.")
    args["input_track_type"] = next(iter(input_types))

    if args["chrom_sizes_file"] is None:
        raise ValueError(
            "A chromosome sizes file must be supplied with `-s/--chrom_sizes_file` when genome defaults are unavailable."
        )
    if (
        args["input_track_type"] == "bam"
        and args["effective_genome_size"] is None
        and args["norm_method"] == "RPGC"
    ):
        raise ValueError(
            "`--effective_genome_size` is required when using `--norm_method RPGC` without genome defaults."
        )
    return args


def _prepare_inputs(args: dict) -> tuple[list, list]:
    bam_files = []
    signal_inputs = []
    expected_input_type = str(args.get("input_track_type", "")).strip().lower()
    for file_ in args["input_files"]:
        if not os.path.exists(file_):
            raise FileNotFoundError(f"File not found: {file_}")
        input_type = _get_input_type(file_)
        if expected_input_type != "" and input_type != expected_input_type:
            raise ValueError("All input files must share the same type.")
        if input_type == "bam":
            bam_files.append(file_)
        signal_inputs.append(file_)

    if args["ignore_for_norm"] is None or len(args["ignore_for_norm"]) == 0:
        args["ignore_for_norm"] = ["chrX", "chrY", "chrM"]

    return bam_files, signal_inputs


def _resolve_chromosomes(args: dict) -> list:
    chroms_to_process = list(get_chroms_and_sizes(args["chrom_sizes_file"]).keys())
    if args["chroms"]:
        chroms_to_process = [
            chrom for chrom in chroms_to_process if chrom in args["chroms"]
        ]
    if args["skip_chroms"]:
        chroms_to_process = [
            chrom for chrom in chroms_to_process if chrom not in args["skip_chroms"]
        ]
    return chroms_to_process


def _resolve_chrom_gamma(
    chrom: str,
    args: dict,
    chrom_scores: np.ndarray,
    budget_rate_meta: dict,
) -> tuple[float, dict | None]:
    if args["gamma"] is not None:
        chrom_gamma = float(args["gamma"])
        if not np.isfinite(chrom_gamma) or chrom_gamma < 0.0:
            raise ValueError("`--gamma` must be finite and non-negative")
        logger.info("%s fixed gamma value=%.6f", chrom, chrom_gamma)
        return float(chrom_gamma), None

    scores_ = np.asarray(chrom_scores, dtype=np.float64)
    positive_scores = scores_[scores_ > 0.0]
    if positive_scores.size == 0:
        positive_scale = 1.0
        positive_count = 0
    else:
        positive_scale = float(np.median(positive_scores))
        positive_count = int(positive_scores.size)
    autocorrelation_time = float(budget_rate_meta.get("autocorrelation_time", 1.0))
    autocorrelation_time = max(1.0, autocorrelation_time)
    characteristic_run = int(np.ceil(autocorrelation_time))
    gamma_raw = 0.5 * float(characteristic_run) * float(positive_scale)
    chrom_gamma = float(np.clip(gamma_raw, 0.5, 10.0))
    gamma_meta = {
        "method": "auto_score_autocorr",
        "autocorrelation_time": float(autocorrelation_time),
        "characteristic_run_length": int(characteristic_run),
        "positive_score_median": float(positive_scale),
        "positive_score_count": int(positive_count),
        "gamma_raw": float(gamma_raw),
        "gamma_clipped": float(chrom_gamma),
        "gamma_clip_min": 0.5,
        "gamma_clip_max": 10.0,
    }
    logger.info("%s auto gamma estimate: %s", chrom, gamma_meta)
    return float(chrom_gamma), gamma_meta


def _resolve_parallel_process_count(
    item_count: int,
    thread_limit: int,
) -> int:
    total_cores = (
        max(1, os.cpu_count() or 1)
        if int(thread_limit) <= 0
        else max(1, int(thread_limit))
    )
    if int(item_count) <= 1 or int(total_cores) <= 1:
        return 1
    if "fork" not in mp.get_all_start_methods():
        return 1
    # Keep this conservative. The arrays are bigger than the solve.
    return int(min(int(item_count), int(total_cores), 4))


def _cpy_narrowpeak_summit_track(
    chrom: str,
    intervals: np.ndarray,
    effect_mean: np.ndarray,
) -> str | None:
    intervals_ = np.asarray(intervals, dtype=np.int64)
    effect_mean_ = np.asarray(effect_mean, dtype=np.float32)
    usable = int(min(max(intervals_.shape[0] - 1, 0), effect_mean_.shape[0]))
    if usable <= 0:
        return None
    starts = intervals_[:usable]
    centers = (
        intervals_[:usable].astype(np.int64)
        + intervals_[1 : usable + 1].astype(np.int64)
    ) // 2
    fd, summit_track_file = tempfile.mkstemp(
        prefix=f"rocco_summit_track_{chrom}_",
        suffix=".npz",
    )
    os.close(fd)
    np.savez(
        summit_track_file,
        starts=starts.astype(np.int64, copy=False),
        centers=centers.astype(np.int64, copy=False),
        mean=effect_mean_[:usable].astype(np.float32, copy=False),
    )
    return summit_track_file


def _write_narrowpeak_summit_offsets(
    peak_file: str,
    chrom_cache: dict,
    output_file: str,
) -> str:
    records, _ = _read_bed_records(peak_file)
    loaded_tracks: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    with open(output_file, "w", encoding="utf-8") as handle:
        for chrom, start, end in records:
            peak_name = f"{chrom}_{start}_{end}"
            summit_offset = -1
            chrom_data = chrom_cache.get(chrom, {})
            summit_track_file = chrom_data.get("summit_track_file")
            peak_length = int(end) - int(start)
            if summit_track_file is not None and peak_length > 0:
                if chrom not in loaded_tracks:
                    with np.load(summit_track_file) as summit_track:
                        loaded_tracks[chrom] = (
                            np.asarray(summit_track["starts"], dtype=np.int64),
                            np.asarray(summit_track["centers"], dtype=np.int64),
                            np.asarray(summit_track["mean"], dtype=np.float64),
                        )
                starts, centers, mean_track = loaded_tracks[chrom]
                left = int(np.searchsorted(starts, int(start), side="left"))
                right = int(np.searchsorted(starts, int(end), side="left"))
                if right > left:
                    local_mean = mean_track[left:right]
                    if np.any(np.isfinite(local_mean)):
                        local_idx = int(np.nanargmax(local_mean))
                        summit_bp = int(centers[left + local_idx])
                        summit_offset = int(
                            np.clip(summit_bp - int(start), 0, max(peak_length - 1, 0))
                        )
            handle.write(f"{peak_name}\t{summit_offset}\n")
    return output_file


def _cleanup_narrowpeak_tempfiles(chrom_cache: dict):
    for chrom_data in chrom_cache.values():
        summit_track_file = chrom_data.pop("summit_track_file", None)
        if summit_track_file is None:
            continue
        try:
            os.remove(summit_track_file)
        except Exception as exc:
            logger.info(
                "Could not remove narrowPeak summit temp. file %s\n%s",
                summit_track_file,
                exc,
            )


def _solve_cached_chromosome(chrom_: str) -> tuple[str, float, dict, str]:
    state = _CHROM_SOLVE_PROCESS_STATE
    if state is None:
        raise RuntimeError("Chromosome solve state is not initialized")
    chrom_data = state["chrom_cache"][chrom_]
    chrom_budget = state["chrom_budgets"][chrom_]
    chrom_gamma = chrom_data["gamma"]
    if not np.all(np.isfinite(chrom_data["scores"])):
        raise ValueError(f"{chrom_} scores contain non-finite values")
    try:
        chrom_budget = float(chrom_budget)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{chrom_} budget could not be read as a finite number"
        ) from exc
    try:
        chrom_gamma = float(chrom_gamma)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{chrom_} gamma could not be read as a finite number"
        ) from exc
    if not np.isfinite(chrom_budget) or chrom_budget < 0.0:
        raise ValueError(f"{chrom_} budget must be finite and non-negative")
    if not np.isfinite(chrom_gamma) or chrom_gamma < 0.0:
        raise ValueError(f"{chrom_} gamma must be finite and non-negative")
    chrom_sol, chrom_obj, chrom_meta = solve_chrom_exact(
        chrom_data["scores"],
        budget=chrom_budget,
        gamma=chrom_gamma,
        selection_penalty=state["selection_penalty"],
        return_details=True,
    )
    chrom_outfile = chrom_solution_to_bed(
        chrom_,
        chrom_data["intervals"],
        chrom_sol,
        state["run_id"],
        check_gaps_intervals=True,
        min_length_bp=state["min_length_bp"],
    )
    return chrom_, float(chrom_obj), chrom_meta, chrom_outfile


def _build_chrom_cache(
    chroms_to_process: list,
    signal_inputs: list,
    args: dict,
) -> dict:
    chrom_cache = {}
    low_memory = bool(args.get("low_memory", False))
    budget_null_processes = (
        1
        if low_memory
        else _resolve_parallel_process_count(
            int(args["budget_null_draws"]),
            int(args["threads"]),
        )
    )
    for chrom_ in chroms_to_process:
        logger.info("Generating chromosome matrix: %s", chrom_)
        chrom_intervals, chrom_matrix = generate_chrom_matrix(
            chrom_,
            signal_inputs,
            args["chrom_sizes_file"],
            args["step"],
            round_digits=args["round_digits"],
            effective_genome_size=args["effective_genome_size"],
            norm_method=args["norm_method"],
            min_mapping_score=args["min_mapping_score"],
            flag_include=args["flag_include"],
            flag_exclude=args["flag_exclude"],
            extend_reads=args["extend_reads"],
            center_reads=args["center_reads"],
            ignore_for_norm=args["ignore_for_norm"],
            scale_factor=args["scale_factor"],
            num_processors=args["threads"],
            low_memory=low_memory,
        )

        if chrom_intervals is None or chrom_matrix is None:
            logger.warning("Skipping chromosome %s... no data found.", chrom_)
            continue

        logger.info("Chromosome %s matrix: %s", chrom_, chrom_matrix.shape)
        if not np.all(np.isfinite(chrom_matrix)):
            raise ValueError(f"{chrom_} matrix contains non-finite values")
        # skip WLS, note that multiple bigwigs really should _not_ be supplied unless it makes sense to aggregate via central tendency
        if args["input_track_type"] == "bigwig":
            if chrom_matrix.shape[0] > 1:
                logger.warning(
                    "Multiple bigwig tracks detected for %s. ROCCO will aggregate these via column-wise central tendency rather than WLS. If this is not the intended behavior, please supply a single bigwig track per chromosome or use BAM inputs.",
                    chrom_,
                )
            chrom_scores = np.asarray(
                score_central_tendency_chrom(
                    chrom_matrix,
                    method="quantile",
                    quantile=0.50,
                    power=1.0,
                ),
                dtype=np.float64,
            )
            if not np.all(np.isfinite(chrom_scores)):
                raise ValueError(f"{chrom_} direct scores contain non-finite values")
            score_details = {
                "mean": chrom_scores.astype(np.float64, copy=False),
            }
            budget_fraction_hat, budget_rate_meta = (
                estimate_budget_nonnull_fraction_from_score_track(
                    chrom_scores,
                    num_null_draws=args["budget_null_draws"],
                    progress_label=f"Budget null {chrom_}",
                    num_processes=min(
                        int(args["budget_null_draws"]),
                        int(budget_null_processes),
                    ),
                    return_details=True,
                )
            )
        else:
            chrom_scores, score_details = score_loci_wls(
                chrom_matrix,
                lower_bound_z=args["score_lower_bound_z"],
                prior_df=args["score_prior_df"],
                min_effect=args.get("score_min_effect"),
                precision_floor_ratio=args["score_precision_floor_ratio"],
                low_memory=low_memory,
                return_details=True,
            )
            if not np.all(np.isfinite(chrom_scores)):
                raise ValueError(f"{chrom_} scores contain non-finite values")
            centered_matrix = np.asarray(
                score_details.pop("centered_matrix"),
                dtype=np.float32 if low_memory else np.float64,
            )
            if not np.all(np.isfinite(centered_matrix)):
                raise ValueError(f"{chrom_} centered matrix contains non-finite values")
            budget_fraction_hat, budget_rate_meta = (
                estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
                    centered_matrix,
                    observed_scores=chrom_scores,
                    lower_bound_z=args["score_lower_bound_z"],
                    prior_df=args["score_prior_df"],
                    min_effect=args.get("score_min_effect"),
                    precision_floor_ratio=args["score_precision_floor_ratio"],
                    dependence_lag_hint=max(
                        25,
                        int(score_details.get("local_baseline_window", 101)),
                    ),
                    num_null_draws=args["budget_null_draws"],
                    progress_label=f"Budget null {chrom_}",
                    num_processes=min(
                        int(args["budget_null_draws"]),
                        int(budget_null_processes),
                    ),
                    return_details=True,
                )
            )
            del centered_matrix
        del chrom_matrix
        if not np.isfinite(budget_fraction_hat):
            raise ValueError(f"{chrom_} budget estimate is not finite")
        budget_total_count_hat = float(
            np.clip(
                budget_rate_meta.get("effective_total_count", chrom_scores.shape[0]),
                1.0,
                chrom_scores.shape[0],
            )
        )
        budget_count_hat = float(
            np.clip(
                budget_fraction_hat * budget_total_count_hat,
                0.0,
                budget_total_count_hat,
            )
        )
        score_label = (
            "direct input scores"
            if args["input_track_type"] == "bigwig"
            else "WLS scores"
        )
        logger.info("%s %s:%s", chrom_, score_label, cscores_quantiles(chrom_scores))
        logger.info(
            "%s raw budget estimate: %s",
            chrom_,
            budget_rate_meta,
        )
        chrom_gamma, gamma_meta = _resolve_chrom_gamma(
            chrom_,
            args,
            chrom_scores,
            budget_rate_meta,
        )
        chrom_effect_mean = np.asarray(
            score_details.get("mean", chrom_scores),
            dtype=np.float64,
        )
        chrom_cache[chrom_] = {
            "intervals": chrom_intervals,
            "scores": chrom_scores,
            "effect_mean": chrom_effect_mean,
            "gamma": chrom_gamma,
            "gamma_meta": gamma_meta,
            "budget_count_hat": float(budget_count_hat),
            "budget_fraction_hat": float(budget_fraction_hat),
            "budget_rate_meta": budget_rate_meta,
            "total_count": float(budget_total_count_hat),
            "num_loci": int(chrom_scores.shape[0]),
        }

    if args.get("narrowPeak", False) and args["input_track_type"] == "bam":
        for chrom_, chrom_data in chrom_cache.items():
            chrom_data["summit_track_file"] = _cpy_narrowpeak_summit_track(
                chrom_,
                chrom_data["intervals"],
                chrom_data["effect_mean"],
            )

    for chrom_data in chrom_cache.values():
        chrom_data.pop("effect_mean", None)
    return chrom_cache


def _resolve_budgets(
    chrom_cache: dict,
    args: dict,
) -> tuple[dict, dict]:
    chrom_budget_counts = {
        chrom: chrom_cache[chrom]["budget_count_hat"] for chrom in chrom_cache
    }
    chrom_total_counts = {
        chrom: chrom_cache[chrom]["total_count"] for chrom in chrom_cache
    }
    chrom_budgets, budget_meta = estimate_empirical_bayes_budgets(
        chrom_budget_counts,
        chrom_total_counts,
        posterior_quantile=args["budget_posterior_quantile"],
    )
    if args["budget"] is not None and budget_meta["genome_wide_budget"] > 0:
        rescale = float(args["budget"]) / budget_meta["genome_wide_budget"]
    else:
        rescale = 1.0
    chrom_budgets = {
        chrom: min(
            max(
                chrom_budgets[chrom] * rescale * float(args["scale_chrom_budgets"]),
                0.005,
            ),
            0.1,
        )
        for chrom in chrom_budgets
    }
    logger.info("Empirical-Bayes budget prior: %s", budget_meta)
    return chrom_budgets, budget_meta


def _solve_cached_chromosomes(
    chrom_cache: dict,
    chrom_budgets: dict,
    args: dict,
    run_id: str,
) -> list:
    tmp_chrom_bed_files = []
    solve_processes = _resolve_parallel_process_count(
        len(chrom_cache),
        int(args["threads"]),
    )
    for chrom_, chrom_data in chrom_cache.items():
        chrom_budget = chrom_budgets[chrom_]
        chrom_gamma = chrom_data["gamma"]
        logger.info(
            "%s: budget=%s gamma=%s",
            chrom_,
            round(chrom_budget, 6),
            round(chrom_gamma, 6),
        )
    global _CHROM_SOLVE_PROCESS_STATE
    _CHROM_SOLVE_PROCESS_STATE = {
        "chrom_cache": chrom_cache,
        "chrom_budgets": chrom_budgets,
        "selection_penalty": args["selection_penalty"],
        "min_length_bp": args["min_length_bp"],
        "run_id": run_id,
    }

    try:
        if solve_processes > 1:
            # Use fork here so the chromosome score arrays do not get copied up front.
            ctx = mp.get_context("fork")
            with ctx.Pool(processes=solve_processes) as pool:
                solve_results = pool.map(_solve_cached_chromosome, list(chrom_cache))
        else:
            solve_results = [_solve_cached_chromosome(chrom_) for chrom_ in chrom_cache]
    finally:
        _CHROM_SOLVE_PROCESS_STATE = None

    for chrom_, chrom_obj, chrom_meta, chrom_outfile in solve_results:
        logger.info(
            "%s solve: selected=%s (%.6f), selection_penalty=%.6f, objective=%.4f",
            chrom_,
            chrom_meta["selected_count"],
            chrom_meta["selected_fraction"],
            chrom_meta["selection_penalty"],
            chrom_obj,
        )
        tmp_chrom_bed_files.append(chrom_outfile)
    return tmp_chrom_bed_files


def _generate_narrowpeak_if_requested(
    args: dict,
    final_output: str,
    chrom_cache: dict,
):
    if not args["narrowPeak"]:
        return
    if args.get("input_track_type") != "bam":
        logger.info(
            "Skipping narrowPeak generation because posthoc peak scoring requires BAM inputs."
        )
        return
    summit_offsets_file = None
    try:
        output_root, output_ext = os.path.splitext(final_output)
        if output_ext.lower() == ".bed":
            sidecar_root = output_root
        else:
            sidecar_root = final_output
        narrowpeak_filepath = f"{sidecar_root}.narrowPeak"
        fd, summit_offsets_file = tempfile.mkstemp(
            prefix="rocco_pointsource_",
            suffix=".tsv",
        )
        os.close(fd)
        _write_narrowpeak_summit_offsets(
            final_output,
            chrom_cache,
            summit_offsets_file,
        )
        posthoc_scores.score_peaks(
            args["input_files"],
            args["chrom_sizes_file"],
            final_output,
            count_matrix_file=f"{sidecar_root}.counts.tsv",
            output_file=narrowpeak_filepath,
            ecdf_nsamples=args["ecdf_samples"],
            seed=args["ecdf_seed"],
            proc=args["ecdf_proc"],
            summit_offsets_file=summit_offsets_file,
        )
        logger.info("Final narrowPeak output: %s", narrowpeak_filepath)
    except Exception as exc:
        logger.info(
            "\nCould not generate narrowPeak-formatted output\n%s",
            exc,
        )
    finally:
        if summit_offsets_file is not None:
            try:
                os.remove(summit_offsets_file)
            except Exception as exc:
                logger.info(
                    "Could not remove narrowPeak pointSource temp. file %s\n%s",
                    summit_offsets_file,
                    exc,
                )


def main():
    run_id = str(int(uuid.uuid4().hex[:5], base=16))
    logger.info("\nID: %s", run_id)
    parser = _build_parser()
    args = _prepare_args(parser)

    bam_files, signal_inputs = _prepare_inputs(args)
    logger.info("Signal inputs: %s", signal_inputs)

    chroms_to_process = _resolve_chromosomes(args)
    logger.info("Chromosomes: %s", chroms_to_process)
    chrom_cache = _build_chrom_cache(
        chroms_to_process,
        signal_inputs,
        args,
    )
    chrom_budgets, _ = _resolve_budgets(chrom_cache, args)
    tmp_chrom_bed_files = _solve_cached_chromosomes(
        chrom_cache,
        chrom_budgets,
        args,
        run_id,
    )

    logger.info("Combining chromosome solutions")
    final_output = combine_chrom_results(
        tmp_chrom_bed_files,
        args["output"],
        name_features=False,
    )
    if os.path.exists(final_output):
        logger.info("Final BED output: %s", final_output)

    logger.info("Cleaning up temporary files")
    for tmp_file in tmp_chrom_bed_files:
        try:
            os.remove(tmp_file)
        except Exception as exc:
            logger.info(
                "Could not remove chromosome-specific temp. file %s\n%s",
                tmp_file,
                exc,
            )

    try:
        _generate_narrowpeak_if_requested(args, final_output, chrom_cache)
    finally:
        _cleanup_narrowpeak_tempfiles(chrom_cache)


if __name__ == "__main__":
    main()
