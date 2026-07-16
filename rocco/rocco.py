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
from collections.abc import Iterator, Mapping
from pprint import pformat

import numpy as np
import scipy.stats as stats

from rocco.constants import GENOME_DICT
from rocco.dp import solve_chrom_exact
from rocco.inference import (
    estimate_budget_nonnull_fraction_from_score_track,
    estimate_budget_nonnull_fraction_from_wild_bootstrap_null,
    estimate_empirical_bayes_budgets,
    score_loci_wls,
)
from rocco.dependence import choose_dependence_span
from rocco._version import __version__
from rocco.readtracks import (
    generate_chrom_matrix,
    get_chroms_and_sizes,
    get_track_type,
)
import rocco.scores as posthoc_scores

logger = logging.getLogger(__name__)
_CHROM_SOLVE_PROCESS_STATE: dict | None = None


class _SpoolMatrixMapping(Mapping[str, np.ndarray]):
    def __init__(self, matrix_paths: dict[str, str]):
        self._matrix_paths = dict(matrix_paths)

    def __getitem__(self, chromosome: str) -> np.ndarray:
        return np.load(
            self._matrix_paths[chromosome],
            mmap_mode="r",
            allow_pickle=False,
        )

    def __iter__(self) -> Iterator[str]:
        return iter(self._matrix_paths)

    def __len__(self) -> int:
        return len(self._matrix_paths)


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
    step_ = int(intervals[1] - intervals[0]) if len(intervals) > 1 else 1
    if ID is None:
        output_file = f"rocco_{chromosome}.bed"
    else:
        output_file = f"rocco_{ID}_{chromosome}.bed"

    selected_records: list[tuple[str, int, int]] = []
    for i in range(len(intervals)):
        # At this point, solutions should be binary. Keep a 0.50 cutoff
        # here so tied or float-valued solutions still behave sensibly.
        if solution[i] > 0.50:
            selected_records.append(
                (str(chromosome), int(intervals[i]), int(intervals[i] + step_))
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

    central_tendency = None

    if method == "quantile":
        if not 0.0 <= quantile <= 1.0:
            logger.warning("`quantile` must be in [0, 1]. Using the median instead.")
            quantile = 0.50
        if quantile == 0.50:
            central_tendency = np.nanmedian(chrom_matrix, axis=0)
        else:
            central_tendency = np.nanquantile(
                chrom_matrix,
                quantile,
                axis=0,
                method="nearest",
            )
    elif method == "tmean":
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
    elif method == "mean":
        central_tendency = np.mean(chrom_matrix, axis=0)

    if central_tendency is None:
        raise ValueError(f"Central tendency method not recognized: {method}")

    return np.power(central_tendency, power)


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
        allow_abbrev=False,
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
        default=64,
        help="Maximum number of null draws used when initializing chromosome budgets. Default is 64.",
    )
    parser.add_argument(
        "--num_null_blocks",
        type=int,
        default=4,
        help="Number of contiguous blocks used for budget initialization. A value of 1 uses chromosome-wide budgets.",
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
        default=0.1,
        help="Lower beta-posterior quantile used to summarize EB chromosome budgets. Smaller values are more conservative.",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=1.0,
        help="Boundary penalty used by the exact DP.",
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
        default=2.0,
        help="Direct floor subtracted from the moderated standardized score when `--score_min_effect` is not supplied.",
    )
    parser.add_argument(
        "--broad_score_lower_bound_z",
        type=float,
        default=1.0,
        help="Weaker score floor used to build broad peak parents.",
    )
    parser.add_argument(
        "--score_prior_df",
        type=float,
        default=10.0,
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
        default=100,
        help="Bin width used for BAM inputs. Ignored for bigWig inputs, which use their native binning scheme.",
    )
    parser.add_argument(
        "--norm_method",
        default="RPGC",
        choices=["RPGC", "CPM", "RPKM", "BPM"],
        help="Normalization method for BAM inputs. Default is RPGC (Reads Per Genomic Content), for which the `--effective_genome_size` argument is required (default EGS values supplied automatically for supported genomes).",
    )
    parser.add_argument(
        "--min_mapping_score",
        type=int,
        default=20,
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
    parser.add_argument(
        "--peak_mode",
        choices=["narrow", "broad", "both"],
        default=None,
    )
    parser.add_argument("--min_peak_score", type=float, default=0.1)
    parser.add_argument("--broad_max_gap_bp", type=int, default=None)
    parser.add_argument("--broad_min_peak_bp", type=int, default=1000)
    parser.add_argument("--window_bp", type=int, default=50000)
    parser.add_argument("--window_count", type=int, default=256)
    parser.add_argument("--working_quantile", type=float, default=0.95)
    parser.add_argument("--bootstrap_draws", type=int, default=500)
    parser.add_argument(
        "--insufficient_data_policy",
        choices=["error", "priorOnly"],
        default="error",
    )
    parser.add_argument("--prior_radius_bp", type=float, default=None)
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
        help="Number of processes used to fit the binned empirical null for peak sidecars.",
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

    if args["norm_method"] not in {"RPGC", "CPM", "RPKM", "BPM"}:
        raise ValueError(
            f"`--norm_method` must be one of RPGC, CPM, RPKM, or BPM, not `{args['norm_method']}`"
        )
    if args["low_memory"]:
        if int(args["threads"]) <= 0:
            total_cores = max(1, os.cpu_count() or 1)
            args["threads"] = int(min(4, max(1, total_cores // 4)))
        if "--budget_null_draws" not in sys.argv and int(
            args["budget_null_draws"]
        ) == int(parser.get_default("budget_null_draws")):
            args["budget_null_draws"] = 16

    args["num_null_blocks"] = int(args["num_null_blocks"])
    if args["num_null_blocks"] <= 0:
        raise ValueError("`--num_null_blocks` must be positive")

    if args["genome"] is not None:
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

    input_types = {get_track_type(file_) for file_ in args["input_files"]}
    if len(input_types) != 1:
        raise ValueError("All input files must share the same type.")
    args["input_track_type"] = next(iter(input_types))
    if not np.isfinite(float(args["gamma"])) or float(args["gamma"]) < 0.0:
        raise ValueError("`--gamma` must be finite and non-negative")
    if (
        not np.isfinite(float(args["score_lower_bound_z"]))
        or float(args["score_lower_bound_z"]) < 0.0
    ):
        raise ValueError("`--score_lower_bound_z` must be finite and non-negative")
    if (
        not np.isfinite(float(args["broad_score_lower_bound_z"]))
        or float(args["broad_score_lower_bound_z"]) < 0.0
    ):
        raise ValueError(
            "`--broad_score_lower_bound_z` must be finite and non-negative"
        )
    if args["peak_mode"] in {"broad", "both"} and float(
        args["broad_score_lower_bound_z"]
    ) > float(args["score_lower_bound_z"]):
        raise ValueError(
            "`--broad_score_lower_bound_z` cannot exceed `--score_lower_bound_z`"
        )
    if args["min_peak_score"] is not None and (
        not np.isfinite(float(args["min_peak_score"]))
        or float(args["min_peak_score"]) < 0.0
    ):
        raise ValueError("`--min_peak_score` must be finite and non-negative")
    if args["broad_max_gap_bp"] is not None and int(args["broad_max_gap_bp"]) <= 0:
        raise ValueError("`--broad_max_gap_bp` must be positive")
    if int(args["broad_min_peak_bp"]) <= 0:
        raise ValueError("`--broad_min_peak_bp` must be positive")
    if int(args["window_bp"]) <= 0:
        raise ValueError("`--window_bp` must be positive")
    if int(args["window_count"]) < 20:
        raise ValueError("`--window_count` must be at least 20")
    if not 0.5 < float(args["working_quantile"]) < 1.0:
        raise ValueError("`--working_quantile` must be strictly between 0.5 and one")
    if int(args["bootstrap_draws"]) < 20:
        raise ValueError("`--bootstrap_draws` must be at least 20")
    if args["insufficient_data_policy"] not in {"error", "priorOnly"}:
        raise ValueError(
            "`--insufficient_data_policy` must be either `error` or `priorOnly`"
        )
    prior_radius_bp = args.get("prior_radius_bp")
    if prior_radius_bp is not None and (
        not np.isfinite(float(prior_radius_bp)) or float(prior_radius_bp) <= 0.0
    ):
        raise ValueError("`--prior_radius_bp` must be positive and finite")
    if args["insufficient_data_policy"] == "priorOnly" and prior_radius_bp is None:
        raise ValueError(
            "`--prior_radius_bp` is required with `insufficient_data_policy=priorOnly`"
        )
    if args["peak_mode"] is not None and args["input_track_type"] != "bam":
        raise ValueError("`--peak_mode` sidecars require BAM inputs")
    if (
        args["peak_mode"] in {"broad", "both"}
        and args.get("score_min_effect") is not None
    ):
        raise ValueError("Broad peak mode requires z-score optimization")

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


def _prepare_inputs(args: dict) -> list:
    signal_inputs = []
    expected_input_type = args.get("input_track_type", "")
    for file_ in args["input_files"]:
        if not os.path.exists(file_):
            raise FileNotFoundError(f"File not found: {file_}")
        input_type = get_track_type(file_)
        if expected_input_type != "" and input_type != expected_input_type:
            raise ValueError("All input files must share the same type.")
        signal_inputs.append(file_)

    if args["ignore_for_norm"] is None or len(args["ignore_for_norm"]) == 0:
        args["ignore_for_norm"] = ["chrX", "chrY", "chrM"]

    return signal_inputs


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


def _make_null_blocks(num_loci: int, num_blocks: int) -> list[tuple[int, int]]:
    num_loci_ = int(num_loci)
    num_blocks_ = int(num_blocks)
    if num_loci_ <= 0:
        raise ValueError("`num_loci` must be positive")
    if num_blocks_ <= 0:
        raise ValueError("`num_null_blocks` must be positive")
    block_count = int(min(num_blocks_, num_loci_))
    edges = np.linspace(0, num_loci_, block_count + 1, dtype=np.int64)
    return [
        (int(edges[idx]), int(edges[idx + 1]))
        for idx in range(block_count)
        if int(edges[idx + 1]) > int(edges[idx])
    ]


def _budget_block_record(
    block_id: int,
    start_idx: int,
    stop_idx: int,
    chrom_intervals: np.ndarray,
    interval_bp: int,
    budget_fraction_hat: float,
    budget_rate_meta: dict,
    working_span: int,
    budget_null_draws: int,
) -> dict:
    block_length = int(stop_idx) - int(start_idx)
    if block_length <= 0:
        raise ValueError("Budget blocks must be non-empty")
    if not np.isfinite(float(budget_fraction_hat)):
        raise ValueError("Budget estimate is not finite")
    budget_total_count_hat = float(
        np.clip(
            budget_rate_meta.get("effective_total_count", block_length),
            1.0,
            block_length,
        )
    )
    budget_count_hat = float(
        np.clip(
            budget_rate_meta.get(
                "budget_count_hat",
                float(budget_fraction_hat) * budget_total_count_hat,
            ),
            0.0,
            budget_total_count_hat,
        )
    )
    return {
        "block_id": int(block_id),
        "start_idx": int(start_idx),
        "stop_idx": int(stop_idx),
        "start_bp": int(chrom_intervals[int(start_idx)]),
        "stop_bp": int(chrom_intervals[int(stop_idx) - 1] + int(interval_bp)),
        "num_loci": int(block_length),
        "budget_fraction_hat": float(budget_fraction_hat),
        "budget_count_hat": float(budget_count_hat),
        "total_count": float(budget_total_count_hat),
        "effective_total_count": float(budget_total_count_hat),
        "working_span_intervals": int(
            budget_rate_meta.get("correlation_length_intervals", working_span)
        ),
        "dwb_bandwidth": int(budget_rate_meta.get("dwb_bandwidth", working_span)),
        "num_null_draws": int(
            budget_rate_meta.get("num_null_draws", budget_null_draws)
        ),
        "requested_num_null_draws": int(budget_null_draws),
    }


def _budget_value(chrom_budget) -> float:
    if isinstance(chrom_budget, dict):
        return float(chrom_budget["budget"])
    budget_arr = np.asarray(chrom_budget, dtype=np.float64)
    if budget_arr.ndim == 0:
        return float(budget_arr)
    if budget_arr.ndim == 2 and budget_arr.shape[1] == 3:
        lengths = budget_arr[:, 1] - budget_arr[:, 0]
        return float(np.average(budget_arr[:, 2], weights=lengths))
    raise ValueError("Chromosome budget has an unsupported shape")


def _solver_budget_blocks(chrom_budget) -> np.ndarray | None:
    if isinstance(chrom_budget, dict):
        if chrom_budget.get("mode") != "block":
            return None
        bounds = np.asarray(chrom_budget["block_bounds"], dtype=np.float64)
        budgets = np.asarray(chrom_budget["block_budgets"], dtype=np.float64)
    else:
        budget_arr = np.asarray(chrom_budget, dtype=np.float64)
        if budget_arr.ndim != 2 or budget_arr.shape[1] != 3:
            return None
        return np.ascontiguousarray(budget_arr, dtype=np.float64)
    if bounds.ndim != 2 or bounds.shape[1] != 2:
        raise ValueError("Block budget bounds must have start and stop columns")
    if budgets.ndim != 1 or budgets.shape[0] != bounds.shape[0]:
        raise ValueError("Block budget values must match block bounds")
    return np.column_stack((bounds, budgets))


def _solve_cached_chromosome(chrom_: str) -> tuple[str, float, dict, np.ndarray, str]:
    state = _CHROM_SOLVE_PROCESS_STATE
    if state is None:
        raise RuntimeError("Chromosome solve state is not initialized")
    chrom_data = state["chrom_cache"][chrom_]
    chrom_gamma = chrom_data["gamma"]
    if not np.all(np.isfinite(chrom_data["scores"])):
        raise ValueError(f"{chrom_} scores contain non-finite values")
    try:
        chrom_gamma = float(chrom_gamma)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{chrom_} gamma could not be read as a finite number"
        ) from exc
    if not np.isfinite(chrom_gamma) or chrom_gamma < 0.0:
        raise ValueError(f"{chrom_} gamma must be finite and non-negative")
    chrom_budget = None
    chrom_budget_blocks = None
    if state["selection_penalty"] is None:
        chrom_budget_spec = state["chrom_budgets"][chrom_]
        chrom_budget_blocks = _solver_budget_blocks(chrom_budget_spec)
        if chrom_budget_blocks is None:
            chrom_budget = _budget_value(chrom_budget_spec)
            if not np.isfinite(chrom_budget) or chrom_budget < 0.0:
                raise ValueError(f"{chrom_} budget must be finite and non-negative")
        else:
            if np.any(~np.isfinite(chrom_budget_blocks[:, 2])) or np.any(
                chrom_budget_blocks[:, 2] < 0.0
            ):
                raise ValueError(
                    f"{chrom_} block budgets must be finite and non-negative"
                )
    chrom_sol, chrom_obj, chrom_meta = solve_chrom_exact(
        chrom_data["scores"],
        budget=chrom_budget,
        budget_blocks=chrom_budget_blocks,
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
    return chrom_, float(chrom_obj), chrom_meta, chrom_sol, chrom_outfile


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
    with tempfile.TemporaryDirectory(prefix="rocco_dependence_") as spool_dir:
        for chrom_index, chrom_ in enumerate(chroms_to_process):
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

            chrom_intervals = np.asarray(chrom_intervals, dtype=np.int64)
            chrom_matrix = np.asarray(chrom_matrix)
            logger.info("Chromosome %s matrix: %s", chrom_, chrom_matrix.shape)
            if chrom_matrix.ndim != 2 or chrom_matrix.shape[0] != len(signal_inputs):
                raise ValueError(f"{chrom_} matrix rows must match input tracks")
            if chrom_matrix.shape[1] != chrom_intervals.size:
                raise ValueError(f"{chrom_} matrix and coordinates do not align")
            interval_diffs = np.diff(chrom_intervals)
            if interval_diffs.size > 0:
                if np.any(interval_diffs <= 0) or np.unique(interval_diffs).size != 1:
                    raise ValueError(
                        f"{chrom_} coordinates must use one fixed bin width"
                    )
                interval_bp = int(interval_diffs[0])
            else:
                interval_bp = int(args["step"])
            if interval_bp <= 0:
                raise ValueError(f"{chrom_} bin width must be positive")

            preparation_meta = {}
            if args["input_track_type"] == "bigwig":
                if np.any(np.isinf(chrom_matrix)):
                    raise ValueError(f"{chrom_} bigWig matrix contains infinite values")
                if chrom_matrix.shape[0] > 1:
                    logger.warning(
                        "Multiple bigWig tracks detected for %s. Region scores use their column median while dependence estimation retains every supplied row.",
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
                    raise ValueError(
                        f"{chrom_} has loci with no supplied bigWig score in any track"
                    )
                score_details = {
                    "mean": chrom_scores.astype(np.float64, copy=False),
                }
                dependence_matrix = chrom_matrix
                preparation_meta = {
                    "signalPreparation": "asProvided",
                    "scoringMethod": (
                        "direct" if chrom_matrix.shape[0] == 1 else "directMedian"
                    ),
                }
            else:
                if not np.all(np.isfinite(chrom_matrix)):
                    raise ValueError(f"{chrom_} BAM matrix contains non-finite values")
                chrom_scores, score_details = score_loci_wls(
                    chrom_matrix,
                    interval_bp=interval_bp,
                    lower_bound_z=args["score_lower_bound_z"],
                    prior_df=args["score_prior_df"],
                    min_effect=args.get("score_min_effect"),
                    precision_floor_ratio=args["score_precision_floor_ratio"],
                    return_details=True,
                )
                if not np.all(np.isfinite(chrom_scores)):
                    raise ValueError(f"{chrom_} scores contain non-finite values")
                dependence_matrix = np.asarray(score_details.pop("centered_matrix"))
                if not np.all(np.isfinite(dependence_matrix)):
                    raise ValueError(
                        f"{chrom_} centered matrix contains non-finite values"
                    )
                preparation_meta = {
                    "signalPreparation": "log2p1SavgolOrder0",
                    "scoringMethod": "WLS",
                    "centeringMethod": score_details["centeringMethod"],
                    "centeringWindowBP": int(score_details["centeringWindowBP"]),
                    "centeringWindowBins": int(score_details["centeringWindowBins"]),
                }

            matrix_path = os.path.join(spool_dir, f"matrix_{chrom_index}.npy")
            np.save(matrix_path, dependence_matrix, allow_pickle=False)
            chrom_effect_mean = np.asarray(
                score_details.get("mean", chrom_scores),
                dtype=np.float64,
            )
            chrom_cache[chrom_] = {
                "intervals": chrom_intervals,
                "scores": chrom_scores,
                "effect_mean": chrom_effect_mean,
                "z_scores": score_details.get("z_scores"),
                "gamma": float(args["gamma"]),
                "interval_bp": int(interval_bp),
                "num_loci": int(chrom_scores.shape[0]),
                "_dependence_matrix_path": matrix_path,
                **preparation_meta,
            }
            score_label = (
                "direct input scores"
                if args["input_track_type"] == "bigwig"
                else "WLS scores"
            )
            logger.info(
                "%s %s:%s", chrom_, score_label, cscores_quantiles(chrom_scores)
            )

        if len(chrom_cache) == 0:
            raise ValueError("No chromosome matrices were available for inference")
        interval_bps = {
            int(chrom_data["interval_bp"]) for chrom_data in chrom_cache.values()
        }
        if len(interval_bps) != 1:
            raise ValueError(
                "Genome-level dependence estimation requires one bin width across chromosomes"
            )
        dependence_step_bp = int(next(iter(interval_bps)))
        chromosome_matrices = _SpoolMatrixMapping(
            {
                chrom_: chrom_data["_dependence_matrix_path"]
                for chrom_, chrom_data in chrom_cache.items()
            }
        )
        chromosome_coordinates = {
            chrom_: chrom_data["intervals"]
            for chrom_, chrom_data in chrom_cache.items()
        }
        (
            correlation_radius,
            correlation_radius_lower,
            correlation_radius_upper,
            dependence_diagnostics,
        ) = choose_dependence_span(
            chromosome_matrices,
            chromosome_coordinates,
            dependence_step_bp,
            window_bp=int(args.get("window_bp", 50000)),
            window_count=int(args.get("window_count", 256)),
            working_quantile=float(args.get("working_quantile", 0.95)),
            bootstrap_draws=int(args.get("bootstrap_draws", 500)),
            insufficient_data_policy=args.get("insufficient_data_policy", "error"),
            prior_radius_bp=args.get("prior_radius_bp"),
        )
        del chromosome_matrices
        working_span = int(dependence_diagnostics["workingSpanIntervals"])
        if working_span <= 0:
            raise ValueError("Genome-level dependence working span must be positive")

        num_null_blocks = int(args.get("num_null_blocks", 4))
        if num_null_blocks <= 0:
            raise ValueError("`--num_null_blocks` must be positive")
        for chrom_, chrom_data in chrom_cache.items():
            chrom_scores = chrom_data["scores"]
            interval_bp = int(chrom_data["interval_bp"])
            budget_blocks = []
            block_slices = _make_null_blocks(chrom_scores.shape[0], num_null_blocks)
            chrom_interval_arr = np.asarray(chrom_data["intervals"], dtype=np.int64)
            dependence_matrix = None
            if args["input_track_type"] == "bam":
                dependence_matrix = np.load(
                    chrom_data["_dependence_matrix_path"],
                    mmap_mode="r",
                    allow_pickle=False,
                )
            for block_id, (start_idx, stop_idx) in enumerate(block_slices):
                block_scores = chrom_scores[start_idx:stop_idx]
                progress_label = (
                    f"Budget null {chrom_} block {block_id + 1}/{len(block_slices)}"
                )
                if dependence_matrix is None:
                    block_fraction_hat, block_rate_meta = (
                        estimate_budget_nonnull_fraction_from_score_track(
                            block_scores,
                            correlation_length=working_span,
                            step_bp=interval_bp,
                            num_null_draws=args["budget_null_draws"],
                            random_seed=1009 * int(block_id),
                            progress_label=progress_label,
                            return_details=True,
                        )
                    )
                else:
                    block_fraction_hat, block_rate_meta = (
                        estimate_budget_nonnull_fraction_from_wild_bootstrap_null(
                            dependence_matrix[:, start_idx:stop_idx],
                            observed_scores=block_scores,
                            lower_bound_z=args["score_lower_bound_z"],
                            prior_df=args["score_prior_df"],
                            min_effect=args.get("score_min_effect"),
                            precision_floor_ratio=args["score_precision_floor_ratio"],
                            correlation_length=working_span,
                            step_bp=interval_bp,
                            num_null_draws=args["budget_null_draws"],
                            random_seed=1009 * int(block_id),
                            progress_label=progress_label,
                            num_processes=min(
                                int(args["budget_null_draws"]),
                                int(budget_null_processes),
                            ),
                            return_details=True,
                        )
                    )
                budget_blocks.append(
                    _budget_block_record(
                        block_id,
                        start_idx,
                        stop_idx,
                        chrom_interval_arr,
                        interval_bp,
                        block_fraction_hat,
                        block_rate_meta,
                        working_span,
                        args["budget_null_draws"],
                    )
                )
            if dependence_matrix is not None:
                del dependence_matrix
            budget_total_count_hat = float(
                sum(float(block["total_count"]) for block in budget_blocks)
            )
            budget_count_hat = float(
                sum(float(block["budget_count_hat"]) for block in budget_blocks)
            )
            budget_fraction_hat = float(
                np.clip(
                    budget_count_hat / max(budget_total_count_hat, 1.0),
                    0.0,
                    1.0,
                )
            )
            dwb_bandwidth = int(
                max(int(block["dwb_bandwidth"]) for block in budget_blocks)
            )
            total_null_draws = int(
                sum(int(block["num_null_draws"]) for block in budget_blocks)
            )
            lean_budget_meta = {
                "budget_fraction_hat": float(budget_fraction_hat),
                "budget_count_hat": float(budget_count_hat),
                "effective_total_count": float(budget_total_count_hat),
                "correlation_radius_intervals": int(correlation_radius),
                "correlation_radius_bp": float(dependence_diagnostics["estimateBP"]),
                "working_span_intervals": int(working_span),
                "working_span_bp": float(dependence_diagnostics["workingSpanBP"]),
                "dwb_bandwidth": int(dwb_bandwidth),
                "num_null_draws_per_block": int(args["budget_null_draws"]),
                "total_null_draws": int(total_null_draws),
                "num_null_blocks": int(len(budget_blocks)),
                "requested_num_null_blocks": int(num_null_blocks),
                "budget_blocks": budget_blocks,
            }
            logger.info("%s raw budget estimate: %s", chrom_, lean_budget_meta)
            chrom_data.update(
                {
                    "budget_count_hat": float(budget_count_hat),
                    "budget_fraction_hat": float(budget_fraction_hat),
                    "budget_rate_meta": lean_budget_meta,
                    "correlation_radius_intervals": int(correlation_radius),
                    "correlation_radius_lower_intervals": int(correlation_radius_lower),
                    "correlation_radius_upper_intervals": int(correlation_radius_upper),
                    "correlation_radius_bp": float(
                        dependence_diagnostics["estimateBP"]
                    ),
                    "working_span_intervals": int(working_span),
                    "working_span_bp": float(dependence_diagnostics["workingSpanBP"]),
                    "dependence_diagnostics": dependence_diagnostics,
                    "dwb_bandwidth": int(dwb_bandwidth),
                    "num_null_draws_per_block": int(args["budget_null_draws"]),
                    "total_null_draws": int(total_null_draws),
                    "total_count": float(budget_total_count_hat),
                }
            )
            chrom_data.pop("_dependence_matrix_path")

    if (
        args.get("peak_mode") in {"narrow", "both"}
        and args["input_track_type"] == "bam"
    ):
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
    num_null_blocks = int(args.get("num_null_blocks", 4))
    budget_unit_counts: dict[str, float] = {}
    budget_unit_totals: dict[str, float] = {}
    block_units: list[tuple[str, str, dict]] = []
    for chrom, chrom_data in chrom_cache.items():
        budget_blocks = chrom_data.get("budget_rate_meta", {}).get("budget_blocks", [])
        if len(budget_blocks) == 0:
            budget_unit_counts[chrom] = chrom_data["budget_count_hat"]
            budget_unit_totals[chrom] = chrom_data["total_count"]
            continue
        for block in budget_blocks:
            block_id = int(block["block_id"])
            unit_key = f"{chrom}:{block_id}"
            budget_unit_counts[unit_key] = float(block["budget_count_hat"])
            budget_unit_totals[unit_key] = float(block["total_count"])
            block_units.append((unit_key, chrom, block))

    unit_budgets, budget_meta = estimate_empirical_bayes_budgets(
        budget_unit_counts,
        budget_unit_totals,
        posterior_quantile=args["budget_posterior_quantile"],
    )
    manual_selection_penalty = args.get("selection_penalty") is not None
    if (
        not manual_selection_penalty
        and args["budget"] is not None
        and budget_meta["genome_wide_budget"] > 0
    ):
        rescale = float(args["budget"]) / budget_meta["genome_wide_budget"]
    else:
        rescale = 1.0
    unit_budgets = {
        unit_key: min(
            max(
                unit_budgets[unit_key] * rescale * float(args["scale_chrom_budgets"]),
                0.001,
            ),
            0.25,
        )
        for unit_key in unit_budgets
    }
    chrom_budgets = {}
    for chrom, chrom_data in chrom_cache.items():
        budget_blocks = chrom_data.get("budget_rate_meta", {}).get("budget_blocks", [])
        if len(budget_blocks) == 0:
            chrom_budgets[chrom] = float(unit_budgets[chrom])
            continue
        resolved_blocks = []
        for block in budget_blocks:
            block_id = int(block["block_id"])
            start_idx = int(block["start_idx"])
            stop_idx = int(block["stop_idx"])
            start_bp = int(block["start_bp"])
            stop_bp = int(block["stop_bp"])
            unit_key = f"{chrom}:{block_id}"
            resolved_budget = float(unit_budgets[unit_key])
            resolved_blocks.append(
                {
                    "block_id": int(block_id),
                    "start_idx": int(start_idx),
                    "stop_idx": int(stop_idx),
                    "start_bp": int(start_bp),
                    "stop_bp": int(stop_bp),
                    "num_loci": int(block["num_loci"]),
                    "raw_budget_fraction_hat": float(block["budget_fraction_hat"]),
                    "raw_budget_count_hat": float(block["budget_count_hat"]),
                    "total_count": float(block["total_count"]),
                    "budget": resolved_budget,
                }
            )
        block_bounds = np.asarray(
            [
                (int(block["start_idx"]), int(block["stop_idx"]))
                for block in resolved_blocks
            ],
            dtype=np.float64,
        )
        block_budgets = np.asarray(
            [float(block["budget"]) for block in resolved_blocks],
            dtype=np.float64,
        )
        block_lengths = block_bounds[:, 1] - block_bounds[:, 0]
        chrom_budgets[chrom] = {
            "mode": "block",
            "budget": float(np.average(block_budgets, weights=block_lengths)),
            "block_bounds": block_bounds,
            "block_budgets": block_budgets,
            "blocks": resolved_blocks,
        }
    budget_meta = dict(budget_meta)
    budget_meta["num_null_blocks"] = int(num_null_blocks)
    budget_meta["manual_selection_penalty"] = bool(manual_selection_penalty)
    budget_meta["budget_rescale"] = float(rescale)
    budget_meta["budget_unit_count"] = int(len(budget_unit_counts))
    budget_meta["budget_unit_scope"] = (
        "chromosome_blocks" if len(block_units) > 0 else "chromosomes"
    )
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
            round(_budget_value(chrom_budget), 6),
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

    for chrom_, chrom_obj, chrom_meta, chrom_solution, chrom_outfile in solve_results:
        chrom_cache[chrom_]["solution"] = np.asarray(chrom_solution, dtype=np.uint8)
        solve_meta = {
            "peak_mode": args.get("peak_mode"),
            "score_lower_bound_z": float(args["score_lower_bound_z"]),
            "broad_score_lower_bound_z": float(args["broad_score_lower_bound_z"]),
            "min_peak_score": float(args["min_peak_score"]),
            "gamma": float(chrom_cache[chrom_]["gamma"]),
            "budget_mode": str(chrom_meta["budget_mode"]),
            "soft_budget_penalty": float(chrom_meta["soft_budget_penalty"]),
            "budget": float(_budget_value(chrom_budgets[chrom_])),
            "selected_count": int(chrom_meta["selected_count"]),
            "correlation_radius_intervals": int(
                chrom_cache[chrom_]["correlation_radius_intervals"]
            ),
            "correlation_radius_bp": float(
                chrom_cache[chrom_]["correlation_radius_bp"]
            ),
            "working_span_intervals": int(
                chrom_cache[chrom_]["working_span_intervals"]
            ),
            "working_span_bp": float(chrom_cache[chrom_]["working_span_bp"]),
            "dwb_bandwidth": int(chrom_cache[chrom_]["dwb_bandwidth"]),
            "num_null_draws_per_block": int(
                chrom_cache[chrom_]["num_null_draws_per_block"]
            ),
            "total_null_draws": int(chrom_cache[chrom_]["total_null_draws"]),
        }
        for preparation_key in (
            "signalPreparation",
            "scoringMethod",
            "centeringMethod",
            "centeringWindowBP",
            "centeringWindowBins",
        ):
            if preparation_key in chrom_cache[chrom_]:
                solve_meta[preparation_key] = chrom_cache[chrom_][preparation_key]
        if "budget_block_count" in chrom_meta:
            solve_meta["budget_blocks"] = [
                {
                    "start_idx": int(start_idx),
                    "stop_idx": int(stop_idx),
                    "budget": float(block_budget),
                    "selection_penalty": float(block_penalty),
                }
                for start_idx, stop_idx, block_budget, block_penalty in zip(
                    chrom_meta["budget_block_starts"],
                    chrom_meta["budget_block_ends"],
                    chrom_meta["budget_block_fractions"],
                    chrom_meta["budget_block_penalties"],
                )
            ]
        chrom_cache[chrom_]["solve_meta"] = solve_meta
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


def _peak_sidecar_root(final_output: str) -> str:
    output_root, output_ext = os.path.splitext(final_output)
    if output_ext == ".bed":
        return output_root
    return final_output


def _filter_scored_peak_file_by_signal(
    peak_file: str,
    min_peak_score: float,
) -> int:
    kept_lines = []
    with open(peak_file, "r", encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, start=1):
            line_ = line.rstrip("\n")
            if line_.strip() == "":
                continue
            fields = line_.split("\t")
            if len(fields) < 7:
                raise ValueError(
                    f"Scored peak row {line_num} in {peak_file} has fewer than 7 columns."
                )
            if float(fields[6]) >= float(min_peak_score):
                kept_lines.append(line_)
    with open(peak_file, "w", encoding="utf-8") as handle:
        for line_ in kept_lines:
            handle.write(f"{line_}\n")
    return int(len(kept_lines))


def _solution_runs(
    solution: np.ndarray,
    intervals: np.ndarray,
) -> list[tuple[int, int, int, int]]:
    intervals_ = np.asarray(intervals, dtype=np.int64)
    solution_ = np.asarray(solution)
    usable = int(min(solution_.shape[0], intervals_.shape[0]))
    if usable <= 0:
        return []
    if usable == 1:
        step_bp = 1
    else:
        interval_diffs = np.diff(intervals_[:usable])
        positive_diffs = interval_diffs[interval_diffs > 0]
        if positive_diffs.size == 0:
            raise ValueError("Intervals must contain increasing bin starts")
        step_bp = int(np.median(positive_diffs))
    selected = np.asarray(solution_[:usable] > 0.5, dtype=np.int8)
    if not np.any(selected):
        return []
    padded = np.zeros(usable + 2, dtype=np.int8)
    padded[1 : usable + 1] = selected
    edges = np.diff(padded)
    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1)
    return [
        (
            int(start_idx),
            int(stop_idx),
            int(intervals_[start_idx]),
            int(intervals_[stop_idx - 1] + step_bp),
        )
        for start_idx, stop_idx in zip(starts, stops)
        if int(intervals_[stop_idx - 1] + step_bp) > int(intervals_[start_idx])
    ]


def _build_broad_parent_records(
    chrom_cache: dict,
    chrom_budgets: dict,
    args: dict,
) -> tuple[
    list[tuple[str, int, int]], dict[tuple[str, int, int], list[tuple[int, int]]]
]:
    records: list[tuple[str, int, int]] = []
    block_map: dict[tuple[str, int, int], list[tuple[int, int]]] = {}
    for chrom_, chrom_data in chrom_cache.items():
        strong_runs = _solution_runs(
            chrom_data["solution"],
            chrom_data["intervals"],
        )
        if len(strong_runs) == 0:
            continue
        z_scores = chrom_data.get("z_scores")
        if z_scores is None:
            raise ValueError("Broad peak mode requires stored z-scores")
        weak_scores = np.asarray(z_scores, dtype=np.float64) - float(
            args["broad_score_lower_bound_z"]
        )
        weak_solution, _, weak_meta = solve_chrom_exact(
            weak_scores,
            budget=None,
            gamma=float(chrom_data["gamma"]),
            selection_penalty=0.0,
            return_details=True,
        )
        chrom_data["broad_solution"] = weak_solution.astype(np.uint8, copy=False)
        chrom_data["broad_solve_meta"] = {
            "budget_mode": str(weak_meta["budget_mode"]),
            "soft_budget_penalty": float(weak_meta["soft_budget_penalty"]),
            "selected_count": int(weak_meta["selected_count"]),
            "budget": float(_budget_value(chrom_budgets[chrom_])),
        }
        weak_runs = _solution_runs(weak_solution, chrom_data["intervals"])
        if len(weak_runs) == 0:
            continue

        strong_starts = np.asarray([run[2] for run in strong_runs], dtype=np.int64)
        strong_ends = np.asarray([run[3] for run in strong_runs], dtype=np.int64)
        prefix_strong_ends = np.maximum.accumulate(strong_ends)
        candidate_runs = []
        for run in weak_runs:
            _, _, run_start, run_end = run
            strong_idx = int(np.searchsorted(strong_starts, run_end, side="left") - 1)
            if strong_idx >= 0 and int(prefix_strong_ends[strong_idx]) > int(run_start):
                candidate_runs.append(run)
        if len(candidate_runs) == 0:
            continue

        interval_bp = int(max(1, chrom_data["interval_bp"]))
        if args["broad_max_gap_bp"] is None:
            max_gap_bp = int(
                max(
                    interval_bp,
                    2 * int(chrom_data["working_span_intervals"]) * interval_bp,
                )
            )
        else:
            max_gap_bp = int(args["broad_max_gap_bp"])
        merged: list[list[int]] = []
        for _, _, run_start, run_end in candidate_runs:
            if len(merged) == 0 or int(run_start) - int(merged[-1][1]) > max_gap_bp:
                merged.append([int(run_start), int(run_end)])
            else:
                merged[-1][1] = max(int(merged[-1][1]), int(run_end))

        for parent_start, parent_end in merged:
            if int(parent_end) - int(parent_start) < int(args["broad_min_peak_bp"]):
                continue
            blocks = []
            for _, _, strong_start, strong_end in strong_runs:
                block_start = int(max(strong_start, parent_start))
                block_end = int(min(strong_end, parent_end))
                if block_end > block_start:
                    blocks.append((block_start, block_end))
            if len(blocks) == 0:
                continue
            key = (chrom_, int(parent_start), int(parent_end))
            records.append(key)
            block_map[key] = blocks
    return records, block_map


def _write_gapped_peak_file(
    scored_parent_file: str,
    output_file: str,
    block_map: dict[tuple[str, int, int], list[tuple[int, int]]],
    min_peak_score: float,
) -> int:
    written = 0
    with open(scored_parent_file, "r", encoding="utf-8") as src, open(
        output_file,
        "w",
        encoding="utf-8",
    ) as dst:
        for line_num, line in enumerate(src, start=1):
            line_ = line.strip()
            if line_ == "":
                continue
            fields = line_.split("\t")
            if len(fields) < 9:
                raise ValueError(
                    f"Scored broad row {line_num} in {scored_parent_file} has fewer than 9 columns."
                )
            chrom = str(fields[0])
            start = int(fields[1])
            end = int(fields[2])
            signal_value = float(fields[6])
            if signal_value < float(min_peak_score):
                continue
            blocks = block_map[(chrom, start, end)]
            block_sizes = ",".join(
                str(int(block_end) - int(block_start))
                for block_start, block_end in blocks
            )
            block_starts = ",".join(
                str(int(block_start) - start) for block_start, _ in blocks
            )
            thick_start = int(blocks[0][0])
            thick_end = int(blocks[-1][1])
            dst.write(
                "\t".join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        str(fields[3]),
                        str(fields[4]),
                        str(fields[5]),
                        str(thick_start),
                        str(thick_end),
                        "0",
                        str(len(blocks)),
                        block_sizes,
                        block_starts,
                        str(fields[6]),
                        str(fields[7]),
                        str(fields[8]),
                    ]
                )
                + "\n"
            )
            written += 1
    return int(written)


def _generate_peak_mode_outputs(
    args: dict,
    final_output: str,
    chrom_cache: dict,
    chrom_budgets: dict,
):
    peak_mode = args.get("peak_mode")
    if peak_mode is None:
        return
    if args.get("input_track_type") != "bam":
        raise ValueError("Peak sidecars require BAM inputs")

    peak_threads = int(args["threads"])
    if peak_threads <= 0:
        peak_threads = None

    sidecar_root = _peak_sidecar_root(final_output)
    if peak_mode in {"narrow", "both"}:
        summit_offsets_file = None
        try:
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
            narrowpeak_filepath = f"{sidecar_root}.narrowPeak"
            posthoc_scores.score_peaks(
                args["input_files"],
                args["chrom_sizes_file"],
                final_output,
                count_matrix_file=f"{sidecar_root}.counts.tsv",
                output_file=narrowpeak_filepath,
                ecdf_nsamples=args["ecdf_samples"],
                seed=args["ecdf_seed"],
                proc=args["ecdf_proc"],
                threads=peak_threads,
                summit_offsets_file=summit_offsets_file,
            )
            kept_count = _filter_scored_peak_file_by_signal(
                narrowpeak_filepath,
                float(args["min_peak_score"]),
            )
            logger.info(
                "Final narrowPeak output: %s (%s peaks)",
                narrowpeak_filepath,
                kept_count,
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

    if peak_mode in {"broad", "both"}:
        broad_records, block_map = _build_broad_parent_records(
            chrom_cache,
            chrom_budgets,
            args,
        )
        gappedpeak_filepath = f"{sidecar_root}.gappedPeak"
        if len(broad_records) == 0:
            open(gappedpeak_filepath, "w", encoding="utf-8").close()
            logger.info("Final gappedPeak output: %s (0 peaks)", gappedpeak_filepath)
            return

        parent_bed_file = None
        scored_parent_file = None
        try:
            fd, parent_bed_file = tempfile.mkstemp(
                prefix="rocco_broad_parent_",
                suffix=".bed",
            )
            os.close(fd)
            fd, scored_parent_file = tempfile.mkstemp(
                prefix="rocco_broad_scored_",
                suffix=".narrowPeak",
            )
            os.close(fd)
            _write_bed_records(broad_records, parent_bed_file, name_features=True)
            posthoc_scores.score_peaks(
                args["input_files"],
                args["chrom_sizes_file"],
                parent_bed_file,
                count_matrix_file=f"{sidecar_root}.broad.counts.tsv",
                output_file=scored_parent_file,
                ecdf_nsamples=args["ecdf_samples"],
                seed=args["ecdf_seed"],
                proc=args["ecdf_proc"],
                threads=peak_threads,
            )
            kept_count = _write_gapped_peak_file(
                scored_parent_file,
                gappedpeak_filepath,
                block_map,
                float(args["min_peak_score"]),
            )
            logger.info(
                "Final gappedPeak output: %s (%s peaks)",
                gappedpeak_filepath,
                kept_count,
            )
        finally:
            for tmp_file in (parent_bed_file, scored_parent_file):
                if tmp_file is None:
                    continue
                try:
                    os.remove(tmp_file)
                except Exception as exc:
                    logger.info(
                        "Could not remove peak sidecar temp. file %s\n%s",
                        tmp_file,
                        exc,
                    )


def main():
    parser = _build_parser()
    args = _prepare_args(parser)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(module)s.%(funcName)s -  %(levelname)s - %(message)s",
    )
    run_id = str(int(uuid.uuid4().hex[:5], base=16))
    logger.info("\nID: %s", run_id)

    signal_inputs = _prepare_inputs(args)
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
        _generate_peak_mode_outputs(args, final_output, chrom_cache, chrom_budgets)
    finally:
        _cleanup_narrowpeak_tempfiles(chrom_cache)


if __name__ == "__main__":
    main()
