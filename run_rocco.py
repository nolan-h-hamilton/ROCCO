#!/usr/bin/env python
import os
import argparse
import subprocess
import ROCCO_chrom
import ROCCO

def subcommand_chrom(args):
    ROCCO_chrom.main(args)
    return True

def subcommand_gw(args):
    ROCCO.main(args)
    return True

def subcommand_prep_bams(args):
    # run prep_bams.py
    return None

def subcommand_est_budgets(args):
    # run est_budgets.py
    return None


def main():
    parser = argparse.ArgumentParser(description="ROCCO CLI")
    subparsers = parser.add_subparsers(dest="command")

    # Add subcommands here
    parser_subcommand_gw = subparsers.add_parser("gwide", help="Run ROCCO genome-wide")
    parser_subcommand_gw.add_argument('-p', '--param_file',
                        default='params.csv',
                        help="CSV param file w/ row for each chromosome")
    # if some parameters are left `NULL` in `param_file`,
    # the values given by the following arguments will be
    # used as defaults.
    parser_subcommand_gw.add_argument('-b', '--budget', default=.035)
    parser_subcommand_gw.add_argument('-g', '--gamma', default=1.0)
    parser_subcommand_gw.add_argument('-t', '--tau', default=0.0)
    parser_subcommand_gw.add_argument('--c1', default=1.0)
    parser_subcommand_gw.add_argument('--c2', default=1.0)
    parser_subcommand_gw.add_argument('--c3', default=1.0)
    parser_subcommand_gw.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_gw.add_argument('--solver', default="ECOS",
                        help="Optimization software used to solve the \
                        LP underlying ROCCO. `ECOS` is used by default \
                        and is a viable open-source option. `MOSEK`\
                        offers significantly greater speed and is free\
                        for academic use. Free trial commerical licenses\
                        are also available. See\
                        https://www.mosek.com/products/academic-licenses/")
    parser_subcommand_gw.add_argument('--bed_format', type=int, default=3,
                        help="Specifies BED3 or BED6 format.\
                        Default is BED6. Generate BED3 output with \
                        --bed_format 3")
    parser_subcommand_gw.add_argument('--identifiers', default=None,
                        help="an optional filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a substring of the `.wig` sample")
    parser_subcommand_gw.add_argument('--outdir', default='.',
                        help="directory in which to store ROCCO's output\
                          files.")
    parser_subcommand_gw.add_argument('--combine', default=None, help="if not None, combine\
                        output bed files and store in the file specified\
                        with this parameter. ex: `--combine bedname.bed` con-\
                        catenates the chromosome-specific bedfiles into `bedname.bed`.")
    parser_subcommand_gw.add_argument('-j', '--jobs', type=int, default=1,
                        help="deprecated: included for backwards compatibility. use `--multi` instead.")
    parser_subcommand_gw.add_argument('-m', '--mem', default=None,
                        help="deprecated: included for backwards compatibility. use `--multi` instead.")
    parser_subcommand_gw.add_argument(
        '--parlog',
        default='ROCCO_parlog.txt',
        help='deprecated: included for backwards compatibility. use `--multi` instead.')
    parser_subcommand_gw.add_argument('--multi', default=False, action='store_true', help='run ROCCO_chrom.py jobs\
        simultaneously to improve speed. increases peak memory expense.')
    parser_subcommand_gw.add_argument('--verbose', default=False, action="store_true")

    parser_subcommand_chrom = subparsers.add_parser("chrom", help="Run ROCCO on a single chromosome")
    parser_subcommand_chrom.add_argument('--start', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--end', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--locus_size', type=int, default=-1,
                        help="this must match the constant step-size \
                        in the wiggle files used as input")
    parser_subcommand_chrom.add_argument('-b', '--budget', default=.035)
    parser_subcommand_chrom.add_argument('--chrom', help="e.g., --chrom chr1")
    parser_subcommand_chrom.add_argument('--wig_path', type=str, default=os.getcwd())
    parser_subcommand_chrom.add_argument('-g', '--gamma', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('-t', '--tau', type=float, default=0.0)
    parser_subcommand_chrom.add_argument('--c1', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--c2', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--c3', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--solver', default="ECOS", help='solver software\
        used to solve LP. Both "ECOS" and "PDLP" are free/open-source')
    parser_subcommand_chrom.add_argument('--bed_format', type=int, default=6,
                        help="`3` for BED3 format and `6` for BED6 format")
    parser_subcommand_chrom.add_argument('--identifiers', default=None,
                        help="(optional) a filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a substring of the `.wig` sample.")
    parser_subcommand_chrom.add_argument('--outdir', type=str, default='.')
    parser_subcommand_chrom.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_chrom.add_argument('--verbose', default=False, action="store_true")
    args = vars(parser.parse_args())

    if args['command'] == "gwide":
        subcommand_gw(args)
        pass
    elif args['command'] == "chrom":
        subcommand_chrom(args)
        pass
    else:
        parser.print_help()

if __name__ == "__main__":
    main()



