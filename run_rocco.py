#!/usr/bin/env python
import os
import argparse
import subprocess
import ROCCO_chrom

def subcommand_chrom(args):
    # run ROCCO_chrom.py
    ROCCO_chrom.main(args)
    return None

def subcommand_gw(args):
    # run ROCCO.py
    return None

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
    parser_subcommand_gw = subparsers.add_parser("subcommand1", help="Description for subcommand1")

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

    if args['command'] == "gw":
        # Call the relevant function for subcommand1
        pass
    elif args['command'] == "chrom":
        subcommand_chrom(args)
        pass
    else:
        parser.print_help()

if __name__ == "__main__":
    main()



