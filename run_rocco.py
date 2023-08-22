#!/usr/bin/env python
import os
import argparse
import ROCCO_chrom
import ROCCO
import prep_bams
import est_budgets

def subcommand_chrom(args):
    ROCCO_chrom.main(args)
    return True

def subcommand_gwide(args):
    ROCCO.main(args)
    return True

def subcommand_prep(args):
    prep_bams.main(args)
    return True

def subcommand_budgets(args):
    est_budgets.main(args)
    return True


def main():
    parser = argparse.ArgumentParser(description="Run ROCCO at the command line")
    subparsers = parser.add_subparsers(dest="command")

    # 'gwide' subcommand parameters
    parser_subcommand_gwide = subparsers.add_parser("gwide", help="Run ROCCO genome-wide")
    parser_subcommand_gwide.add_argument('-p', '--param_file',
                        default='params.csv',
                        help="CSV param file w/ row for each chromosome")
    parser_subcommand_gwide.add_argument('-b', '--budget', default=.035)
    parser_subcommand_gwide.add_argument('-g', '--gamma', default=1.0)
    parser_subcommand_gwide.add_argument('-t', '--tau', default=0.0)
    parser_subcommand_gwide.add_argument('--c1', default=1.0)
    parser_subcommand_gwide.add_argument('--c2', default=1.0)
    parser_subcommand_gwide.add_argument('--c3', default=1.0)
    parser_subcommand_gwide.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_gwide.add_argument('--solver', default="ECOS",
                        help="Optimization software used to solve the \
                        LP underlying ROCCO. `ECOS` is used by default \
                        and is a viable open-source option. `MOSEK`\
                        offers significantly greater speed and is free\
                        for academic use. Free trial commerical licenses\
                        are also available. See\
                        https://www.mosek.com/products/academic-licenses/")
    parser_subcommand_gwide.add_argument('--bed_format', type=int, default=3,
                        help="Specifies BED3 or BED6 format.\
                        Default is BED6. Generate BED3 output with \
                        --bed_format 3")
    parser_subcommand_gwide.add_argument('--identifiers', default=None,
                        help="an optional filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a substring of the `.wig` sample")
    parser_subcommand_gwide.add_argument('--outdir', default='.',
                        help="directory in which to store ROCCO's output\
                          files.")
    parser_subcommand_gwide.add_argument('--combine', default=None, help="if not None, combine\
                        output bed files and store in the file specified\
                        with this parameter. ex: `--combine bedname.bed` con-\
                        catenates the chromosome-specific bedfiles into `bedname.bed`.")
    parser_subcommand_gwide.add_argument('--multi', default=False, action='store_true', help='run ROCCO_chrom.py jobs\
        simultaneously to improve speed. increases peak memory expense.')
    parser_subcommand_gwide.add_argument('--verbose', default=False, action="store_true")


    # 'chrom' subcommand parameters
    parser_subcommand_chrom = subparsers.add_parser("chrom", help="Run ROCCO on a single chromosome")
    parser_subcommand_chrom.add_argument('--start', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--end', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--locus_size', type=int, default=-1,
                        help="this must match the constant step-size \
                        observed in the input wiggle files")
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

    # 'prep' subcommand parameters
    parser_subcommand_prep = subparsers.add_parser("prep", help="preprocess BAM files")
    parser_subcommand_prep.add_argument('-i','--bamdir', default='.', type=str, help='path to directory containing BAM files')
    parser_subcommand_prep.add_argument('-o', '--outdir', default= '.', help='output directory')
    parser_subcommand_prep.add_argument('-s', '--sizes',
                        default='hg38',
                        help="A path to a chromosome sizes file. OR\
                        an assembly name")
    parser_subcommand_prep.add_argument('-L', '--interval_length',
                        default=50,
                        help="wiggle track fixed interval length")
    parser_subcommand_prep.add_argument('-c', '--cores',
                        type=int,
                        default=1,
                        help="`bamSitesToWig.py`'s cores  parameter.\
                        Altering this parameter to use >1 core\
                        may cause issues on Mac OS")
    parser_subcommand_prep.add_argument('--multi',
                         default=True,
                         help='`True` to run `bamSitesToWig` jobs\
                         simultaneously.')
    parser_subcommand_prep.add_argument('--bstw_path',
                            default='pepatac/bamSitesToWig.py',
                            help="path to bamSitesToWig.py script")

    # 'budgets' subcommand parameters
    parser_subcommand_budgets = subparsers.add_parser("budgets",
                            help="estimate budget parameter using read densities")
    parser_subcommand_budgets.add_argument('-d', '--bamdir', type=str, help='Path to the directory containing .bam and .bai files')
    parser_subcommand_budgets.add_argument('-s', '--sizes', type=str, help='chromosome sizes file')
    parser_subcommand_budgets.add_argument('-a', type=float, default=0.0, help='Minimum allowed budget. This bound is ignored if\
        `--desired_avg` is non-negative.')
    parser_subcommand_budgets.add_argument('-b', type=float, default=0.05, help='Maximum allowed budget. This bound is ignored if\
        `--desired_avg` is non-negative.')
    parser_subcommand_budgets.add_argument('--desired_avg', type=float, default=-1.0, help='Scaled read densities (i.e., budgets) will\
        average to this value if non-negative. Defaults to -1.')

    args = vars(parser.parse_args())

    if args['command'] == "gwide":
        subcommand_gwide(args)
    elif args['command'] == "chrom":
        subcommand_chrom(args)
    elif args['command'] == 'prep':
        subcommand_prep(args)
    elif args['command'] == 'budgets':
        subcommand_budgets(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()



