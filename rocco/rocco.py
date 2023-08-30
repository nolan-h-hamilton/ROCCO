"""
```
ROCCO: [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization.

PyPI : https://pypi.org/project/rocco/

GitHub: https://github.com/nolan-h-hamilton/ROCCO/

Demo: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb

Documentation: https://nolan-h-hamilton.github.io/ROCCO/

Paper: https://doi.org/10.1101/2023.05.24.542132


usage: rocco [-h] {gwide,chrom,prep,budgets,get_sizes} ...
  {gwide,chrom,prep,budgets,get_sizes}
    gwide               run rocco genome-wide/on multiple chromosomes (ROCCO_gwide.py)
    chrom               run ROCCO on a single chromosome (ROCCO_chrom.py)
    prep                Preprocess BAM files (prep_bams.py)
    budgets             Compute a budget (maximum fraction of basepairs that can be selected as 'open') for
                        each chromosome ranked by average read-density observed in the supplied BAM files (est_budgets.py)
    get_sizes           download sizes file for a genome in the ucsc genome registry
options:
  -h, --help            show this help message and exit

```

Subcommand Documentation:
    [`gwide`](https://nolan-h-hamilton.github.io/ROCCO/rocco/ROCCO_gwide.html)

    [`chrom`](https://nolan-h-hamilton.github.io/ROCCO/rocco/ROCCO_chrom.html)

    [`prep`](https://nolan-h-hamilton.github.io/ROCCO/rocco/prep_bams.html)

    [`budgets`](https://nolan-h-hamilton.github.io/ROCCO/rocco/est_budgets.html)
    

"""
#!/usr/bin/env python
import os
import argparse
from . import rocco_aux
from . import ROCCO_chrom
from . import ROCCO_gwide
from . import prep_bams
from . import est_budgets
from . import locus
from . import loci


def subcommand_chrom(args):
    ROCCO_chrom.main(args)
    return True

def subcommand_gwide(args):
    ROCCO_gwide.main(args)
    return True

def subcommand_prep(args):
    prep_bams.main(args)
    return True

def subcommand_budgets(args):
    est_budgets.main(args)
    return True


def main():
    parser = argparse.ArgumentParser(description="ROCCO: [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization.\n\nPyPI : https://pypi.org/project/rocco/ \nGitHub: https://github.com/nolan-h-hamilton/ROCCO/ \nDemo: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb \nPaper: https://doi.org/10.1101/2023.05.24.542132\n", add_help=True, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest="command")

    # 'gwide' subcommand parameters
    parser_subcommand_gwide = subparsers.add_parser("gwide", help='run rocco genome-wide/on multiple chromosomes (ROCCO_gwide.py)')
    parser_subcommand_gwide.add_argument('-p', '--param_file', required=True,
                                         help='A CSV file specifying the parameters for each chromosome. See: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo_params.csv')
    parser_subcommand_gwide.add_argument('-b', '--budget', type=float, default=.035)
    parser_subcommand_gwide.add_argument('-g', '--gamma', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('-t', '--tau', type=float, default=0.0)
    parser_subcommand_gwide.add_argument('--c1', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('--c2', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('--c3', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_gwide.add_argument('--solver', type=str, default="ECOS")
    parser_subcommand_gwide.add_argument('--bed_format', type=int, default=3)
    parser_subcommand_gwide.add_argument('--identifiers', default=None)
    parser_subcommand_gwide.add_argument('--outdir', default='.')
    parser_subcommand_gwide.add_argument('--combine', default=None)
    parser_subcommand_gwide.add_argument('--multi', default=1, type=int, help='number of simultaneous `rocco chrom` jobs to execute')
    parser_subcommand_gwide.add_argument('--verbose', default=False, action="store_true")

    # 'chrom' subcommand parameters
    parser_subcommand_chrom = subparsers.add_parser("chrom", help='run ROCCO on a single chromosome (ROCCO_chrom.py)')
    parser_subcommand_chrom.add_argument('--start', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--end', type=int, default=-1)
    parser_subcommand_chrom.add_argument('--locus_size', type=int, default=-1)
    parser_subcommand_chrom.add_argument('-b', '--budget', type=float, default=.035)
    parser_subcommand_chrom.add_argument('--chrom', type=str, help='')
    parser_subcommand_chrom.add_argument('--wig_path', type=str, default=os.getcwd())
    parser_subcommand_chrom.add_argument('-g', '--gamma', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('-t', '--tau', type=float, default=0.0)
    parser_subcommand_chrom.add_argument('--c1', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--c2', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--c3', type=float, default=1.0)
    parser_subcommand_chrom.add_argument('--solver', type=str,default="ECOS")
    parser_subcommand_chrom.add_argument('--bed_format', type=int, default=6)
    parser_subcommand_chrom.add_argument('--identifiers', default=None)
    parser_subcommand_chrom.add_argument('--outdir', type=str, default='.')
    parser_subcommand_chrom.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_chrom.add_argument('--verbose', default=False, action="store_true")

    # 'prep' subcommand parameters
    parser_subcommand_prep = subparsers.add_parser("prep", help='Preprocess BAM files (prep_bams.py)')
    parser_subcommand_prep.add_argument('-i', '--bamdir', default='.', type=str)
    parser_subcommand_prep.add_argument('-o', '--outdir', default='.')
    parser_subcommand_prep.add_argument('-s', '--sizes', default='hg38')
    parser_subcommand_prep.add_argument('-L', '--interval_length', default=50)
    parser_subcommand_prep.add_argument('-c', '--cores', type=int, default=1)
    parser_subcommand_prep.add_argument('--multi', default=True)
    parser_subcommand_prep.add_argument('--bstw_path', default=None, help='deprecated')

    # 'budgets' subcommand parameters
    parser_subcommand_budgets = subparsers.add_parser("budgets", help='Compute a budget for each chromosome ranked by average read-density observed in the supplied BAM files (est_budgets.py)')
    parser_subcommand_budgets.add_argument('-i', '--bamdir', type=str)
    parser_subcommand_budgets.add_argument('-s', '--sizes', type=str)
    parser_subcommand_budgets.add_argument('-a', type=float, default=0.0)
    parser_subcommand_budgets.add_argument('-b', type=float, default=0.05)
    parser_subcommand_budgets.add_argument('--desired_avg', type=float, default=-1.0)
    parser_subcommand_budgets.add_argument('--index', default=False, action='store_true')

    # 'get_sizes' subcommand parameters
    parser_subcommand_get_sizes = subparsers.add_parser("get_sizes", help="download sizes file for a genome in the ucsc genome registry")
    parser_subcommand_get_sizes.add_argument('-g', '--genome', type=str, default='hg38')

    args = vars(parser.parse_args())
    if args['command'] == "gwide":
        subcommand_gwide(args)
    elif args['command'] == "chrom":
        subcommand_chrom(args)
    elif args['command'] == 'prep':
        subcommand_prep(args)
    elif args['command'] == 'budgets':
        subcommand_budgets(args)
    elif args['command'] == 'get_sizes':
        rocco_aux.get_size_file(args['genome'])
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
