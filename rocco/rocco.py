r"""
# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization.
Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin across many samples.**

#### GitHub (Homepage): https://github.com/nolan-h-hamilton/ROCCO/

#### Demo: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb

#### Paper: https://doi.org/10.1093/bioinformatics/btad725


Subcommand Documentation:
    [`gwide`](https://nolan-h-hamilton.github.io/ROCCO/rocco/gwide.html)

    [`chrom`](https://nolan-h-hamilton.github.io/ROCCO/rocco/chrom.html)

    [`prep`](https://nolan-h-hamilton.github.io/ROCCO/rocco/prep.html)

    [`budgets`](https://nolan-h-hamilton.github.io/ROCCO/rocco/budgets.html)

General Notation and Terminology:
    - $\mathscr{L}$: genomic region, e.g., a chromosome
    - $L$: size of loci in bp, i.e., fixed step size in wiggle tracks
    - $n$: number of loci in the given chromosome, defined by $n \approx |\mathscr{L}|/L$.
    - $\mathcal{S}(i)$: score for $i$th locus:
        $$c_1 g_1(i) - c_2 g_2(i) + c_3 g_3(i)$$
    - $\ell \in [0,1]^n$ (relaxed) or $\ell \in \mathbb{Z}^n_{0,1}$ (unrelaxed): solutions to the optimization problem, $\ell$, specify which loci are accessible ($\ell_i = 1$) and which are closed ($\ell_i = 0$)
    - objective function $f$: $f(\mathbf{\ell}) = -\mathcal{S}^{T}\mathbf{\ell} + \gamma \sum_{i=1}^{n-1}|\ell_i - \ell_{i+1}|$
        - the first term represents the sum of scores, $\sum_i \mathcal{S}(i)$, over selected loci
        - the second term induces sparsity and controls fragmentation of selected regions
    - budget constraint: $\sum_{i=1}^{n}\ell_i \leq [nb]$. Upper bounds the proportion of loci selected as accessible.
    - Relaxed optimization problem to obtain initial solution:
        $$\text{Minimize: } \ell \in [0,1]^n, f(\mathbf{\ell}) = -\mathcal{S}^{T}\mathbf{\ell} + \gamma \sum_{i=1}^{n-1}|\ell_i - \ell_{i+1}|$$
        $$\text{Subject to: }\sum_{i=1}^{n}\ell_i \leq [nb]$$
    - `RR`: iterative randomization procedure described in the paper to derive integral solutions (open chromatin annotations) from the relaxed solution as:
        $$\ell^{\textsf{rand}}_i \sim \text{Bernoulli}(\ell_i)$$
        `--rr_iter` such solutions are generated, after which the best feasible solution is selected to determine the final annotation.

"""

#!/usr/bin/env python
import os
import argparse
import subprocess
from . import rocco_aux
from . import chrom
from . import gwide
from . import prep
from . import budgets
from . import locus
from . import loci
import pandas as pd

def parse_coldata(coldata_file: str, group_column: str, sample_column: str, split_sex: str = None, delimiter='\t'):
    """
    Parses a metadata file for the samples. Useful if running on subgroups separately.

    Args:
        coldata_file (str): path to the colData file in CSV format
        group_column (str):  column containing group label
        sample_column (str): column containing sample names
        split_sex (str): if not None, create separate text files for each sex within each group
            based on the specified column name.
        delimiter (str): delimiter used in coldata file

    Returns:
        List: a list of the file names created
    """
    created_files = []
    try:
        coldata_df = pd.read_csv(coldata_file, delimiter=delimiter)
        group_names = coldata_df[group_column].unique()
        for group_name in group_names:
            group_samples = coldata_df[coldata_df[group_column] == group_name][sample_column]
            file_name = f'group_{group_name}'
            if split_sex is None:
                print(f'rocco.parse_coldata(): writing {file_name}')
                with open(file_name, 'w') as file:
                    file.write('\n'.join(group_samples))
                    file.write('\n')
                created_files.append(file_name)
            if split_sex is not None:
                for sex in coldata_df[split_sex].unique():
                    sex_samples = coldata_df[
                        (coldata_df[group_column] == group_name) &
                        (coldata_df[split_sex] == sex)
                    ][sample_column]
                    sex_file_name = f'group_{group_name}_sex_{sex}'
                    print(f'rocco.parse_coldata(): writing {sex_file_name}')
                    with open(sex_file_name, 'w') as file:
                        file.write('\n'.join(sex_samples))
                        file.write('\n')
                    created_files.append(sex_file_name)
    except Exception as e:
        print(f"rocco.parse_coldata(): Parsing coldata failed. Ensure column names are correctly specified.")
        raise
    return created_files

def subcommand_chrom(args):
    chrom.main(args)
    return True

def subcommand_gwide(args):
    if args['coldata'] is not None:
        ind_files = parse_coldata(args['coldata'],
                                  args['group_column'],
                                  args['sample_column'],
                                  args['split_sex'])
        bed_files = []
        for file_ in ind_files:
            temp_args = args.copy()
            print(f'rocco gwide: {file_}')
            temp_args['identifiers'] = file_
            temp_args['outdir'] = file_ + '_dir'
            temp_args['combine'] = file_ + '.bed'
            bed_files.append(temp_args['combine'])
            gwide.main(temp_args)

        print(f'rocco gwide completed: {bed_files}')
    else:
        gwide.main(args)

    return True

def subcommand_prep(args):
    prep.main(args)
    return True

def subcommand_budgets(args):
    budgets.main(args)
    return True


def main():
    parser = argparse.ArgumentParser(description="ROCCO: [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization.\n\nPyPI : https://pypi.org/project/rocco/ \nGitHub: https://github.com/nolan-h-hamilton/ROCCO/ \nDemo: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb \nPaper: https://doi.org/10.1093/bioinformatics/btad725 \n", add_help=True, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest="command")

    # 'gwide' subcommand parameters
    parser_subcommand_gwide = subparsers.add_parser("gwide", help='run rocco genome-wide/on multiple chromosomes (gwide.py)')
    parser_subcommand_gwide.add_argument('-p', '--param_file', help='Path to the parameter file containing per-chromosome parameters.\
                        \nSee https://github.com/nolan-h-hamilton/ROCCO/blob/main/rocco/hg38_params.csv for an example with chromosome-specific budgets on hg38.\
                        \n`--param_file hg_params` can be used for suggested chromosome-specific parameters in human samples.\n')
    parser_subcommand_gwide.add_argument('-b', '--budget', type=float, default=.035)
    parser_subcommand_gwide.add_argument('-g', '--gamma', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('-t', '--tau', type=float, default=0.0)
    parser_subcommand_gwide.add_argument('--c1', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('--c2', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('--c3', type=float, default=1.0)
    parser_subcommand_gwide.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_gwide.add_argument('--solver', type=str, default="CLARABEL")
    parser_subcommand_gwide.add_argument('--bed_format', type=int, default=3)
    parser_subcommand_gwide.add_argument('--identifiers', default=None)
    parser_subcommand_gwide.add_argument('--outdir', default=None, help='if `None`, a temporary directory will be created')
    parser_subcommand_gwide.add_argument('--combine', default=None)
    parser_subcommand_gwide.add_argument('--multi', default=1, type=int, help='number of simultaneous `rocco chrom` jobs to execute')
    parser_subcommand_gwide.add_argument('--verbose', default=False, action="store_true")
    parser_subcommand_gwide.add_argument('--coldata', default=None, help='if not None, parse coldata file to create group-specific `--identifiers` files on which to run rocco')
    parser_subcommand_gwide.add_argument('--group_column', default='group', help='column in coldata file containing group labels. Only used if --coldata is not None')
    parser_subcommand_gwide.add_argument('--sample_column', default='sample', help='column in coldata file containing sample labels. Only used if --coldata is not None')
    parser_subcommand_gwide.add_argument('--split_sex', default=None, help='column in coldata file containing sex labels. Only used if --coldata is not None')
    parser_subcommand_gwide.add_argument('--tracks_path', default=None)
    parser_subcommand_gwide.add_argument('--exclude_chroms', default=None)
    parser_subcommand_gwide.add_argument('--fixedStep', default=False, action="store_true")
    
    # 'chrom' subcommand parameters
    parser_subcommand_chrom = subparsers.add_parser("chrom", help='run ROCCO on a single chromosome (chrom.py)')
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
    parser_subcommand_chrom.add_argument('--solver', type=str,default="CLARABEL")
    parser_subcommand_chrom.add_argument('--bed_format', type=int, default=6)
    parser_subcommand_chrom.add_argument('--scale_bedscores', action='store_true', default=False)
    parser_subcommand_chrom.add_argument('--identifiers', default=None)
    parser_subcommand_chrom.add_argument('--outdir', type=str, default='.')
    parser_subcommand_chrom.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_chrom.add_argument('--verbose', default=False, action="store_true")
    parser_subcommand_chrom.add_argument('--fixedStep', default=False, action="store_true")

    # 'prep' subcommand parameters
    parser_subcommand_prep = subparsers.add_parser("prep", help='Preprocess BAM files (prep.py)')
    parser_subcommand_prep.add_argument('-i', '--bamdir', default='.', type=str)
    parser_subcommand_prep.add_argument('-o', '--outdir', default='.')
    parser_subcommand_prep.add_argument('-s', '--sizes', default='hg38')
    parser_subcommand_prep.add_argument('-L', '--interval_length', default=50)
    parser_subcommand_prep.add_argument('-c', '--cores', type=int, default=1)
    parser_subcommand_prep.add_argument('--multi', default=True)
    parser_subcommand_prep.add_argument('--bstw_path', default=None, help='deprecated')

    # 'budgets' subcommand parameters
    parser_subcommand_budgets = subparsers.add_parser("budgets", help='Compute a budget for each chromosome ranked by read density (budgets.py). ')
    parser_subcommand_budgets.add_argument('--wigdir', type=str)
    parser_subcommand_budgets.add_argument('-s', '--sizes', type=str)
    parser_subcommand_budgets.add_argument('--smean', type=float, default=.035)
    parser_subcommand_budgets.add_argument('--samp_rate', type=float, default=0.10)
    parser_subcommand_budgets.add_argument('-o', '--outfile', type=str, default='params.csv')
    parser_subcommand_budgets.add_argument('-L','--locsize',type=int,default=50)

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
