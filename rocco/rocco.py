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
    gwide               run rocco genome-wide/on multiple chromosomes (gwide.py)
    chrom               run ROCCO on a single chromosome (chrom.py)
    prep                Preprocess BAM files (prep.py)
    budgets             Compute a budget (upper-bound on the fraction of basepairs that can be selected as 'open') for
                        each chromosome ranked by average read-density observed in samples (budgets.py)
    get_sizes           download sizes file for a genome in the ucsc genome registry
options:
  -h, --help            show this help message and exit

```

Subcommand Documentation:
    [`gwide`](https://nolan-h-hamilton.github.io/ROCCO/rocco/gwide.html)

    [`chrom`](https://nolan-h-hamilton.github.io/ROCCO/rocco/chrom.html)

    [`prep`](https://nolan-h-hamilton.github.io/ROCCO/rocco/prep.html)

    [`budgets`](https://nolan-h-hamilton.github.io/ROCCO/rocco/budgets.html)

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
    Given a metadata file for the samples, create text files containing sample names for each group

    These created textfiles are then used for the `--identifiers` argument of `rocco gwide`/`rocco chrom`

    *The entries in `sample_column` should uniquely identify the corresponding `.wig` files* generated with `rocco prep`.
        Ideally, they are just the stripped names of the corresponding BAM files supplied to `rocco prep`.

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
        print(coldata_df)
        group_names = coldata_df[group_column].unique()
        for group_name in group_names:
            group_samples = coldata_df[coldata_df[group_column] == group_name][sample_column]
            print(group_samples)
            file_name = f'group_{group_name}.txt'
            if split_sex is None:
                print(f'writing {file_name}')
                with open(file_name, 'w') as file:
                    file.write('\n'.join(group_samples))
                created_files.append(file_name)
            if split_sex is not None:
                for sex in coldata_df[split_sex].unique():
                    sex_samples = coldata_df[
                        (coldata_df[group_column] == group_name) &
                        (coldata_df[split_sex] == sex)
                    ][sample_column]
                    sex_file_name = f'group_{group_name}_sex_{sex}.txt'
                    print(f'writing {sex_file_name}')
                    with open(sex_file_name, 'w') as file:
                        file.write('\n'.join(sex_samples))
                    created_files.append(sex_file_name)
    except Exception as e:
        print(f"parsing coldata failed:\n{e}.\nEnsure column names are correctly specified.")
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
        
        print(f'completed: {bed_files}')
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
    parser = argparse.ArgumentParser(description="ROCCO: [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization.\n\nPyPI : https://pypi.org/project/rocco/ \nGitHub: https://github.com/nolan-h-hamilton/ROCCO/ \nDemo: https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb \nPaper: https://doi.org/10.1101/2023.05.24.542132\n", add_help=True, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest="command")

    # 'gwide' subcommand parameters
    parser_subcommand_gwide = subparsers.add_parser("gwide", help='run rocco genome-wide/on multiple chromosomes (gwide.py)')
    parser_subcommand_gwide.add_argument('-p', '--param_file', required=True,
                                         help='A CSV file specifying the parameters for each chromosome. See: https://github.com/nolan-h-hamilton/ROCCO/blob/main/hg38_params.csv')
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
    parser_subcommand_gwide.add_argument('--coldata', default=None, help='if not None, parse coldata file to create group-specific `--identifiers` files on which to run rocco')
    parser_subcommand_gwide.add_argument('--group_column', default='group', help='column in coldata file containing group labels. Only used if --coldata is not None')
    parser_subcommand_gwide.add_argument('--sample_column', default='sample', help='column in coldata file containing sample labels. Only used if --coldata is not None')
    parser_subcommand_gwide.add_argument('--split_sex', default=None, help='column in coldata file containing sex labels. Only used if --coldata is not None')

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
    parser_subcommand_chrom.add_argument('--solver', type=str,default="ECOS")
    parser_subcommand_chrom.add_argument('--bed_format', type=int, default=6)
    parser_subcommand_chrom.add_argument('--identifiers', default=None)
    parser_subcommand_chrom.add_argument('--outdir', type=str, default='.')
    parser_subcommand_chrom.add_argument('-N', '--rr_iter', type=int, default=50)
    parser_subcommand_chrom.add_argument('--verbose', default=False, action="store_true")

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
