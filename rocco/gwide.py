"""
Runs multiple [`rocco chrom`](https://nolan-h-hamilton.github.io/ROCCO/rocco/chrom.html) jobs,
which execute ROCCO on individual chromosomes. Chromosome-specific parameters are collected in a
CSV file specified with argument `-p --param_file`. See example `--param_file` below for hg38.

If any parameter entry in the CSV file `--param_file` is set to `NULL`, the script uses the corresponding default value
specified with the command-line arguments.

Parameters:
    -p, --param_file (str, required): Path to the parameter file containing per-chromosome parameters. Required.
        See https://github.com/nolan-h-hamilton/ROCCO/blob/main/hg38_params.csv for an example with chromosome-specific budgets on hg38.
        `rocco budgets` can also be used to compute chromosome-specific budgets and produces a valid `-p/--param_file`.
    -b, --budget (float): Budget parameter (upper-bound on the fraction of base pairs selected as 'open' in a given chrom.) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 0.035).
    -g, --gamma (float): Gamma parameter (discontinuity penalty weight) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 1.0).
    -t, --tau (float): Tau parameter (enrichment threshold) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 0.0).
    --c1 (float): g_1 coefficient in score function (enrichment reward) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 1.0).
    --c2 (float): g_2 coefficient in score function (dispersion penalty) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 1.0).
    --c3 (float): g_3 coefficient in score function (local shift) used for each chromosome with a `NULL` entry observed in `--param_file` (default: 1.0).
    -N, --rr_iter (int): Number of RR iterations (default: 50).
    --solver (str): Optimization software used to solve the main LP. `ECOS` is used by default (default: "ECOS").
    --bed_format (int): Specifies BED3 or BED6 format. Default is BED6. Generate BED3 output with --bed_format 3 (default: 6).
    --identifiers (str): (Optional) a filename containing identifiers for samples to include in the experiment. Each identifier should be a uniquely-identifying substring of the respective `.wig` sample. If not specified, all samples are used (default: None).
    --outdir (str): Directory in which to store output bed files from the calls to rocco chrom (default: current directory).
    --combine (str): If not None, combine output bed files and store in the file specified with this parameter. For example, `--combine bedname.bed` concatenates the chromosome-specific bedfiles into `bedname.bed` (default: None).
    --multi (int): Run `--multi` rocco chrom jobs simultaneously to improve speed. May increase peak memory use (default: 1).
    --verbose (bool): Set to `True` for verbose logging (default: False).
    --coldata (str): If not None, parse coldata file to create group-specific `--identifiers` files on which to run rocco (default: None).
    --group_column (str): Column in coldata file containing group labels. Only used if --coldata is not None (default: 'group').
    --sample_column (str): Column in coldata file containing sample labels. Only used if --coldata is not None (default: 'sample').
    --split_sex (str): Column in coldata file containing sex labels. Only used if --coldata is not None (default: None).


`hg38_params.csv` with chromosome-specific budgets:
    ```
    chromosome,input_path,budget,gamma,tau,c1,c2,c3
    chr1,tracks_chr1,0.035,NULL,NULL,NULL,NULL,NULL
    chr2,tracks_chr2,0.03,NULL,NULL,NULL,NULL,NULL
    chr3,tracks_chr3,0.03,NULL,NULL,NULL,NULL,NULL
    chr4,tracks_chr4,0.02,NULL,NULL,NULL,NULL,NULL
    chr5,tracks_chr5,0.03,NULL,NULL,NULL,NULL,NULL
    chr6,tracks_chr6,0.035,NULL,NULL,NULL,NULL,NULL
    chr7,tracks_chr7,0.035,NULL,NULL,NULL,NULL,NULL
    chr8,tracks_chr8,0.03,NULL,NULL,NULL,NULL,NULL
    chr9,tracks_chr9,0.03,NULL,NULL,NULL,NULL,NULL
    chr10,tracks_chr10,0.035,NULL,NULL,NULL,NULL,NULL
    chr11,tracks_chr11,0.04,NULL,NULL,NULL,NULL,NULL
    chr12,tracks_chr12,0.04,NULL,NULL,NULL,NULL,NULL
    chr13,tracks_chr13,0.025,NULL,NULL,NULL,NULL,NULL
    chr14,tracks_chr14,0.03,NULL,NULL,NULL,NULL,NULL
    chr15,tracks_chr15,0.035,NULL,NULL,NULL,NULL,NULL
    chr16,tracks_chr16,0.04,NULL,NULL,NULL,NULL,NULL
    chr17,tracks_chr17,0.055,NULL,NULL,NULL,NULL,NULL
    chr18,tracks_chr18,0.025,NULL,NULL,NULL,NULL,NULL
    chr19,tracks_chr19,0.06,NULL,NULL,NULL,NULL,NULL
    chr20,tracks_chr20,0.04,NULL,NULL,NULL,NULL,NULL
    chr21,tracks_chr21,0.03,NULL,NULL,NULL,NULL,NULL
    chr22,tracks_chr22,0.04,NULL,NULL,NULL,NULL,NULL
    chrX,tracks_chrX,0.02,NULL,NULL,NULL,NULL,NULL
    chrY,tracks_chrY,0.01,NULL,NULL,NULL,NULL,NULL
    ```

Examples [from demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb):
    ```
    rocco gwide -p ../hg38_params.csv --outdir demo_outdir --combine demo_out.bed
    rocco gwide -p ../hg38_params.csv --coldata coldata.csv --sample_column sample --group_column group
    ```
"""

import os
import sys
import argparse
import subprocess
import tempfile
import pandas as pd
from . import rocco_aux

def get_params(param_file: str, budget: float, gamma: float, tau: float,
               c1: float, c2: float, c3: float) -> list:
    """
    Grabs parameters for each chromosome from `param_file`

    This function collects parameters from each chromosome-row of
    `params.csv`. If a `NULL` entry is encountered, it is replaced
    by the CLI-specified default (see `args`).

    Args:
        param_file (str) : a CSV file with a row containing parameter vals\
            to use for each chromosome.

    Returns:
        list: a list of lists, with each element containing
          the necessary parameters to run a `rocco chrom` job.
    """
    defaults = [None, None, budget, gamma, tau, c1, c2, c3]
    chrom_params = []
    with open(param_file, mode='r', encoding='utf-8') as par_file:
        header = True
        for line in par_file:
            if header is True:
                header = False
                continue
            if ',' not in line:
                continue
            line = line.strip()
            line = line.split(',')
            for i, entry in enumerate(line):
                if entry == 'NULL':
                    # replace NULL entries with the genome-wide default
                    line[i] = defaults[i]
            chrom_params.append(line)
    return chrom_params


def call_rocco(chrom, wig_path, budget, gamma, tau, c1, c2, c3, solver,
               bed_format, verbose=False, N=50, identifiers=None,
               outdir='.') -> str:
    r"""
    Formats a command to run `rocco chrom` for a given chromosome

    Args:
        chrom (str): chromosome name. Example: `chr1`
        wig_path (str): path to the *directory* containing samples' signal
            tracks for `chrom`
        budget (float): budget constraint
        tau (float): tau in $\mathcal{S}(i)$
        gam (float): gamma parameter in $f(\mathbf{\ell})$
        c1 (float): c1 value in $\mathcal{S}(i)$
        c2 (float): c2 value in $\mathcal{S}(i)$
        c3 (float): c3 value in $\mathcal{S}(i)$
        N (int): RR procedure iterations. If N <= 0,
            $\texttt{floor\_eps}$ procedure is applied...
            See `Loci.run_rr()`
        verbose_ (bool): Verbosity flag for the solver
        solver (str): the solver to use--either "ECOS" or "MOSEK"
    Returns:
        str: a formatted rocco chrom command with the appropriate parameters.
        substituted.
    """
    cli_args = ["rocco", "chrom",
                '--chrom', chrom,
                '--wig_path', wig_path,
                '--budget', budget,
                '--gamma', gamma,
                '--tau', tau,
                '--c1', c1,
                '--c2', c2,
                '--c3', c3,
                '--solver', solver,
                '--bed_format', bed_format,
                '--outdir', outdir,
                '--rr_iter', N,
                '--identifiers', identifiers,
                '--verbose']

    if not verbose:
        cli_args.remove('--verbose')
    if identifiers is None:
        cli_args.remove('--identifiers')
        cli_args.remove(identifiers)
    return ' '.join(cli_args)

def main(args):
    chrom_args = get_params(args['param_file'], args['budget'],
                            args['gamma'], args['tau'], args['c1'],
                            args['c2'], args['c3'])

    args['outdir'] = rocco_aux.trim_path(args['outdir'])
    if not os.path.exists(args['outdir']):
        os.mkdir(args['outdir'])

    tmp = tempfile.NamedTemporaryFile(mode="w+")
    for i, arglist in enumerate(chrom_args):
        arglist = [str(x) for x in arglist]
        if not os.path.exists(arglist[1]) or not len(os.listdir(arglist[1])) > 0:
            print(f'directory {arglist[1]} does not exist or is empty...skipping')
            continue
        cmd = call_rocco(arglist[0], arglist[1], arglist[2], arglist[3],
                         arglist[4], arglist[5], arglist[6], arglist[7],
                         args['solver'], str(args['bed_format']),
                         args['verbose'], str(args['rr_iter']),
                         identifiers=args['identifiers'], outdir=args['outdir'])

        if args['multi'] == 1:
            try:
                seq_process = subprocess.run(cmd.split(' '),
                                         capture_output=True, text=True, check=True)
                print(seq_process.stdout)
            except Exception as ex:
                print(ex)
                raise ex

        tmp.write(str(cmd + '\n'))

    tmp.flush()
    if args['multi'] > 1:
        rocco_aux.run_par(tmp.name, threads=args['multi'], verbose=args['verbose'])
    tmp.close()

    if args['combine'] is not None:
        print('combining output files --> {}'.format(args['combine']))
        rocco_aux.sort_combine_bed(args['combine'], dir_=args['outdir'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--param_file',
                        required=True, help='Path to the parameter file containing per-chromosome parameters. Required.\
                        \nSee https://github.com/nolan-h-hamilton/ROCCO/blob/main/hg38_params.csv for an example with chromosome-specific budgets on hg38.\n\
                        `rocco budgets` can also be used to compute chromosome-specific budgets given the input wig files and produces a valid -p/--param_file.')
    parser.add_argument('-b', '--budget', type=float, default=.035, help='budget parameter (largest allowed fraction of selected bp) used for each chromosome with a `NULL` entry observed in `--param_file`')
    parser.add_argument('-g', '--gamma', type=float, default=1.0, help='gamma parameter (discontig. penalty weight) used for each chromosome  with a `NULL` entry observed in `--param_file`')
    parser.add_argument('-t', '--tau', type=float, default=0.0, help='tau parameter (enrichment threshold) used for each chromosome  with a `NULL` entry observed in `--param_file`')
    parser.add_argument('--c1', type=float, default=1.0, help='g_1 coefficient in score function (enrichment reward) used for each chromosome  with a `NULL` entry observed in `--param_file`')
    parser.add_argument('--c2', type=float, default=1.0, help='g_2 coefficient in score function (dispersion penalty) used for each chromosome  with a `NULL` entry observed in `--param_file`')
    parser.add_argument('--c3', type=float, default=1.0, help='g_3 coefficient in score function (local shift) used for each chromosome  with a `NULL` entry observed in `--param_file`')
    parser.add_argument('-N', '--rr_iter', type=int, default=50, help = 'number of RR iterations')
    parser.add_argument('--solver', default="ECOS",
                        help="Optimization software used to solve the \
                        main LP. `ECOS` is used by default.")
    parser.add_argument('--bed_format', type=int, default=3,
                        help="Specifies BED3 or BED6 format.\
                        Default is BED6. Generate BED3 output with \
                        --bed_format 3")
    parser.add_argument('--identifiers', default=None,
                        help="(optional) a filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a uniquely-identifying substring of\
                          the respective `.wig` sample. If not specified, all\
                          samples are used.")
    parser.add_argument('--outdir', default='.',
                        help="directory in which to store output bed files from the calls to rocco chrom")
    parser.add_argument('--combine', default=None, help="if not None, combine\
                        output bed files and store in the file specified\
                        with this parameter. ex: `--combine bedname.bed` con-\
                        catenates the chromosome-specific bedfiles into `bedname.bed`.")
    parser.add_argument('--multi', type=int, default=1, help='run `--multi` rocco chrom jobs\
        simultaneously to improve speed. May increase peak memory use.')
    parser.add_argument('--verbose', default=False, action="store_true")
    parser.add_argument('--coldata', default=None, help='if not None, parse coldata file to create group-specific `--identifiers` files on which to run rocco')
    parser.add_argument('--group_column', default='group', help='column in coldata file containing group labels. Only used if --coldata is not None')
    parser.add_argument('--sample_column', default='sample', help='column in coldata file containing sample labels. Only used if --coldata is not None')
    parser.add_argument('--split_sex', default=None, help='column in coldata file containing sex labels. Only used if --coldata is not None')

    args = vars(parser.parse_args())
    main(args)