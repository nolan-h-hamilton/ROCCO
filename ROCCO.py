"""
Runs multiple `ROCCO_chrom.py` jobs with chromosome-specific parameters given in
a CSV file, e.g., `ROCCO/params.csv`.
"""
import os
import argparse
import subprocess
import tempfile


def get_params(param_file, budget, gamma, tau, c1, c2, c3):
    """Grabs parameters for each chromosome from `param_file`

    This function collects parameters from each chromosome-row of
    `params.csv`. If a `NULL` entry is encountered, it is replaced
    by the CLI-specified default (see `args`).

    Args:
        sizes_file (str) : a chromosome sizes filepath.


    Returns:
        chrom_params: a list of lists, with each element containing
          the necessary parameters to run a `ROCCO_chrom.py` job.
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
                    line[i] = defaults[i]
            chrom_params.append(line)
    return chrom_params


def call_rocco(chrom, wig_path, budget, gamma, tau, c1, c2, c3, solver,
               bed_format, verbose=False, N=50, identifiers=None,
               outdir='.'):
    """Formats a command to run `ROCCO_chrom.py` for a given chromosome

    Args:
        chrom (str) : chromosome on which to run `ROCCO_chrom.py`
        wig_path (str): directory containing wig files
        verbose (bool): `True` yields detailed output during workflow.
    """
    cli_args = ["python3", os.getcwd() + '/' + "ROCCO_chrom.py",
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


def sort_combine_bed(outfile, dir_='.'):
    """sorts and combines chromosome-specific bedfiles"""
    if dir_[-1] == '/':
        dir_ = dir_[0:-1]
    filenames = [f_ for f_ in os.listdir(dir_)
                 if f_.split('.')[-1] == 'bed' and 'ROCCO_out' in f_]
    filenames = sorted(filenames,
                       key=lambda x: int(val) if (val := x.split('_')[2][3:]).isnumeric() else ord(val))
    filenames = [dir_ + '/' + fname for fname in filenames]
    with open(outfile, mode='w', encoding='utf-8') as outfile_:
        for fname in filenames:
            cat_process = subprocess.Popen(('cat', fname),
                                           stdout=outfile_.fileno())
            cat_process.wait()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--param_file',
                        default='params.csv',
                        help="CSV param file w/ row for each chromosome")
    # if some parameters are left `NULL` in `param_file`,
    # the values given by the following arguments will be
    # used as defaults.
    parser.add_argument('-b', '--budget', default=.035)
    parser.add_argument('-g', '--gamma', default=1.0)
    parser.add_argument('-t', '--tau', default=0.0)
    parser.add_argument('--c1', default=1.0)
    parser.add_argument('--c2', default=1.0)
    parser.add_argument('--c3', default=1.0)
    parser.add_argument('-N', '--rr_iter', type=int, default=50)
    parser.add_argument('--solver', default="ECOS",
                        help="Optimization software used to solve the \
                        LP underlying ROCCO. `ECOS` is used by default \
                        and is a viable open-source option. `MOSEK`\
                        offers significantly greater speed and is free\
                        for academic use. Free trial commerical licenses\
                        are also available. See\
                        https://www.mosek.com/products/academic-licenses/")
    parser.add_argument('--bed_format', type=int, default=3,
                        help="Specifies BED3 or BED6 format.\
                        Default is BED6. Generate BED3 output with \
                        --bed_format 3")
    parser.add_argument('--identifiers', default=None,
                        help="an optional filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a substring of the `.wig` sample")
    parser.add_argument('--outdir', default='.',
                        help="directory in which to store ROCCO's output\
                          files.")
    parser.add_argument('--combine', default=None, help="if not None, combine\
                        output bed files and store in the file specified\
                        with this parameter. ex: `--combine bedname.bed` con-\
                        catenates the chromosome-specific bedfiles into `bedname.bed`.")
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help="number of chromosome-specific jobs to run in parallel.\
                            If set to `0`, the maximum possible jobs are created.")
    parser.add_argument('-m', '--mem', default=None,
                        help="if running in parallel, wait to start a ROCCO_chrom.py job\
                            until there is at least `-m` memory avilable. See GNU Parallel's `--memfree`\
                            parameter.")
    parser.add_argument(
        '--parlog',
        default='ROCCO_parlog.txt',
        help='logfile for GNU Parallel')
    parser.add_argument('--verbose', default=False, action="store_true")
    args = vars(parser.parse_args())

    chrom_args = get_params(args['param_file'], args['budget'],
                            args['gamma'], args['tau'], args['c1'],
                            args['c2'], args['c3'])
    if args['jobs'] == 0:
        print("running maximum possible ROCCO_chrom.py jobs in parallel")
    elif args['jobs'] == 1:
        print("running each ROCCO_chrom.py job sequentially")
    else:
        print("running min({}, num. chroms) jobs in parallel".format(
                args['jobs']))

    if not os.path.exists(args['outdir']):
        os.mkdir(args['outdir'])

    tmp = tempfile.NamedTemporaryFile(mode="w+", delete=False)
    for i, arglist in enumerate(chrom_args):
        arglist = [str(x) for x in arglist]
        cmd = call_rocco(arglist[0], arglist[1], arglist[2], arglist[3],
                         arglist[4], arglist[5], arglist[6], arglist[7],
                         args['solver'], str(args['bed_format']),
                         args['verbose'], str(args['rr_iter']),
                         identifiers=args['identifiers'], outdir=args['outdir'])
        print("job {}: {}".format(i, cmd))

        if args['jobs'] == 1:
            seq_process = subprocess.run(cmd.split(' '),
                                         capture_output=True, text=True, check=True)
            print(seq_process.stdout)

        tmp.write(str(cmd + '\n'))

    tmp.flush()
    cmd = ['parallel',
           '-k',
           '--joblog', args['parlog'],
           '--jobs', str(args['jobs']),
           '-a', tmp.name]

    if args['mem'] is not None:
        cmd += ['--memfree', args['mem']]
    if args['jobs'] != 1:
        print(' '.join(cmd))
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as par_proc:
            if args['verbose']:
                for line in par_proc.stdout:
                    print(line, end='')
            par_proc.wait()
    tmp.close()

    if args['combine'] is not None:
        print('combining output files --> {}'.format(args['combine']))
        sort_combine_bed(args['combine'], dir_=args['outdir'])


if __name__ == "__main__":
    main()
