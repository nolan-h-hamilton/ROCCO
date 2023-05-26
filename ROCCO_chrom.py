"""
Construct $mathbf{S}_{chr}$, score loci, and solve the optimization problem
underlying ROCCO. This script looks in `--wig_path` for smooth signal tracks
for each replicate.
"""
import os
import argparse
import warnings
import subprocess
import numpy as np
import pandas as pd
from locus import Locus
from loci import Loci


def get_locus_size(wig_path):
    """Infers interval size from replicates` wig files

    Args:
        wig_path (str) : path to directory containing wig files


    Returns:
        an integer locus/interval size used in the wig files
    """
    i = 0
    pos1 = 0
    pos2 = 0
    wig_file = [x for x in os.listdir(wig_path) if '.wig' in x][0]
    with open(wig_path + '/' + wig_file, mode='r', encoding='utf-8') as wig:
        for line in wig:
            line = line.strip()
            if line.replace('\t', '').replace(' ', '').isnumeric():
                if i == 0:
                    pos1 = int(line.split('\t')[0])
                if i == 1:
                    pos2 = int(line.split('\t')[0])
                if pos2 - pos1 > 0:
                    break
                i += 1
    return pos2 - pos1


def get_start_end(wig_path):
    """Infers starting and ending nucleotide position common to the wigs

    To infer a reasonable genomic region over which to run ROCCO,
    this function returns the median starting point and endpoints
    among the multiple replicates' signal files in `wig_path`.

    Args:
        wig_path (str) : path to directory containing wig files

    Returns:
        start (int): median starting position
        end (int): median ending position
    """
    starts = []
    ends = []

    def get_start(track):
        with open(track, mode='r', encoding='utf-8') as wig:
            for line in wig:
                line = line.strip()

                if line.replace('\t', '').replace(' ', '').isnumeric():
                    start = int(line.split('\t')[0])
                    return start
        return None

    def get_end(track):
        with subprocess.Popen(["tail", "-n1", track],
                              stdout=subprocess.PIPE) as proc:
            for line in proc.stdout:
                line = line.decode()
                line = line.strip()
                return int(line.split('\t')[0])
        return None

    wig_files = [x for x in os.listdir(wig_path) if '.wig' in x]
    for wig_file in wig_files:
        starts.append(get_start(wig_path + '/' + wig_file))
        ends.append(get_end(wig_path + '/' + wig_file))
    return sorted(starts)[len(starts) // 2], sorted(ends)[len(ends) // 2]


def read_wig(wig_file, start=0, end=10**10, locus_size=50):
    """Processes a wig file for inclusion in the signal matrix S_chr

    This function prepares signal data for to fit in a constant-sized
    Kxn matrix. If a wig file begins at a position < `start`, it will
    be padded with zeros. Likewise if a wig file ends at a position >
    `end`, the extra values will be discarded. Gaps between start/end
    are replaced with a zero-entry.

    Args:
        wig_file (str): path to wiggle formatted signal track
        start (int): inferred or manually-specified starting nucleotide
          position
        end (int): inferred or manually-specified ending nucleotide
          position
        locus_size (int): interval/locus size. Each wig file will have
          signal values at every start + i*`locus_size` nucleotide pos.
          for i = 1,2,...

    Returns:
        loci: list of starting nucleotide positions for each locus
        signal: enrichment signal value at each locus in `loci`
    """
    log("ROCCO_chrom: reading wig file {}".format(wig_file))
    loci = []
    signal = []
    with open(wig_file, mode='r', encoding='utf-8') as wig:
        for line in wig:
            line = line.strip()
            try:
                line = line.split('\t')
                line[0] = int(line[0])
                line[1] = float(line[1])
                if line[1] < 0:
                    line[1] = 0
                # case: wig begins after specified `start`
                if start < line[0] and len(loci) == 0:
                    for loc in range(start, line[0], locus_size):
                        loci.append(loc)
                        signal.append(0)
                    loci.append(line[0])
                    signal.append(line[1])
                    continue
                if line[0] < start:
                    continue
                if line[0] > end:
                    break
                # case: loci not contiguous
                if len(loci) > 0 and line[0] - loci[-1] > locus_size:
                    loci.append(loci[-1] + locus_size)
                    signal.append(0)
                    continue

                loci.append(line[0])
                signal.append(line[1])

            except ValueError:
                continue
            except IndexError:
                continue
    # case: wig ends prematurely
    if loci[-1] < end:
        for loc in range(loci[-1] + locus_size, end + locus_size, locus_size):
            loci.append(loc)
            signal.append(0)

    return loci, signal


def collect_wigs(wig_files, start=0, end=10**10, locus_size=50):
    """creates a Kxn dataframe of signals derived from wig files"""
    loci = [x for x in range(start, end + locus_size, locus_size)]
    all_wigs = pd.DataFrame(
        np.array([read_wig(x, start, end, locus_size)[1]
                  for x in wig_files]),
        columns=loci, index=None)
    return all_wigs


def emp_cdf(val, arr):
    """Get empirical cdf of `val` in `arr`

    Note: Assumes `arr` is sorted!


    Args:
        val (float): a numeric element in `arr`
        arr (array): a SORTED array of peak scores

    Returns:
        (float) estimate for percentile of `val` in `arr`
    """
    return np.searchsorted(arr, val) / len(arr)


def sig_selected(loci_obj):
    """Return vector of scores from all selected, unmerged loci"""
    sig_ref = []
    for loc in loci_obj:
        if loc.accessible == 1:
            sig_ref.append(np.sum(loc.sig_data) / loc.size)
    return np.array(sig_ref)


def tags(fname):
    """Get sample ID tags from newline-separated file"""
    tags_ = []
    for line in open(fname, mode='r', encoding='utf-8'):
        tags_.append(line.strip())
    return tags_


def log(text, verbose=True):
    """only print if in verbose mode"""
    if verbose:
        print(str(text))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', help="e.g., --chrom chr1")
    parser.add_argument('--start', type=int, default=-1)
    parser.add_argument('--end', type=int, default=-1)
    parser.add_argument('--locus_size', type=int, default=-1,
                        help="this must match the constant step-size \
                        in the wiggle files used as input")
    parser.add_argument('--wig_path', type=str, default=os.getcwd())
    parser.add_argument('-N', '--rr_iter', type=int, default=50)
    parser.add_argument('--verbose', default=False, action="store_true")
    parser.add_argument('--integral', type=bool, default=False)
    parser.add_argument('-b', '--budget', type=float, default=.035)
    parser.add_argument('-g', '--gamma', type=float, default=1.0)
    parser.add_argument('-t', '--tau', type=float, default=0.0)
    parser.add_argument('--c1', type=float, default=1.0)
    parser.add_argument('--c2', type=float, default=1.0)
    parser.add_argument('--c3', type=float, default=1.0)
    parser.add_argument('--solver', default="ECOS")
    parser.add_argument('--bed_format', type=int, default=6,
                        help="`3` for BED3 format and `6` for BED6 format")
    parser.add_argument('--identifiers', default=None,
                        help="(optional) a filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a substring of the `.wig` sample.")
    parser.add_argument('--outdir', type=str, default='.')

    args = vars(parser.parse_args())

    log('ROCCO_chrom: args', args['verbose'])
    log(args, args['verbose'])

    if args['locus_size'] == -1:
        args['locus_size'] = get_locus_size(args['wig_path'])

    if args['outdir'][-1] == '/':
        args['outdir'] = args['outdir'][0:-1]

    if args['start'] == -1 or args['end'] == -1:
        args['start'], args['end'] = get_start_end(args['wig_path'])
    if args['bed_format'] not in [3, 6]:
        warnings.warn('Only BED3 and BED6 are supported--setting to BED3')
        args['bed_format'] = 6

    log("ROCCO_chrom: inferred args", args['verbose'])
    log(args, args['verbose'])

    wig_files = [args['wig_path'] + '/' + fname
                 for fname in os.listdir(args['wig_path'])
                 if 'wig' in fname.split('.')[-1]]
    if args['identifiers'] is not None:
        new_wig_files = []
        for tag in tags(args['identifiers']):
            for wfname in wig_files:
                if tag in wfname:
                    new_wig_files.append(wfname)
        wig_files = new_wig_files

    signal_matrix = collect_wigs(wig_files,
                                 start=args['start'],
                                 end=args['end'],
                                 locus_size=args['locus_size'])

    log("ROCCO_chrom: building Loci object", args['verbose'])
    InitLoci = Loci()
    for i, loc in enumerate(signal_matrix.columns):
        size_ = args['locus_size']
        NewLoc = Locus(position=args['start'] + i * size_,
                       size=size_,
                       sig_data=signal_matrix[loc])
        InitLoci.append_locus(NewLoc)

    log("ROCCO_chrom: solving for optimal solution", args['verbose'])

    if not args['integral']:
        InitLoci.rocco_lp(budget=args['budget'],
                          gam=args['gamma'],
                          tau=args['tau'],
                          c1=args['c1'],
                          c2=args['c2'],
                          c3=args['c3'],
                          verbose_=args['verbose'],
                          solver=args['solver'],
                          N=args['rr_iter'])

    if args['integral']:
        # Optimize with integer constraints
        # integer optimization may be very slow using the open-source
        # default solver, ECOS_BB, consider using the MOSEK solver
        # for problems with integer constraints
        InitLoci.rocco_ip(budget=args['budget'],
                          gam=args['gamma'],
                          tau=args['tau'],
                          c1=args['c1'],
                          c2=args['c2'],
                          c3=args['c3'],
                          verbose_=args['verbose'],
                          solver=args['solver'])

    log("ROCCO_chrom: merging adjacent selections", args['verbose'])
    InitLoci.combine_selected()

    if not os.path.exists(args['outdir']):
        os.mkdir('./' + args['outdir'])

    fname = args['outdir'] + '/' + "ROCCO_out_{}_{}_{}_{}_{}_{}_{}.bed".format(
        args['chrom'],
        args['budget'],
        args['gamma'],
        args['tau'],
        args['c1'],
        args['c2'],
        args['c3'])

    log('ROCCO_chrom: writing output: {}'.format(fname), True)
    outfile = open(fname, mode='w', encoding='utf-8')
    output = ""
    if args['bed_format'] == 6:
        loc_count = 1
        sorted_sel_sig = np.sort(sig_selected(InitLoci))
        for loc in InitLoci:
            assert loc.size > 0
            if loc.accessible == 1:
                # score is used for coloring in ucsc genome browser
                loc_score = int(500 * emp_cdf(np.sum(loc.sig_data) / loc.size,
                                              sorted_sel_sig) + 500)
                loc_name = args['chrom'] + '_' + str(loc_count)
                output += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    args['chrom'], loc.position, loc.position + loc.size,
                    loc_name, loc_score, '.')
                loc_count += 1

    if args['bed_format'] == 3:
        for loc in InitLoci:
            if loc.accessible == 1:
                output += "{}\t{}\t{}\n".format(
                    args['chrom'], loc.position, loc.position + loc.size)

    outfile.write(output)
    outfile.close()


if __name__ == "__main__":
    main()
