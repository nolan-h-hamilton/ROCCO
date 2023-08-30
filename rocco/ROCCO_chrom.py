"""
Run ROCCO on singe particular chromosome `--chrom`.
This script is the workhorse for [`rocco gwide`](https://nolan-h-hamilton.github.io/ROCCO/rocco/ROCCO_gwide.html) which generates and compiles
results for multiple chromosomes

Arguments:
    --chrom (str):
        Chromosome identifier, e.g., --chrom chr1.

    --wig_path (str, default=os.getcwd()):
        Directory containing .wig files for each sample.

    -N, --rr_iter (int, default=50):
        Number of Repetitive Refinement (RR) iterations.

    --verbose (bool, default=False):
        Set to `True` for verbose logging.

    -b, --budget (float, default=0.035):
        Budget parameter (largest allowed fraction of selected bp).

    -g, --gamma (float, default=1.0):
        Gamma parameter (discontiguity penalty weight).

    -t, --tau (float, default=0.0):
        Tau parameter (enrichment threshold).

    --c1 (float, default=1.0):
        g_1 coefficient in score function (enrichment reward).

    --c2 (float, default=1.0):
        g_2 coefficient in score function (dispersion penalty).

    --c3 (float, default=1.0):
        g_3 coefficient in score function (local shift).

    --solver (str, default="ECOS"):
        Optimization software used to solve the main LP. `ECOS` is used by default.

    --bed_format (int, default=6):
        Format for output BED file. `3` for BED3 format and `6` for BED6 format.

    --identifiers (str, default=None):
        (optional) Filename containing identifiers for samples to include in the experiment. Each identifier should be a uniquely-identifying substring of the respective `.wig` sample. If not specified, all samples are used.

    --outdir (str, default='.'):
        Directory in which to store the output BED file.

"""

import os
import gc
import argparse
import warnings
import subprocess
import numpy as np
import pandas as pd
from .locus import Locus
from .loci import Loci



def get_locus_size(wig_path) -> int:
    """
    Infers step size from replicates` wig files.

    Args:
        wig_path (str) : path to *directory* containing wig files


    Returns:
        int: an integer locus/step size used in the wig files

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
    """
    Infers starting and ending nucleotide position common to the
    multiple signal tracks

    To infer a reasonable genomic region over which to run ROCCO,
    this function returns the median starting point and endpoints
    among the multiple replicates' signal files in `wig_path`.

    Args:
        wig_path (str) : path to *directory* containing wig files

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
    r"""
    Processes a wig file for inclusion in the signal matrix $\mathbf{S}_{chr}$

    This function prepares signal data for to fit in a constant-sized
    $K \times n$ matrix. If a wig file begins at a position < `start`,
    it will be padded with zeros. Likewise if a wig file ends at a
    position > `end`, the extra values are ignored. Gaps between `start`
    and `end` are replaced with a zero-entry.

    Args:
        wig_file (str): path to wiggle formatted signal track
        start (int): inferred or manually-specified starting nucleotide
          position
        end (int): inferred or manually-specified ending nucleotide
          position
        locus_size (int): interval/locus size. Each wig file will have
          signal values at every `start + i*locus_size` nucleotide pos.
          for i = 0,1,2,...

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


def collect_wigs(wig_files, start=0, end=10**10, locus_size=50) -> pd.DataFrame:
    r"""
    creates a $K \times n$ dataframe of signals derived from wig files

    Args:
        wig_files (list): list of wig filepaths
        start (int): starting nucleotide position of wig files
        end (int): ending nucleotide position of wig files
        locus_size (int): initial size (in nucleotides) of Locus objects.
            Corresponds to $L$ in paper.

    Returns:
        pd.DataFrame: $K \times n$ dataframe of individuals' signals.
    """
    loci = [x for x in range(start, end + locus_size, locus_size)]
    all_wigs = pd.DataFrame(
        np.array([read_wig(x, start, end, locus_size)[1]
                  for x in wig_files],dtype=np.float32),
        columns=loci, index=None)
    return all_wigs


def emp_cdf(val, arr) -> float:
    """
    Get empirical cdf of `val` in `arr`

    Used to compute 0-1000 scaled score for BED6 files
    and visualization in genome browser.

    Note:
        Assumes `arr` is sorted!

    Args:
        val (float): a numeric element in `arr`
        arr (array): a SORTED array of peak values

    Returns:
        (float) estimate for percentile of `val` in `arr`
    """
    return np.searchsorted(arr, val) / len(arr)


def sig_selected(loci_obj: Loci) -> np.ndarray:
    """
    Get enrichment signals from all selected loci

    Args:
        loci_obj (Loci): Loci object after optimization

    Returns:
        np.ndarray: vector of enrichment signals from all selected
            (but still unmerged) loci
    """
    sig_ref = []
    for loc in loci_obj:
        if loc.accessible == 1:
            sig_ref.append(np.sum(loc.sig_data) / loc.size)
    return np.array(sig_ref)


def tags(fname) -> list:
    """
    Get sample IDs from file specified with command-line argument `--identifiers`

    Args:
        fname (str): path to file containing a subset of sample IDs--one line for each ID

    Returns:
        list: list of sample IDs in `fname`
    """
    tags_ = []
    for line in open(fname, mode='r', encoding='utf-8'):
        tags_.append(line.strip())
    return tags_


def log(text, verbose=True) -> None:
    """
    Only print `text` if `verbose=True`.
    """
    if verbose:
        print(str(text))


def main(args):
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

    # collect wig files for each replicate/sample for `--chrom`
    wig_files = [args['wig_path'] + '/' + fname
                 for fname in os.listdir(args['wig_path'])
                 if 'wig' in fname.split('.')[-1]]

    # select a subset of the wig files if `--identifiers` is defined
    if args['identifiers'] is not None:
        new_wig_files = []
        for tag in tags(args['identifiers']):
            for wfname in wig_files:
                if tag in wfname:
                    new_wig_files.append(wfname)
        wig_files = new_wig_files

    # create a dataframe of signal values for `K` replicates and `n` loci
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
                       sig_data=np.array(signal_matrix[loc],dtype=np.float32))
        InitLoci.append_locus(NewLoc)
    del signal_matrix
    gc.collect()

    log("ROCCO_chrom: solving for optimal solution", args['verbose'])
    InitLoci.rocco_lp(budget=args['budget'],
                        gam=args['gamma'],
                        tau=args['tau'],
                        c1=args['c1'],
                        c2=args['c2'],
                        c3=args['c3'],
                        verbose_=args['verbose'],
                        solver=args['solver'],
                        N=args['rr_iter'])

    log("ROCCO_chrom: merging adjacent selections", args['verbose'])
    InitLoci.combine_selected()

    if not os.path.exists(args['outdir']):
        os.mkdir(args['outdir'])

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
                # the score can be used for coloring in ucsc genome browser
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', type=str, help="e.g., --chrom chr1")
    parser.add_argument('--start', type=int, default=-1, help='beginning nucleotide position (default: infer from wiggle files)')
    parser.add_argument('--end', type=int, default=-1, help='ending nucleotide position (default: infer from wiggle files)')
    parser.add_argument('--locus_size', type=int, default=-1, help='size of loci (L parameter) (default: infer from wiggle files)')
    parser.add_argument('--wig_path', type=str, default=os.getcwd(), help='directory containing .wig files for each sample')
    parser.add_argument('-N', '--rr_iter', type=int, default=50, help = 'number of RR iterations')
    parser.add_argument('--verbose', default=False, action="store_true", help='set to `True` for verbose logging')
    parser.add_argument('-b', '--budget', type=float, default=.035, help='budget parameter (largest allowed fraction of selected bp)')
    parser.add_argument('-g', '--gamma', type=float, default=1.0, help='gamma parameter (discontig. penalty weight)')
    parser.add_argument('-t', '--tau', type=float, default=0.0, help='tau parameter (enrichment threshold)')
    parser.add_argument('--c1', type=float, default=1.0, help='g_1 coefficient in score function (enrichment reward)')
    parser.add_argument('--c2', type=float, default=1.0, help='g_2 coefficient in score function (dispersion penalty) ')
    parser.add_argument('--c3', type=float, default=1.0, help='g_3 coefficient in score function (local shift)')
    parser.add_argument('--solver', default="ECOS", help='Optimization software used to solve the \
                        main LP. `ECOS` is used by default.')
    parser.add_argument('--bed_format', type=int, default=6,
                        help="`3` for BED3 format and `6` for BED6 format")
    parser.add_argument('--identifiers', default=None,
                        help="(optional) a filename containing identifiers\
                          for samples to include in experiment. Each identi\
                          fier should be a uniquely-identifying substring of\
                          the respective `.wig` sample. If not specified, all\
                          samples are used.")
    parser.add_argument('--outdir', type=str, default='.', help="directory in which to store output bed file")

    args = vars(parser.parse_args())
    main(args)
