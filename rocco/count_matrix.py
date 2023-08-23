"""
Create a count matrix/df with rows corresponding to samples
and columns corresponding to genomic regions (i.e., peaks
identified w/ ROCCO). Each `(i,j)` entry represents the number of
of reads mapping to the region in column `j` and sample in row `i`.

Usage:
    `count_matrix.py [-h] [-i PEAKFILE] [-m METADATA] [--bamdir BAMDIR] [-o OUTFILE]`

Arguments:
    peakfile (str): (`-i`) BED-formatted peak file
    metadata (str, optional):  (`-m`) path to sample metadata file. Assumed
            structure is a TSV file with a row for each sample with the
            leftmost column specifying the sample identifier. *If this argument
            is not invoked, this script will consider all samples with BAM files
            in `bamdir` to construct the count matrix*.
    bamdir (str): path to directory containing samples' BAM alignments.
    outfile (str): filename of generated count matrix
"""
import os
import argparse
import pandas as pd
import multiprocessing
from pybedtools import BedTool
from . import rocco_aux

def proc_bed(bedfile, header=False):
    """
    Process a BED file for input to function `count_matrix.peak_count_matrix()`

    Args:
        bedfile (str): Path to the BED file containing peak data
        header (bool): `True` if BED file contains a header,
            `False` if not. Defaults to `False`

    Returns:
        A pybedtols.BedTool object derived from `bedfile` that is formatted
            for `count_matrix.peak_count_matrix()`
    """
    i = 0
    bed_string = ""
    with open(bedfile, 'r') as bf:
        for line in bf:

            if header and i == 0:
                i += 1
                continue

            line = line.strip()
            line = line.split('\t')

            if len(line) < 3:
                raise SyntaxError('Peak/BED file is not formatted\
                    correctly')

            # if already in BED6 format, return BedTool object as is
            if len(line) == 6:
                return BedTool(bedfile)

            # if more than 6 entries, just truncate each line
            if len(line) > 6:
                new_line = "\t".join(line[0:5]) + '\n'
                bed_string += new_line
                continue

            # get first three entries of line
            new_line = "\t".join(line[0:3])

            # BED3 format, add a name, score 1000, strand `.`
            if len(line) == 3:
                # create name for peak: `<chr_start_end>`
                name = '_'.join([line[0],line[1],line[2]])
                # assign 1000 for score and . for strand
                new_line += "\t{}\t{}\t{}".format(name,'1000','.')

            # BED4, add a score 1000, strand `.`
            elif len(line) == 4:
                newline += "\t{}".format(line[3])
                new_line += "\t{}\t{}".format('1000','.')

            # BED5, add a strand `.`
            elif len(line) == 5:
                newline += "\t{}".format(line[3])
                newline += "\t{}".format(line[4])
                new_line += "\t{}".format('.')

            new_line += '\n'
            bed_string += new_line

    return BedTool(bed_string, from_string=True)


def peak_count_matrix(bedfile, metadata, bamdir='.'):
    """
    Wrapper for `pybedtools.multi_bam_coverage()` for a specified
    `bedfile` and samples' BAM files in `bamdir`. Generates a TSV
    file, that can be interpreted as a count matrix suitable for
    differential analyses.

    Args:
        bedfile (str): path to peak data file (bed formatted)
        metadata (str, optional): (`-m`) path to sample metadata file. Assumed
            structure is a TSV file with a row for each sample
            and the leftmost column specifying the sample identifier.
            *If this argument is not specified, this script will use
            all samples with BAM files in `bamdir` to construct the
            count matrix*.
        bamdir (str): path to the directory containing samples'
            BAM files. Defaults to current working directory.
    """

    bamfiles = []
    bdir = rocco_aux.trim_path(bamdir)
    # read sample `metadata` as a pandas dataframe
    try:
        bam_df = pd.read_csv(metadata,
                                header=None,
                                index_col=False,
                                skiprows=[0],
                                sep="\t")

        IDs = [x for x in list(bam_df.iloc[:,0])]

    except (ValueError,FileNotFoundError) as no_md_err:
        print(no_md_err)
        print(f'count_matrix.py: could not read metadata... attempting to\
            use BAM filenames from {os.path.abspath(bamdir)}\
            as sample IDs')
        IDs = [fname for fname in os.listdir(bamdir) if 'bam' in fname.split('.')[-1]]

    # match sample IDs in `metadata` with their corresponding
    # BAM file in `bamdir`.
    init_IDs= [x for x in IDs] # copy of original sample IDs for output
    for i,ID in enumerate(IDs):
        if not os.path.exists(ID):
            for fname in os.listdir(bdir):
                fname = f'{bdir}/{fname}'
                if ID in fname and rocco_aux.is_alignment(fname):
                    if fname in IDs:
                        continue
                    IDs[i] = fname
    bamfiles = IDs
    assert len(set(bamfiles)) == len(bamfiles), (f'\n{init_IDs}\n{IDs}\nsample names\
        should uniquely identify their corresponding bam file')

    # create multiple parallel bedtools calls
    num_processes = multiprocessing.cpu_count() - 1
    m_peaks = [bedfile for i in range(len(bamfiles))]
    count_df = pd.DataFrame()
    with multiprocessing.Pool(processes=num_processes) as pool:
        result = pool.starmap(get_coverage, zip(m_peaks,bamfiles))
        # maintain order of samples observed in metadata file
        sorted_result = []
        for id_ in IDs:
            for res in result:
                if id_ in res[0]:
                    print(f'peak_count_matrix(): adding counts for {id_}')
                    sorted_result.append(res[1])
        count_df = pd.DataFrame(sorted_result)

    count_df.insert(0, 'name', init_IDs)
    count_df = count_df.set_index('name')
    return count_df


def get_coverage(peak_file, bamfile):
    """
    Called by pool.map() to retrieve coverage at each peak position
    """
    names = []
    depths = []
    bed = proc_bed(peak_file)
    mbc_out = bed.multi_bam_coverage(bams=bamfile,f='0.0')
    for line in mbc_out:
        line = str(line).strip().split('\t')
        assert len(line) >= 7, '`pybedtools.multi_bam_coverage` failed'
        names.append(line[3])
        depths.append(float(line[6]))
    assert len(names) == bed.count()
    return bamfile,dict(zip(names,depths))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--peakfile', help="BED file containing peaks")
    parser.add_argument('-m', '--metadata', default=None,
                        help="sample metadata file. leftmost column should list sample\
                          identifiers for which there exists an alignment file named \
                          accordingly in the directory `--bamdir`")
    parser.add_argument('--bamdir', type=str, default='.')
    parser.add_argument('-o','--outfile', default='ROCCO_out_countmat.tsv')
    args = vars(parser.parse_args())
    result_df = peak_count_matrix(args['peakfile'], args['metadata'],
                      bamdir=args['bamdir'])

    result_df.to_csv(args['outfile'], sep='\t')
if __name__ == "__main__":
    main()
