"""
Generate a count matrix over a ROCCO-generated peak file. Included for convenience,
feel free to use an alternative external tool to generate the count matrix for the
ROCCO bed/peak file and samples' alignments.

Usage:
    `count_matrix.py [-h] [-i PEAKFILE] [-m METADATA] [--bamdir BAMDIR] [--sample_column SAMPLE_COLUMN] [-o OUTFILE] [--delimiter DELIMITER] [--samples_by_peaks]`

Arguments:
    peakfile (str): (`-i`) peak file (BED format)
    metadata (str): (`-m`) path to sample metadata file. Assumed
        structure is a TSV/CSV file with a row for each sample
        and names included in a column `--sample_column <default='sample'>`.
    sample_column: name of the column in the metadata file `[-m METADATA]` in which sample IDs are listed. Defaults to 'sample'.
    bamdir (str): path to directory containing samples' BAM alignments and their corresponding `.bai` index files
    outfile (str): filename of generated count matrix
    delimiter (str): delimiter used in the input sample metadata file and output count matrix file
    samples_by_peaks (bool): If invoked, count matrix rows are sample names, columns are peak names. By default, rows are genomic regions named as (`chr_start_end`) and columns are sample names
"""

import os
import argparse
import numpy as np
import pandas as pd
import multiprocessing

try:
    import pybedtools as pybed
except ImportError:
    print('Though not a core dependency of ROCCO, the package `pybedtools` is required for the count_matrix.py script. Install via `pip install pybedtools`.')
    raise

def proc_bed(input_bedfile):
    """
    Standardize input BED/peak file: `chr\tstart\tend\tchr_start_end`
    """
    if not os.path.exists(input_bedfile):
        raise FileNotFoundError(f"Input BED file not found: {input_bedfile}")

    output_bedfile = input_bedfile + '.4col.bed'
    bed_string = ""

    with open(input_bedfile, 'r') as bf:
        for i, line in enumerate(bf):
            line = line.strip().split('\t')
            name = '_'.join(line[:3])
            if len(line) < 3:
                raise SyntaxError(f'Peak/BED file is not formatted correctly: {line}')
            line = "\t".join(line[:3])  # get chrom, start, end
            line += f'\t{name}\n'
            bed_string += line
    with open(output_bedfile, 'w') as output:
        output.write(bed_string)
    return output_bedfile


def peak_count_matrix(bedfile, metadata, bamdir='.', sample_column='sample', delim='\t', strip_extensions_sample_names=True):
    """
    Wrapper for `pybedtools.multi_bam_coverage()` for a specified
    `bedfile` and samples' BAM files in `bamdir`. Generates a TSV
    file, that can be interpreted as a count matrix suitable for
    differential analyses.

    Args:
        bedfile (str): path to peak data file (bed formatted)
        metadata (str): path to sample metadata file. Assumed
            structure is a TSV file with a row for each sample
            and column `sample_column` specifying their IDs.
        bamdir (str): path to the directory containing samples'
            BAM files. Defaults to current working directory.
        sample_column (str): name of column containing sample IDs
        strip_extensions_sample_names (bool): whether to ensure removal of potential '.bam' file extensions
            that may exist in the *samples' names*, given in `sample_column`.
    """
    if not os.path.exists(bedfile):
        raise FileNotFoundError(f"peak file not found: {bedfile}")

    if not os.path.exists(metadata):
        raise FileNotFoundError(f"metadata file not found: {metadata}")

    if not os.path.exists(bamdir) or not os.path.isdir(bamdir):
        raise FileNotFoundError(f"bamdir not found: {bamdir}")
    
    bamfiles = []
    # read sample `metadata` as a pandas dataframe
    try:
        samp_df = pd.read_csv(metadata,
                                index_col=False,
                                sep=delim)
        IDs = [x for x in list(samp_df[sample_column])]

        if strip_extensions_sample_names:
            IDs = [x.replace('.bam','')  for x in IDs]
        
        assert len(IDs) == len(set(IDs)), f'sample names are not unique:\n{IDs}\n{list(set(IDs))}\n'
    except (ValueError,FileNotFoundError, AssertionError):
        raise

    sample_bam_map = {}

    bamfiles = [os.path.join(os.path.normpath(bamdir), x) for x in os.listdir(bamdir) if os.path.splitext(x)[-1].lower() == '.bam']
    for ID in IDs:
        matches = []
        final_match = None
        for bamfile in bamfiles:
            if ID in bamfile:
                matches.append(bamfile)
        if len(matches) == 0:
            continue
        if len(matches) > 1:
            final_match = sorted(matches,key=len)[0]
        if len(matches) == 1:
            final_match = matches[0]
        if final_match is None:
            continue
        if final_match not in sample_bam_map:
            sample_bam_map.update({ID: final_match})
        
    bamfiles = list(sample_bam_map.values())
    # create multiple parallel bedtools calls
    num_processes = multiprocessing.cpu_count() - 1
    m_peaks = [bedfile for i in range(len(bamfiles))]
    count_df = pd.DataFrame()
    with multiprocessing.Pool(processes=num_processes) as pool:
        result = pool.starmap(get_coverage, zip(m_peaks,bamfiles))
        # maintain order of samples observed in metadata file
        sorted_result = []
        for ID in sample_bam_map.keys():
            for res in result:
                if res[0] == sample_bam_map[ID]:
                    sorted_result.append(res[1])
        count_df = pd.DataFrame(sorted_result)

    count_df.insert(0, sample_column, sample_bam_map.keys())
    return count_df

def get_coverage(peak_file, bamfile, min_overlap=1e-9):
    """
    Called by pool.map() to retrieve coverage at each peak position
    
    Args:
        peak_file (str): path of BED/peak file over which to compute coverage
        bamfile (str): path to bam file
        min_overlap (float): `-f` argument in bedtools multicov.
    """
    peaks = []
    counts = []
    bed = pybed.BedTool(proc_bed(peak_file))
    mbc_out = bed.multi_bam_coverage(bams=bamfile,f=str(min_overlap))
    for line in mbc_out:
        line = str(line).strip().split('\t')
        assert len(line) >= 4, f'\nunexpected multicov output: {bamfile},{bed},{line}'
        peak_name = line[3]
        overlap_count = line[4]
        peaks.append(peak_name)
        counts.append(float(overlap_count))
    return bamfile,dict(zip(peaks,counts))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--peakfile', help="BED file containing peaks")
    parser.add_argument('-m', '--metadata', default=None,
                        help="sample metadata file with sample IDs in `--sample_column`")
    parser.add_argument('--bamdir', type=str, default='.', help="path to directory containing samples' BAM files and their corresponding `.bai` index files")
    parser.add_argument('--sample_column', default='sample', type=str, help='metadata column containing unniquely identifying sample IDs')
    parser.add_argument('-o','--outfile', default='ROCCO_out_countmat.tsv')
    parser.add_argument('--delimiter', default='\t')
    parser.add_argument('--samples_by_peaks', default=False, action='store_true', help='If invoked, rows are sample names, columns are peak names. Otherwise rows are genomic regions (`chr_start_end`) and columns are sample names (default)')

    args = vars(parser.parse_args())


    count_mat = peak_count_matrix(args['peakfile'],
                      args['metadata'],
                      bamdir=args['bamdir'],
                      sample_column=args['sample_column'],
                      delim=args['delimiter'])
    counts = count_mat.to_numpy(dtype=int)[:,1:]
    sample_names =  list(count_mat[args['sample_column']])
    peak_names = [x for x in count_mat.columns[1:]]

    if args['samples_by_peaks']:
        pd.DataFrame(data=counts,
                    index=sample_names,
                    columns=peak_names).to_csv(args['outfile'], sep=args['delimiter'], index_label='sample')
    else:
        pd.DataFrame(data=counts.T,
                    index=peak_names,
                    columns=sample_names).to_csv(args['outfile'], sep=args['delimiter'], index_label='peak')

if __name__ == "__main__":
    main()
