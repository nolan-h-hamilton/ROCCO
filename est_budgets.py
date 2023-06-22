"""
Estimates budgets for each chromosome based on average read density observed
across multiple BAM files in directory `--bamdir`. Output (stdout) is formatted
for easy use as a parameter CSV file for ROCCO.py (`--param_file`).

The results offer a fair starting point for users wishing to use chromosome-specific budgets--
some modification may be necessary depending on preferences.

Usage:
    est_budgets.py [-h] [-d BAMDIR] [-s SIZES] [-a A] [-b B] [--desired_avg DESIRED_AVG]

Arguments:
    bamdir (str): Path to the directory containing BAM files.
    size_file (str): path to a chromosome sizes file
    a (float): Minimum allowed budget for a particular chromosome
        this bound is ignored if `desired_avg` is non-negative.
    b (float): Maximum allowed budget for a particular chromosome
        this bound is ignored if `desired_avg` is non-negative.
    desired_avg (float): scaled read densities will average to this value
        if it is modified to be non-negative.
        
Examples:
    ```
    python3 est_budgets.py -d ../bamfiles -s hg38.sizes -a 0 -b .05
        chromosome,input_path,budget,gamma,tau,c1,c2,c3
        chr1,tracks_chr1,0.037,NULL,NULL,NULL,NULL,NULL
        chr2,tracks_chr2,0.038,NULL,NULL,NULL,NULL,NULL
        chr3,tracks_chr3,0.039,NULL,NULL,NULL,NULL,NULL
        chr4,tracks_chr4,0.036,NULL,NULL,NULL,NULL,NULL
        chr5,tracks_chr5,0.037,NULL,NULL,NULL,NULL,NULL
        chr6,tracks_chr6,0.039,NULL,NULL,NULL,NULL,NULL
        chr7,tracks_chr7,0.037,NULL,NULL,NULL,NULL,NULL
        chr8,tracks_chr8,0.038,NULL,NULL,NULL,NULL,NULL
        chr9,tracks_chr9,0.032,NULL,NULL,NULL,NULL,NULL
        chr10,tracks_chr10,0.038,NULL,NULL,NULL,NULL,NULL
        chr11,tracks_chr11,0.041,NULL,NULL,NULL,NULL,NULL
        chr12,tracks_chr12,0.04,NULL,NULL,NULL,NULL,NULL
        chr13,tracks_chr13,0.03,NULL,NULL,NULL,NULL,NULL
        chr14,tracks_chr14,0.032,NULL,NULL,NULL,NULL,NULL
        chr15,tracks_chr15,0.029,NULL,NULL,NULL,NULL,NULL
        chr16,tracks_chr16,0.036,NULL,NULL,NULL,NULL,NULL
        chr17,tracks_chr17,0.043,NULL,NULL,NULL,NULL,NULL
        chr18,tracks_chr18,0.035,NULL,NULL,NULL,NULL,NULL
        chr19,tracks_chr19,0.05,NULL,NULL,NULL,NULL,NULL
        chr20,tracks_chr20,0.04,NULL,NULL,NULL,NULL,NULL
        chr21,tracks_chr21,0.028,NULL,NULL,NULL,NULL,NULL
        chr22,tracks_chr22,0.031,NULL,NULL,NULL,NULL,NULL
        chrX,tracks_chrX,0.023,NULL,NULL,NULL,NULL,NULL
        chrY,tracks_chrY,0.0,NULL,NULL,NULL,NULL,NULL
    ```

    ```
    python3 est_budgets.py -d ../bamfiles -s hg38.sizes --desired_avg .035
        chromosome,input_path,budget,gamma,tau,c1,c2,c3
        chr1,tracks_chr1,0.037,NULL,NULL,NULL,NULL,NULL
        chr2,tracks_chr2,0.039,NULL,NULL,NULL,NULL,NULL
        chr3,tracks_chr3,0.039,NULL,NULL,NULL,NULL,NULL
        chr4,tracks_chr4,0.036,NULL,NULL,NULL,NULL,NULL
        chr5,tracks_chr5,0.037,NULL,NULL,NULL,NULL,NULL
        chr6,tracks_chr6,0.039,NULL,NULL,NULL,NULL,NULL
        chr7,tracks_chr7,0.038,NULL,NULL,NULL,NULL,NULL
        chr8,tracks_chr8,0.038,NULL,NULL,NULL,NULL,NULL
        chr9,tracks_chr9,0.033,NULL,NULL,NULL,NULL,NULL
        chr10,tracks_chr10,0.039,NULL,NULL,NULL,NULL,NULL
        chr11,tracks_chr11,0.041,NULL,NULL,NULL,NULL,NULL
        chr12,tracks_chr12,0.04,NULL,NULL,NULL,NULL,NULL
        chr13,tracks_chr13,0.031,NULL,NULL,NULL,NULL,NULL
        chr14,tracks_chr14,0.033,NULL,NULL,NULL,NULL,NULL
        chr15,tracks_chr15,0.03,NULL,NULL,NULL,NULL,NULL
        chr16,tracks_chr16,0.037,NULL,NULL,NULL,NULL,NULL
        chr17,tracks_chr17,0.043,NULL,NULL,NULL,NULL,NULL
        chr18,tracks_chr18,0.035,NULL,NULL,NULL,NULL,NULL
        chr19,tracks_chr19,0.05,NULL,NULL,NULL,NULL,NULL
        chr20,tracks_chr20,0.04,NULL,NULL,NULL,NULL,NULL
        chr21,tracks_chr21,0.029,NULL,NULL,NULL,NULL,NULL
        chr22,tracks_chr22,0.032,NULL,NULL,NULL,NULL,NULL
        chrX,tracks_chrX,0.024,NULL,NULL,NULL,NULL,NULL
        chrY,tracks_chrY,0.002,NULL,NULL,NULL,NULL,NULL
    ```
"""

import argparse
import numpy as np
import os
from typing import Dict
import pysam
import rocco_aux


def rd_dens(bamfile: str, size_file: str, a: float = 0.0, b: float = 0.05, desired_avg: float = -1.0) -> Dict[str, float]:
    """
    Calculate chromosome-specific budgets as scaled read density for `bamfile`

    Args:
        bamfile (str): Path to the BAM file.
        rdlen (int): Read length
        size_file (str): path to a chromosome sizes file
        a (float): Minimum allowed budget for a particular chromosome
            this bound is ignored if `desired_avg` is non-negative.
        b (float): Maximum allowed budget for a particular chromosome
            this bound is ignored if `desired_avg` is non-negative.
        desired_avg (float): scaled read densities will average to this value
            if it is modified to be non-negative.

    Returns:
        dict: A dictionary {chrom: budget}
    """
    size_dict = rocco_aux.parse_size_file(size_file)
    aln = pysam.AlignmentFile(bamfile)
    chr_reads = {}

    for chrom_record in aln.get_index_statistics():
        if chrom_record[0] not in size_dict.keys():
            continue
        chr_reads.update({chrom_record[0]: chrom_record[1]})

    for key in chr_reads.keys():
        chr_reads[key] /= size_dict[key]
    
    # scale mean to determine budgets if desired_avg >= 0
    current_avg = np.mean(list(chr_reads.values()))
    if desired_avg >= 0:
        for key in chr_reads.keys():
            chr_reads[key] *= desired_avg/current_avg
        return chr_reads
    
    # (default) apply min-max normalization to scale read densities in [a,b] if desired_avg < 0 
    min_val = min(chr_reads.values())
    max_val = max(chr_reads.values())
    for key in chr_reads.keys():
        chr_reads[key] = (((chr_reads[key] - min_val) / (max_val - min_val)) * (b - a) + a)
    return chr_reads


def avg_rd(bamdir: str, size_file: str, a: float = 0.0, b: float = 0.05, desired_avg: float = -1) -> Dict[str, float]:
    """
    Calculate the average read density across *multiple* BAM files in `bamdir`.

    Args:
        bamdir (str): Path to the directory containing BAM files.
        size_file (str): path to a chromosome sizes file
        a (float): Minimum allowed budget for a particular chromosome
            this bound is ignored if `desired_avg` is non-negative.
        b (float): Maximum allowed budget for a particular chromosome
            this bound is ignored if `desired_avg` is non-negative.
        desired_avg (float): scaled read densities will average to this value
            if it is modified to be non-negative.

    Returns:
        dict: A dictionary {chrom: budget}
    """
    bamdir = rocco_aux.trim_path(bamdir)
    bamfiles = [f'{bamdir}/{x}' for x in os.listdir(bamdir)
                if x.split('.')[-1] == 'bam']
    bam_rd_dicts = [rd_dens(x, size_file=size_file, a=a, b=b, desired_avg=desired_avg) for x in bamfiles]
    final_dict = dict.fromkeys(bam_rd_dicts[0].keys(), 0)

    for key in final_dict.keys():
        for bam_rd_dict in bam_rd_dicts:
            final_dict[key] += bam_rd_dict[key] / len(bam_rd_dicts)
    return final_dict


def main():
    parser = argparse.ArgumentParser(description='Compute chromosome-specific budget parameters based on observed read densities.\
        Uses min-max normalization (default) on the read density vals to yield budgets in interval [a,b] OR scales the values by a\
        constant such that their mean is `desired_avg`.')
    parser.add_argument('-d', '--bamdir', type=str, help='Path to the directory containing .bam and .bai files')
    parser.add_argument('-s', '--sizes', type=str, help='chromosome sizes file')
    parser.add_argument('-a', type=float, default=0.0, help='Minimum allowed budget. This bound is ignored if\
        `--desired_avg` is non-negative.')
    parser.add_argument('-b', type=float, default=0.05, help='Maximum allowed budget. This bound is ignored if\
        `--desired_avg` is non-negative.')
    parser.add_argument('--desired_avg', type=float, default=-1.0, help='Scaled read densities (i.e., budgets) will\
        average to this value if non-negative. Defaults to -1.')
    args = parser.parse_args()
    
    dict_output = avg_rd(args.bamdir, args.sizes,  a=args.a, b=args.b, desired_avg=args.desired_avg)
    print('chromosome,input_path,budget,gamma,tau,c1,c2,c3')
    for key, val in dict_output.items():
        print(f'{key},tracks_{key},{round(val, 3)},NULL,NULL,NULL,NULL,NULL')


if __name__ == '__main__':
     main()
