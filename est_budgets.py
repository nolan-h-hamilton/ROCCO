"""
Estimate budgets for each chromosome using data from `--bamdir`

This script estimates budgets for each chromosome based on read density
across multiple BAM files in `--bamdir`. Output (stdout) is formatted for
easy use as a parameter CSV file for ROCCO.py.

Usage:
    python est_budgets.py [-d <bamdir>] [-s <sizes>] [-a <min_value>] [-b <max_value>]

Arguments:
    bamdir: Path to the directory containing BAM files.
    sizes: Path to a chromosome sizes file
    a: minimum budget for any chromosome
    b: maximum budget for any chromosome

Example:
```
$ python3 est_budgets.py -s hg38.sizes -d ../atac_bamfiles
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
"""

import os
import argparse
import pysam
from collections import OrderedDict
from typing import Dict


def parse_size_file(size_file: str) -> Dict[str, int]:
    """Parse a size file and return a dictionary {chr: size}.

    Args:
        size_file: Path to the size file.

    Returns:
        A dictionary where chromosome names are keys and their sizes are values.
    """
    chroms = []
    sizes = []

    with open(size_file, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip().split('\t')
            if '_' not in line[0]:
                chroms.append(line[0])
                sizes.append(int(line[1]))
    chroms = [x for x in chroms
            if '_' not in x
            and x[3:] not in ['M', 'MT']]
    return OrderedDict(zip(chroms, sizes))


def rd_dens(bamfile: str, size_file: str = "hg38.sizes", a: float = 0.01, b: float = 0.1) -> Dict[str, float]:
    """Calculate read density for a single BAM file and apply min-max normalization
      to compute a budget for each chromosome.

    Args:
        bamfile: Path to the BAM file.
        size_file: path to a chromosome sizes file
        a: Minimum allowed budget for a particular chromosome
        b: Maximum allowed budget for a particular chromosome.

    Returns:
        A dictionary of estimated budgets for each chromosome.
    """
    size_dict = parse_size_file(size_file)
    aln = pysam.AlignmentFile(bamfile)
    chr_reads = {}

    for chrom_record in aln.get_index_statistics():
        if chrom_record[0] not in size_dict.keys():
            continue
        chr_reads.update({chrom_record[0]: chrom_record[1]})

    for key in chr_reads.keys():
        chr_reads[key] /= size_dict[key]

    min_val = min(chr_reads.values())
    max_val = max(chr_reads.values())

    for key in chr_reads.keys():
        chr_reads[key] = (((chr_reads[key] - min_val) / (max_val - min_val)) * (b - a) + a)

    return chr_reads


def avg_rd(bamdir: str, size_file: str = "hg38.sizes", a: float = 0, b: float = 0.05) -> Dict[str, float]:
    """Calculate the average read density across multiple BAM files.

    Args:
        bamdir: Path to the directory containing BAM files.
        a: Minimum value for read density normalization.
        b: Maximum value for read density normalization.

    Returns:
        Dictionary with estimated budgets for each chromosome.
    """
    if bamdir[-1] == '/':
        bamdir = bamdir[0:-1]
    bamfiles = [f'{bamdir}/{x}' for x in os.listdir(bamdir)
                if x.split('.')[-1] == 'bam']
    bam_rd_dicts = [rd_dens(x, size_file=size_file, a=a, b=b) for x in bamfiles]
    final_dict = dict.fromkeys(bam_rd_dicts[0].keys(), 0)

    for key in final_dict.keys():
        for bam_rd_dict in bam_rd_dicts:
            final_dict[key] += bam_rd_dict[key] / len(bam_rd_dicts)

    return final_dict


def main():
    parser = argparse.ArgumentParser(description='Estimate budgets')
    parser.add_argument('-d', '--bamdir', type=str, help='Path to the directory containing .bam and .bai files')
    parser.add_argument('-s', '--sizes', type=str, help='chromosome sizes file')
    parser.add_argument('-a', type=float, default=0, help='Minimum allowed budget')
    parser.add_argument('-b', type=float, default=0.05, help='Maximum allowed budget')
    args = parser.parse_args()
    dict_output = avg_rd(args.bamdir, args.sizes,  a=args.a, b=args.b)
    print('chromosome,input_path,budget,gamma,tau,c1,c2,c3')
    for key, val in dict_output.items():
        print(f'{key},tracks_{key},{round(val, 3)},NULL,NULL,NULL,NULL,NULL')


if __name__ == '__main__':
     main()
