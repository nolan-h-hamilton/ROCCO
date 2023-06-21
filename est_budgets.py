"""
Estimates budgets for each chromosome based on avg. read density
across multiple BAM files in `--bamdir`. Output (stdout) is formatted for
easy use as a parameter CSV file for ROCCO.py (`--param_file`).

The results may serve as a good starting point for users wishing to apply chromosome-specific budgets.

Usage:
    python est_budgets.py [-d <bamdir>] [-s <sizes>] [-a <min_value>] [-b <max_value>] [--read_length <read_length>]

Arguments:
    bamdir: Path to the directory containing BAM files.
    sizes: Path to a chromosome sizes file
    a: minimum budget for any chromosome
    b: maximum budget for any chromosome
    read_length: read length, estimated by default

"""

from typing import Dict
import argparse
import os
import subprocess
import pysam
import rocco_aux

def est_rdlen(bamfile: str, chromosome: str, total_reads: int, sample_frac: float = .0001, min_reads: int = 50) -> float:
    """
    estimate mean read length by subsampling `bamfile`.

    Args:
        bamfile (str): filepath to BAM alignment
        chromosome (str): name of chromosome, e.g., 'chr1'
        total_reads (int): total number of reads in `bamfile`
        sample_frac (float): fraction of reads to randomly sample for the estimate. Defaults to .0001.
        min_reads (int): minimum number of reads to sample for the estimate. Defaults to 50.

    Returns:
        float: estimated mean length of reads in `bamfile` in `chromosome`
    """
    sample_frac = max(min_reads/total_reads, sample_frac)
    cmd = f"samtools view --subsample {sample_frac} {bamfile} {chromosome} | awk '{{print length($10)}}' | awk '{{ sum += $1; count++ }} END {{ avg = sum / count; print avg }}'"
    output = subprocess.check_output(cmd, shell=True, text=True)
    avg_rdlen = float(output.strip())
    return avg_rdlen


def rd_dens(bamfile: str, size_file: str = "hg38.sizes", a: float = 0.0, b: float = 0.05, rdlen: int = -1) -> Dict[str, float]:
    """
    Calculate read density for a single BAM file and apply min-max normalization
      to compute a budget for each chromosome.

    Args:
        bamfile (str): Path to the BAM file.
        size_file (str): path to a chromosome sizes file
        a (float): Minimum allowed budget for a particular chromosome
        b (float): Maximum allowed budget for a particular chromosome.
        rdlen (int): Read length. If unspecified this value is estimated
            with `est_budgets.est_rdlen()`

    Returns:
        dict: A dictionary of estimated budgets for each chromosome.
    """
    size_dict = rocco_aux.parse_size_file(size_file)
    aln = pysam.AlignmentFile(bamfile)
    chr_reads = {}

    for chrom_record in aln.get_index_statistics():
        if chrom_record[0] not in size_dict.keys():
            continue
        chr_reads.update({chrom_record[0]: chrom_record[3]})

    for key in chr_reads.keys():
        if rdlen > 0:
            chr_reads[key] *= rdlen
        else:
            chr_reads[key] *= est_rdlen(bamfile, key, chr_reads[key])
        chr_reads[key] /= size_dict[key]

    min_val = min(chr_reads.values())
    max_val = max(chr_reads.values())

    for key in chr_reads.keys():
        chr_reads[key] = (((chr_reads[key] - min_val) / (max_val - min_val)) * (b - a) + a)

    return chr_reads


def avg_rd(bamdir: str, size_file: str = "hg38.sizes", a: float = 0.0, b: float = 0.05, rdlen:int = -1) -> Dict[str, float]:
    """
    Calculate the average read density across multiple BAM files.

    Args:
        bamfile (str): Path to the BAM file.
        size_file (str): path to a chromosome sizes file
        a (float): Minimum allowed budget for a particular chromosome
        b (float): Maximum allowed budget for a particular chromosome.
        rdlen (int): Read length. If unspecified this value is estimated
            with `est_budgets.est_rdlen()`

    Returns:
        dict: Dictionary with estimated budgets for each chromosome.
    """
    bamdir = rocco_aux.trim_path(bamdir)
    bamfiles = [f'{bamdir}/{x}' for x in os.listdir(bamdir)
                if x.split('.')[-1] == 'bam']
    print(bamfiles)
    bam_rd_dicts = [rd_dens(x, size_file=size_file, a=a, b=b, rdlen=rdlen) for x in bamfiles]
    final_dict = dict.fromkeys(bam_rd_dicts[0].keys(), 0)

    for key in final_dict.keys():
        for bam_rd_dict in bam_rd_dicts:
            final_dict[key] += bam_rd_dict[key] / len(bam_rd_dicts)

    return final_dict


def main():
    parser = argparse.ArgumentParser(description='Estimate budgets')
    parser.add_argument('-d', '--bamdir', type=str, help='Path to the directory containing .bam and .bai files')
    parser.add_argument('-s', '--sizes', type=str, help='chromosome sizes file')
    parser.add_argument('-a', type=float, default=0.0, help='Minimum allowed budget')
    parser.add_argument('-b', type=float, default=0.05, help='Maximum allowed budget')
    parser.add_argument('--read_length', type=int, default=-1, help="if you know the average read length, specify it with this parameter")
    args = parser.parse_args()
    dict_output = avg_rd(args.bamdir, args.sizes,  a=args.a, b=args.b, rdlen=args.read_length)
    print('chromosome,input_path,budget,gamma,tau,c1,c2,c3')
    for key, val in dict_output.items():
        print(f'{key},tracks_{key},{round(val, 3)},NULL,NULL,NULL,NULL,NULL')


if __name__ == '__main__':
     main()
