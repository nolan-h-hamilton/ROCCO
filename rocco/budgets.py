"""
Compute chromosome-specific budgets from sample data based on relative read densities.

Procedure:
    ```
    For each chromosome:
        Randomly select a fraction, `--samp_rate`, of loci in the chromosome (i.e., randomly choose a set of indices/nucleotide positions)
        Compute the average signal value at each selected locus with respect to each sample/replicate's wig file.
            --> Denote the set of such values as `chrom_loci_avgs`
        Calculate the average signal value across all sampled loci, i.e., take the mean of values in `chrom_loci_avgs`
            --> Denote as `est_chrom_density`

    With `est_chrom_density` recorded for each chromosome,
        Scale the `est_chrom_density` values such that they average to `--smean`
    ```

Parameters:
    --wigdir (str): Specifies the path to the directory containing the tracks_chr[] subdirectories generated by `rocco prep`.
    --smean (float): The chromosome-specific budgets will scale to this value. Default is 0.035.
    --samp_rate (float): Fraction of loci to sample for computing the mean signal value.
    -s, --sizes (str): Chromosome sizes file. Default is 'hg38.sizes'.
    -o, --outfile (str): Name of the output file. Output CSV is structured for use with `rocco gwide --param_file`. Default is 'params.csv'.
    -L, --locsize (int): Step-size in wig files. Default is 50.

Input:
    Path to directory containing the `tracks_chr[]` output from `rocco prep`

Output:
    CSV file that can be supplied as the parameter file for `rocco gwide --param_file`

Example [from demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb):
    ```
    rocco budgets --wigdir . -s hg38.sizes --smean .035 --samp_rate .20 -o params.csv
    ```

"""

import os
import argparse
import random
import numpy as np
from typing import Dict
import pysam
from . import rocco_aux
from . import chrom

def sample_wigs(chrom_wigdir, samp_rate=1.0, locsize=50) -> float:
    """
    Sample loci in `chrom`, compute mean signal values at each, compute mean over all sampled loci

    Args:
        chrom_wigdir (str): The path to the chromosome-specific wig directory. e.g., `tracks_chr1`
        samp_rate (float): Fraction of loci to sample for computing the mean signal value.
        locsize (int): Step-size in wig files.

    Returns:
        float: mean signal value normalized by the sampling rate and the number of loci.

    """
    chrom_wigdir = rocco_aux.trim_path(chrom_wigdir)
    wig_files = [chrom_wigdir + '/' + fname
                 for fname in os.listdir(chrom_wigdir) if 'wig' in fname.split('.')[-1]]
    start_, end_ = chrom.get_start_end(chrom_wigdir)

    signal_matrix = chrom.collect_wigs(wig_files,
                                 start=start_,
                                 end=end_,
                                 locus_size=locsize)
    sum_=0
    for i, loc in enumerate(signal_matrix.columns):
        if random.random() <= samp_rate:
            sigs = np.array(signal_matrix[loc],dtype=np.float32)
            sum_ += np.mean(sigs)
    return sum_/(samp_rate*len(signal_matrix.columns))

def scale_densities(chrom_dict, smean=.035) -> dict:
    """
    Scale chromosome-specific densities to mean `smean`.

    Args:
        chrom_dict (Dict): A dictionary mapping chromosome names to density values.
        smean (float): The target mean value to which densities will be scaled.

    Returns:
        Dict: A dictionary with chromosome names as keys and budget estimates averaging to `smean`

    """
    curr_mean = np.mean([float(x) for x in chrom_dict.values()])
    for key in chrom_dict.keys():
        chrom_dict[key] *= smean/curr_mean
    return chrom_dict

def main(args):

    chrom_dict = rocco_aux.parse_size_file(args['sizes'])
    args['wigdir'] = rocco_aux.trim_path(args['wigdir'])
    for key in chrom_dict.keys():
        chrom_dict[key] = sample_wigs(args['wigdir'] + "/tracks_" + str(key), samp_rate=args['samp_rate'],locsize=args['locsize'])
    chrom_dict = scale_densities(chrom_dict,smean=args['smean'])
    out = open(args['outfile'],'w')
    out.write('chromosome,input_path,budget,gamma,tau,c1,c2,c3\n')
    for key, val in chrom_dict.items():
        out.write(f'{key},{args["wigdir"] + "/tracks_" + str(key)},{round(val, 3)},NULL,NULL,NULL,NULL,NULL\n')
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute chromosome-specific budget parameters based on read densities')
    parser.add_argument('--wigdir',type=str, help="specifies the path to the directory containing the tracks_chr[] subdirectories generated by `rocco prep`" )
    parser.add_argument('--smean', type=float, default=.035, help="the chromosome-specific budgets will scale to this value")
    parser.add_argument('--samp_rate', type=float, help="fraction of loci to sample for computing mean signal value")
    parser.add_argument('-s','--sizes', type=str, default='hg38.sizes', help='chromosome sizes file')
    parser.add_argument('-o', '--outfile', type=str, default='params.csv', help='name of output file. Output CSV is structured for use with `rocco gwide --param_file`')
    parser.add_argument('-L','--locsize', type=int, default=50, help='step-size in wig files')
    args = vars(parser.parse_args())
    main(args)