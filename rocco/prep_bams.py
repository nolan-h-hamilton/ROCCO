"""
Obtain ROCCO conformable input from BAM files

This script is a ROCCO-specific wrapper for the PEPATAC (github.com/databio/pepatac)
tool `bamSitesToWig.py`. BAM files for each sample/replicate in `bamdir` are used to
create smooth, fixed-step bigwig signal tracks. These signal tracks are split by
chromosome, converted to human-readable .wig format, and then placed into directories
à la [this flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png).

The resulting `tracks_chr[]` directories can then be supplied to rocco chrom
via `--wig_path` to construct the signal matrix $\mathbf{S}_{chr}$.

```

Example:

    Run on toy data in `tests/data`:
    ```
    rocco prep --bamdir tests/data -s tests/data/test_sizes.sizes
    ```
Arguments:
    -i, --bamdir (str, default='.'):
        Path to the directory containing BAM files.

    -o, --outdir (str, default='.'):
        Output directory.

    -s, --sizes (str, default='hg38'):
        A path to a chromosome sizes file, or an assembly name.

    -L, --interval_length (int, default=50):
        Wiggle track fixed interval length.

    -c, --cores (int, default=1):
        The `bamSitesToWig.py` cores parameter. Altering this parameter to use >1 core may cause issues on Mac OS.

    --multi (bool, default=True):
        Set to `False` to run `bamSitesToWig` jobs sequentially. Note that this definition of `--multi` differs from `ROCCO_gwide.py`. 

    --index (int, default=1):
        Deprecated - backwards compatibility. Now, if BAM files are not indexed, pysam.index() is called by default. Before, this behavior was specified with this argument.

"""

import os
import argparse
import subprocess
import pysam
from . import rocco_aux
import tempfile


def main(args):
    # process command line arguments
    cwd = os.getcwd()
    args['bamdir'] = rocco_aux.trim_path(os.path.abspath(args['bamdir']))
    args['outdir'] = rocco_aux.trim_path(os.path.abspath(args['outdir']))

    if not os.path.exists(args['sizes']):
        args['sizes'] = rocco_aux.get_size_file(args['sizes'])

    bigwig_files = []
    tf = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for fname in os.listdir(args['bamdir']):
        fname = args['bamdir'] + '/' + fname
        fname = os.path.abspath(fname)
        if rocco_aux.is_alignment(fname):
            aln = pysam.AlignmentFile(fname)
            try:
                aln.check_index()
            except ValueError:
                print(f'no index file available for {fname}, calling pysam.index()')
                pysam.index(fname)
            print('{}: running bamSitesToWig.py'.format(aln.filename.decode()))
            bstw_path = os.path.join(os.path.dirname(__file__), "bamSitesToWig.py")
            bstw_cmd = ['python3', bstw_path,
                        '-i', aln.filename.decode(),
                        '-c', args['sizes'],
                        '-w', args['outdir'] + '/' + aln.filename.decode().split('/')[-1] + '.bw',
                        '-r', str(args['interval_length']),
                        '-m', 'atac',
                        '-p', str(args['cores']),
                        '--variable-step']
            if not args['multi']:
                # if `--multi False
                subprocess.run(bstw_cmd,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL,
                               shell=False, check=True)
            else:
                tf.write(' '.join(bstw_cmd))
                tf.write('\n')

            bigwig_file = args['outdir'] + '/' + aln.filename.decode().split('/')[-1] + '.bw'
            bigwig_files.append(bigwig_file)
    tf.close()
    if args['multi']:
        rocco_aux.run_par(tf.name)

    chroms = list(rocco_aux.parse_size_file(args['sizes']).keys())
    os.chdir(args['outdir'])
    for chrom in chroms:
        try:
            os.mkdir(f'tracks_{chrom}')
        except FileExistsError:
            pass
        for bigwig_file in [bw_file.split('/')[-1] for bw_file in bigwig_files]:
            subprocess.run(['bigWigToWig', bigwig_file,
                        chrom + '_' + bigwig_file + '.wig',
                        '-chrom={}'.format(chrom)], check=True)
            subprocess.run(["mv", chrom + '_' + bigwig_file + '.wig',
                        "tracks_{}".format(chrom)], check=True)
    os.chdir(cwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BAM preprocessing')
    parser.add_argument('-i','--bamdir', default='.', type=str, help='path to directory containing BAM files')
    parser.add_argument('-o', '--outdir', default= '.', help='output directory')
    parser.add_argument('-s', '--sizes',
                        default='hg38',
                        help="A path to a chromosome sizes file. OR\
                        an assembly name")
    parser.add_argument('-L', '--interval_length',
                        default=50,
                        help="wiggle track fixed interval length")
    parser.add_argument('-c', '--cores',
                        type=int,
                        default=1,
                        help="`bamSitesToWig.py`'s cores  parameter.\
                        Altering this parameter to use >1 core\
                        may cause issues on Mac OS")
    parser.add_argument('--multi',
                         default=True,
                         help='Set to `False` to run `bamSitesToWig` jobs\
                         sequentially.')
    parser.add_argument('--index',
                        type=int,
                        default=1,
                        help="Deprecated -- backwards compatibility")
    parser.add_argument('--bstw_path', default=None, help="Deprecated -- backwards compatibility")
    args = vars(parser.parse_args())
    main(args)
