"""
Obtain ROCCO conformable input from BAM files

This script is a ROCCO-specific wrapper for the PEPATAC (github.com/databio/pepatac)
tool `bamSitesToWig.py`. BAM files for each sample/replicate in `bamdir` are used to
create smooth, fixed-step bigwig signal tracks. These signal tracks are split by
chromosome, converted to human-readable .wig format, and then placed into directories
à la [this flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png).

The resulting `tracks_chr[]` directories can then be supplied to ROCCO_chrom.py
via `--wig_path` to construct the signal matrix $\mathbf{S}_{chr}$.

```
usage: prep_bams.py [-h] [-i BAMDIR] [-o OUTDIR] [-s SIZES]
                    [-L INTERVAL_LENGTH] [-c CORES] [--multi] [--index INDEX]
                    [--bstw_path BSTW_PATH]
options:
  -h, --help            show this help message and exit
  -i BAMDIR, --bamdir BAMDIR
                        path to directory containing BAM files
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -s SIZES, --sizes SIZES
                        A path to a chromosome sizes file. OR an assembly name
  -L INTERVAL_LENGTH, --interval_length INTERVAL_LENGTH
                        wiggle track fixed interval length
  -c CORES, --cores CORES
                        `bamSitesToWig.py`'s cores parameter. Altering this
                        parameter to use >1 core may cause issues on Mac OS
  --multi               Invoke to run `bamSitesToWig` jobs simultaneously.
  --index INDEX         deprecated--included for backwards compatibility
  --bstw_path BSTW_PATH
                        path to bamSitesToWig.py script, included in
                        ROCCO/pepatac by default
```

Example:

    Run on toy data in `tests/data`:
    ```
    python3 prep_bams.py --bamdir tests/data -s tests/data/test_sizes.sizes --multi
    ```

Notes:
    - Several alternative tools exist for creating signal tracks in the desired format,\
        but have not yet been tested.
"""

import os
import argparse
import subprocess
import pysam
import rocco_aux
import tempfile


def main():
    parser = argparse.ArgumentParser(description='Description of your program.')
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
                         default=False,
                         action='store_true',
                         help='Invoke to run `bamSitesToWig` jobs\
                         simultaneously.')
    parser.add_argument('--index',
                        type=int,
                        default=1,
                        help="deprecated--included for backwards compatibility")
    parser.add_argument('--bstw_path', default='pepatac/bamSitesToWig.py', help="path to bamSitesToWig.py script, included in ROCCO/pepatac by default")
    args = vars(parser.parse_args())

    # process command line arguments
    cwd = os.getcwd()
    args['bamdir'] = rocco_aux.trim_path(os.path.abspath(args['bamdir']))
    args['outdir'] = rocco_aux.trim_path(os.path.abspath(args['outdir']))
    args['bstw_path'] = rocco_aux.trim_path(os.path.abspath(args['bstw_path']))
    if not args['multi']:
        print('prep_bams.py(): consider running with `--multi` for\
        parallel execution of bamSitesToWig.py jobs')

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
            bstw_cmd = ['python3', args['bstw_path'],
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
    main()
