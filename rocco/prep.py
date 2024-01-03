"""
# prep.py

Obtain ROCCO conformable input from samples' BAM files.

ROCCO accepts multiple samples' coverage tracks that are assumed to have been generated following comparable
experimental, alignment, QC, etc. protocols. Though ROCCO is generally robust to deviations from this assumption,
we divide the BAM->WIG step and ROCCO itself so that users can, if necessary, apply their own preprocessing
protocol to the alignment files for computing coverage, normalization, etc. as their experimental context
and goals dictate. This general utility for generating ROCCO conformable input (coverage tracks) from
BAM files is offered as part of the ROCCO package as a subcommand `prep` for convenience and has produced
sufficient results across various tested scenarios.

BAM files for each sample in `--bamdir` are used to create fixed-step bigwig signal tracks.
These signal tracks are then split by chromosome into wig files in the subdirectories `tracks_chr[]`.

The resulting `tracks_chr[]` directories can be listed in the `--param_file` parsed by [rocco gwide](https://nolan-h-hamilton.github.io/ROCCO/rocco/gwide.html)

Alternatively, if running [rocco chrom](https://nolan-h-hamilton.github.io/ROCCO/rocco/chrom.html) manually, e.g.,
    ```
    rocco chrom --chrom chr11 --wig_path tracks_chr11
    ```

Parameters:
    -i, --bamdir (str, default='.'):
        Path to the directory containing BAM files.
    -o, --outdir (str, default='.'):
        Output directory.
    -s, --sizes (str, default='hg38'):
        A path to a chromosome sizes file OR an assembly name.
    -L, --interval_length (int, default=50):
        Wiggle track fixed step size
    -c, --cores (int, default=1):
        PEPATAC's `bamSitesToWig.py` cores parameter. Altering this parameter to use >1 core may cause issues on Mac OS.
    --scale_factor_file (str,default=None):
        Path to a tab-separated file with two columns: samples' alignment filenames and the other for their respective scale factors as computed with an external method.
        For instance, an output file from [deepTools](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html) `multiBamSummary --scalingFactors` command.

Examples:
    - create input `tracks_chr[]` directories for the hg38 assembly, no sizes file available:
        ```
        rocco prep -s hg38 --bamdir [/path/to/bamfiles]
        ```
    - create input `tracks_chr[]` directories from BAM files using a locally available sizes file:
        ```
        rocco prep -s [/path/to/sizefile] --bamdir [path/to/bamfiles]
        ```
    - create input`tracks_chr[]` directories from hg19 BAM files in the current directory, no sizes file available:
        ```
        rocco prep -s hg19 --bamdir .
        ```

Note:
    External utilities such as [`deepTools bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) can also be used to
    generate track data from BAM files, e.g.,
        ```
        bamCoverage -b sample.bam --binSize 50 -o sample.bam.bw [...]
        ```
        Then, to structure input using UCSC's `bigWigToWig` :
        ```
        bigWigToWig sample.bam.bw tracks_chr[]/chr[]_sample.bam.bw.wig -chrom=chr[]
        ```
    Such utilities offer several additional features that users may find helpful for
    normalization, sample comparison, visualization/plotting, etc. but ROCCO has not
    been tested rigorously on the signals generated from these tools.

"""

import os
import argparse
import subprocess
import pandas as pd
import pysam
from . import rocco_aux
import tempfile
import warnings


def main(args):
    # process command line arguments
    cwd = os.getcwd()
    args['bamdir'] = os.path.abspath(args['bamdir'])
    args['outdir'] = os.path.abspath(args['outdir'])
    sf_dict = None
    if args['scale_factor_file'] is not None:
        aln_names = pd.read_csv(args['scale_factor_file'],sep='\t',header=None)[0]
        scale_factors = pd.read_csv(args['scale_factor_file'],sep='\t',header=None)[1]
        sf_dict = dict(zip(aln_names,scale_factors))

    if not os.path.exists(args['sizes']):
        # if args['sizes'] is not a locally available file,
        # check UCSC for the assembly's sizes file
        try:
            downloaded_size_file = rocco_aux.get_size_file(args['sizes'])
        except Exception as ex:
            print(f'prep: {args["sizes"]} could not be found locally as a file (e.g. `hg38.sizes`) or on UCSC server as an assembly name (e.g. `hg38`).')
            raise
        args['sizes'] = downloaded_size_file

    chroms = list(rocco_aux.parse_size_file(args['sizes']).keys())
    bigwig_files = []
    tf = tempfile.NamedTemporaryFile(mode='w', delete=False)

    for fname in os.listdir(args['bamdir']):
        fname_sf = 1.0
        try:
            if args['scale_factor_file'] is not None:
                fname_sf = float(sf_dict[fname])
        except:
            warnings.warn(
                f'prep: unable to parse scale factor for {fname} in {args["scale_factor_file"]}...setting to 1.0')
            fname_sf = 1.0

        fname = os.path.join(args['bamdir'], fname)
        fname = os.path.abspath(fname)

        # create an index for bam files if not present already
        if rocco_aux.is_alignment(fname):
            aln = pysam.AlignmentFile(fname)
            try:
                aln.check_index()
            except ValueError:
                print(f'prep: no index file available for {fname}, calling pysam.index()')
                pysam.index(fname)

            print('prep: {}'.format(aln.filename.decode()))
            bstw_path = os.path.join(os.path.dirname(__file__), "bamSitesToWig.py") 
            bstw_cmd = ['python3', bstw_path,
                        '-i', aln.filename.decode(),
                        '--scale', str(fname_sf),
                        '-c', args['sizes'],
                        '-w', args['outdir'] + '/' + aln.filename.decode().split('/')[-1] + '.bw',
                        '-r', str(args['interval_length']),
                        '-m', 'atac',
                        '-l', str(int(args['interval_length'])//2),
                        '-p', str(args['cores']),
                        '--variable-step',
                        '--limit', " ".join([x for x in chroms if rocco_aux.has_reads(bamfile=aln.filename.decode(),min_reads=1,chrom=x)])]

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

    # after running bamSitesToWig, split each bigwig file by chromosome into tracks_chr[]
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

            # if wig file is nonempty, place in tracks_chr[]
            if os.path.getsize(chrom + '_' + bigwig_file + '.wig') > 1024:
                subprocess.run(["mv", chrom + '_' + bigwig_file + '.wig',
                        "tracks_{}".format(chrom)], check=True)
            else:
                subprocess.run(["rm", chrom + '_' + bigwig_file + '.wig'], check=True)

    os.chdir(cwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BAM preprocessing')
    parser.add_argument('-i','--bamdir', default='.', type=str, help='Path to the directory containing BAM files.')
    parser.add_argument('-o', '--outdir', default= '.', help='Output directory')
    parser.add_argument('-s', '--sizes',
                        default='hg38',
                        help="A path to a chromosome sizes file. OR\
                        an assembly name")
    parser.add_argument('-L', '--interval_length',
                        default=50,
                        help="wiggle track fixed step size")
    parser.add_argument('-c', '--cores',
                        type=int,
                        default=1,
                        help="`bamSitesToWig.py`'s cores  parameter.\
                        Altering this parameter to use >1 core\
                        may cause issues on Mac OS")
    parser.add_argument('--scale_factor_file', default=None, help="Path to a tab-separated file with two columns:\
                        the first for samples' alignment filenames and the other for their respective scale factors.\
                        e.g., the output file from deepTools `multiBamSummary --scalingFactors` command.")
    parser.add_argument('--multi',
                         default=True)
    args = vars(parser.parse_args())
    main(args)
