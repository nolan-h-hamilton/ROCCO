"""
Obtain ROCCO conformable input from BAM files
"""

import os
import time
import os.path
import argparse
import pysam
import subprocess
from shutil import which
from subprocess import Popen, run, CalledProcessError, DEVNULL, check_call

"""
chromosome data for hg38 and mm39 included for convenience,
but different assemblies can be used by specifying a sizes
file with `--sizes`
"""
hg38_sizes = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415}

mm39_sizes = {
    'chr1': 195154279,
    'chr2': 181755017,
    'chr3': 159745316,
    'chr4': 156860686,
    'chr5': 151758149,
    'chr6': 149588044,
    'chr7': 144995196,
    'chr8': 130127694,
    'chr9': 124359700,
    'chr10': 130530862,
    'chr11': 121973369,
    'chr12': 120092757,
    'chr13': 120883175,
    'chr14': 125139656,
    'chr15': 104073951,
    'chr16': 98008968,
    'chr17': 95294699,
    'chr18': 90720763,
    'chr19': 61420004,
    'chrX': 169476592,
    'chrY': 91455967}

def tmp_sizes_file(assembly):
    if assembly == 'hg38':
        with open ('hg38.sizes', "w") as f:
            for key, value in hg38_sizes.items():
                f.write("{}\t{}\n".format(key, value))
            return 'hg38.sizes'
    elif assembly == 'mm39':
        with open ('mm39.sizes', "w") as f:
            for key, value in mm39_sizes.items():
                f.write("{}\t{}\n".format(key, value))
            return 'mm39.sizes'
    return assembly


def get_chroms(sizes):
    """Gathers a list of chromosomes present in the given sizes file

    Args:
        sizes (str) : a chromosome sizes filepath
        
    Returns:
        a list of chromosome names.

    Raises:
        FileNotFoundError: If `sizes` does not exist in the cwd.
    """
    
    chroms = []
    
    for line in open(sizes,'r'):
        line = line.strip()
        line = line.split('\t')
        chroms.append(line[0])
    return chroms


def create_chrom_dirs(sizes, names=None):
    """Create directories for each chromosome's signal files

    Args:
        sizes (str) : a chromosome sizes filepath.
        names (list): allows users to manually specify chromosome names but
          default behavior is to infer chromosome names from `sizes`.
    """
    # names parameter allows user to specify chromosome names manually,
    # but the default behavior is to infer them from the sizes file.
    if names is None:
        chroms = get_chroms(sizes)
        for chrom in chroms:
            chrom = chrom.replace('chromosome','chr')
            chrom = chrom.replace('chrom','chr')
            if not os.path.exists('tracks_' + chrom):
                print("creating directory `tracks_{}`".format(chrom))
                run(["mkdir", "tracks_{}".format(chrom)])

    else:
        chroms = names
        for chrom in chroms:
            if not os.path.exists('tracks_' + chrom):
                print("creating directory `tracks_{}`".format(chrom))
                run(["mkdir", "tracks_{}".format(chrom)])
        
            
def divide_by_chrom(bigwig_file, sizes,names=None):
    """Split genome-wide bigwig file into chromosome-specific wig files

    Divides larger bigwig by chromosome and converts to wiggle format.
    The resulting chromosome-specific wiggle files are moved to their
    appropriate `tracks_` directory.

    Args:
        bigwig_file (str) : path to a bigwig file.
        sizes (str) : a chromosome sizes filepath.
        names (list): allows users to manually specify chromosome names but
          default behavior is to infer chromosome names from `sizes`.
    """

    if names is None:
        chroms = get_chroms(sizes)
    else:
        chroms = names
        
    for chrom in chroms:
        chrom = chrom.replace('chromosome','chr')
        chrom = chrom.replace('chrom','chr')
        
        cmd_success = run(['bigWigToWig', bigwig_file,
                        chrom + '_' + bigwig_file +'.wig',
                        '-chrom={}'.format(chrom)]).returncode
        
        run(["mv", chrom + '_' + bigwig_file +'.wig',
                        "tracks_{}".format(chrom)])

        
def pepatac_script(scriptname, dload='wget'):
    # download necessary PEPATAC tools using wget
    if not os.path.exists(scriptname):
        print("\nCannot find {} in working directory...downloading\
        via wget".format(scriptname))
        # link to PEPATAC repo
        cmd = "https://raw.githubusercontent.com/databio/pepatac/master/tools/{}".format(scriptname)
        download_success = run([dload,cmd]).returncode
        if download_success != 0:
            raise Exception("Lacking dependencies...unable to download \
              via wget")
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--bamdir',
                        default='.',
                        help="path to directory containing BAM files")
    parser.add_argument('-s','--sizes',
                        default='hg38',
                        help="A path to a chromosome sizes file. OR\
                        one of the assembly names `hg38`, `mm39`.\
                        NOTE: if no file or assembly is specified, the\
                        hg38 assembly will be used by default.")
    parser.add_argument('-L','--interval_length',
                        default=50,
                        help="wiggle track fixed interval length")
    parser.add_argument('-c', '--cores',
                        type=int,
                        default=1,
                        help="Number of cores to process BAM files with.\
                        Altering this parameter to use >1 core can speed \
                        things up considerably but may cause issues on\
                        Mac OS")
    parser.add_argument('--index',
                        action="store_true",
                        default=False,
                        help="invoke to create index files for each BAM")
    parser.add_argument('--retain_links', action="store_true",
                        default=False,
                        help="invoke to preserve hard-linked files")
    
    args = vars(parser.parse_args())
    print(args)

    if args['bamdir'][-1] == '/':
        args['bamdir'] = args['bamdir'][0:-1]
        
    args['sizes'] = tmp_sizes_file(args['sizes'])

    pepatac_script('bamSitesToWig.py')
    pepatac_script('cutsToWig.pl')
    pepatac_script('smoothWig.pl')

    # create links to bam/index files if they are not in cwd
    initial_files = os.listdir()
    links = []
    if args['bamdir'] != '.':
        for fname in os.listdir(args['bamdir']):
            if fname.split('.')[-1] == 'bam':
                print('creating link for alignment {}/{}'.format(
                    args['bamdir'],fname))
                try:
                    run("ln -f {}/{} .".format(args['bamdir'], fname),
                        check=True,shell=True)
                    links.append(fname)
                except CalledProcessError as ex:
                    raise ex                
            if fname.split('.')[-1] == 'bai':
                print('creating link for index {}/{}'.format(
                    args['bamdir'],fname))
                try:
                    run("ln -f {}/{} .".format(args['bamdir'], fname),
                        check=True,shell=True)
                    links.append(fname)
                except CalledProcessError as ex:
                    raise ex
    
    # create directory for each chromosome
    create_chrom_dirs(args['sizes'])
    
    bamfiles = []
    # index/sort bam files
    print('running samtools and PEPATAC tools with {} cores'.format(
        args['cores']))
    for bamfile in os.listdir():
        if bamfile.split('.')[-1] == 'bam':
            print("\n\nprocessing: {}".format(bamfile))
            bamfiles.append(bamfile)
            print("reading alignment file with pysam...")
            sam_parsed = pysam.AlignmentFile(bamfile)

            # sort if needed
            if "SO:coordinate" not in str(sam_parsed.header):
                 print("sorting {} by coordinate...".format(bamfile))
                 pysam.sort(bamfile,'-o',bamfile,'-@',str(args['cores']))

            # index if the `--index` parameter was invoked
            if args['index'] == True:
                print('indexing {}...'.format(bamfile))
                pysam.index(bamfile)

    for bf in bamfiles:
        print('{}: running bamSitesToWig.py'.format(bf))
        run(['python3','bamSitesToWig.py','-i',
                          bf,'-c',
                          args['sizes'], '-w',
                          bf + '.bw',
                          '-r', str(args['interval_length']),
                          '-m', 'atac', '-p',
                          str(args['cores']),'--variable-step'],
                          shell=False, stderr=DEVNULL)
        print("Splitting bigWig file by chromosome")   
        divide_by_chrom(bf + '.bw', args['sizes'])
        
    if not args['retain_links']:
        for linked_file in links:
            if linked_file not in initial_files:
                print('deleting linked file: {}'.format(linked_file))
                run(['rm', linked_file])
            
if __name__ == "__main__":
    main()
