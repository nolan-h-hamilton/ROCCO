"""
Obtain ROCCO conformable input from BAM files
"""

import os
import argparse
import subprocess
import pysam
import pybedtools


def _is_link(fname):
    if os.path.islink(fname):
        return True
    if os.stat(fname).st_nlink > 1:
        return True
    return False


def sizes_file(assembly='hg38'):
    if not os.path.exists(assembly):
        fname = assembly + '.sizes'
        assembly_dict = pybedtools.chromsizes(assembly)
        keys = [x for x in assembly_dict.keys()
                if '_' not in x
                and x[3:] not in ['M', 'MT']]
        keys = sorted(keys, key=lambda x: int(val)
                      if (val := x[3:]).isnumeric() else ord(val))
        with open(fname, "w", encoding='utf-8') as assembly_file:
            for name in keys:
                assembly_file.write(
                    "{}\t{}\n".format(
                        name, assembly_dict[name][1]))
            return fname
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

    for line in open(sizes, 'r', encoding='utf-8'):
        line = line.strip()
        line = line.split('\t')
        chroms.append(line[0])
    return chroms


def create_chrom_dirs(sizes, names=None):
    """Create directories for each chromosome's signal files

    Args:
        sizes (str) : a chromosome sizes filepath.
    """
    chroms = get_chroms(sizes)
    for chrom in chroms:
        if not os.path.isdir('tracks_' + chrom):
            print("creating directory `tracks_{}`".format(chrom))
            os.mkdir("tracks_{}".format(chrom))


def divide_by_chrom(bigwig_file, sizes, names=None):
    """Split genome-wide bigwig file into chromosome-specific wig files

    Divides larger bigwig by chromosome and converts to wiggle format.
    The resulting chromosome-specific wiggle files are moved to their
    appropriate `tracks_` directory.

    Args:
        bigwig_file (str) : path to a bigwig file.
        sizes (str) : a chromosome sizes filepath.
    """

    chroms = get_chroms(sizes)
    for chrom in chroms:
        subprocess.run(['bigWigToWig', bigwig_file,
                        chrom + '_' + bigwig_file + '.wig',
                        '-chrom={}'.format(chrom)], check=True)

        subprocess.run(["mv", chrom + '_' + bigwig_file + '.wig',
                        "tracks_{}".format(chrom)], check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--bamdir',
                        default='.',
                        help="path to directory containing BAM files")
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
                        help="Number of cores to process BAM files with.\
                        Altering this parameter to use >1 core can speed \
                        things up considerably but may cause issues on\
                        Mac OS")
    parser.add_argument('--index',
                        action="store_true",
                        default=False,
                        help="invoke to create index files for each BAM")
    parser.add_argument('--retain', action="store_true",
                        default=False,
                        help="invoke to preserve created or hard-linked\
                            intermediate files")

    args = vars(parser.parse_args())

    if args['bamdir'][-1] == '/':
        args['bamdir'] = args['bamdir'][0:-1]

    args['sizes'] = sizes_file(args['sizes'])

    print(args)

    # create links to bam/index files if they are not in cwd
    initial_files = os.listdir()
    bamfiles = []
    links = []
    if args['bamdir'] != '.':
        for fname in os.listdir(args['bamdir']):
            if os.path.exists(fname) and not _is_link(fname):
                raise FileExistsError('{}: cannot overwrite'.format(fname))

            if fname.split('.')[-1] == 'bam':
                print('creating link for alignment {}/{}'.format(
                    args['bamdir'], fname))
                try:
                    subprocess.run("ln -f {}/{} .".format(args['bamdir'],
                                                          fname),
                                   check=True, shell=True)
                    bamfiles.append(fname)
                    links.append(fname)
                except subprocess.CalledProcessError as ex:
                    raise ex

            if fname.split('.')[-1] == 'bai':
                print('creating link for index {}/{}'.format(
                    args['bamdir'], fname))
                try:
                    subprocess.run("ln -f {}/{} .".format(args['bamdir'], fname),
                                   check=True, shell=True)
                    links.append(fname)
                except subprocess.CalledProcessError as ex:
                    raise ex
    else:
        bamfiles = [x for x in os.listdir() if x.split('.')[-1] == 'bam']

    # create directory for each chromosome
    create_chrom_dirs(args['sizes'])

    # index/sort bam files
    print('running with {} cores'.format(
        args['cores']))
    for bamfile in os.listdir():
        if bamfile.split('.')[-1] == 'bam':
            print("\n\nprocessing: {}".format(bamfile))

            print("reading alignment file with pysam...")
            sam_parsed = pysam.AlignmentFile(bamfile)

            # sort if needed
            if "SO:coordinate" not in str(sam_parsed.header):
                print("sorting {} by coordinate...".format(bamfile))
                pysam.sort(bamfile, '-o', bamfile, '-@', str(args['cores']))

            # index if the `--index` parameter was invoked
            if args['index'] == True:
                print('indexing {}...'.format(bamfile))
                pysam.index(bamfile)

    for bf in bamfiles:
        print('{}: running bamSitesToWig.py'.format(bf))
        subprocess.run(['python3', 'pepatac/bamSitesToWig.py', '-i',
                        bf, '-c',
                        args['sizes'], '-w',
                        bf + '.bw',
                        '-r', str(args['interval_length']),
                        '-m', 'atac', '-p',
                        str(args['cores']), '--variable-step'],
                       shell=False, check=True)
        print("Splitting bigWig file by chromosome")
        divide_by_chrom(bf + '.bw', args['sizes'])
        if not args['retain']:
            try:
                os.remove(bf + '.bw')
                os.remove(bf + '.bai')
            except FileNotFoundError as fnf_err:
                print(fnf_err)
                pass

    # remove links
    if not args['retain']:
        for linked_file in links:
            if linked_file not in initial_files and _is_link(linked_file):
                os.remove(linked_file)


if __name__ == "__main__":
    main()
