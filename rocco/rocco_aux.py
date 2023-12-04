"""
# rocco_aux.py

Miscellaneous helper functions

#### [Project Homepage](https://github.com/nolan-h-hamilton/ROCCO/)
"""
import os
import subprocess
import pandas as pd
import numpy as np
import pysam
import tempfile
import warnings
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor

def sort_combine_bed(outfile: str, dir_: str = '.', exclude_list: list = ['EBV', 'M', 'MT']):
    """
    Sorts and combines chromosome-specific peak files. Creates a new
    bed file `outfile` to store the combined results.

    Parameters:
        outfile (str): output file name
        dir_ (str): directory containing chromosome-specific bed files
            Default is the current directory ('.').
        exclude_list (list): list of chromosomes to exclude. (default: ['EBV', 'M', 'MT'])

    Notes:
        - this function assumes the bed files are named with the\
            same convention and delimiter '_'
    """
    dir_ = os.path.normpath(dir_)

    filenames = [f_ for f_ in os.listdir(dir_)
                 if f_.split('.')[-1] == 'bed' and 'ROCCO_out' in f_
                 and f_ not in exclude_list]

    name_template = []
    try:
        name_template = filenames[0].split('_')
    except (IndexError, ValueError) as none_found:
        print('sort_combine_bed(): no valid ROCCO output files found in `dir_`')
        raise none_found
    chr_index = [i for i in range(len(name_template)) if 'chr' in name_template[i][:3]][0]
    filenames = sorted(filenames,
                       key=lambda x: int(val) if (val := x.split('_')[chr_index][3:]).isnumeric() else ord(val))
    filenames = [os.path.join(dir_, fname) for fname in filenames]
    with open(outfile, mode='w', encoding='utf-8') as outfile_:
        for fname in filenames:
            cat_process = subprocess.Popen(('cat', fname), stdout=outfile_.fileno())
            cat_process.wait()


def download_file(url, output_path):
    "download file at `url` to `output_path` using either curl or wget"
    try:
        subprocess.run(['curl', '-o', output_path, url], check=True)
    except subprocess.CalledProcessError:
        try:
            subprocess.run(['wget', '-O', output_path, url], check=True)
        except subprocess.CalledProcessError:
            raise RuntimeError("Failed to download the file.")


def get_size_file(assembly='hg38', exclude_list=['EBV', 'M', 'MT']) -> str:
    """
    if `assembly` is included in the UC Santa Cruz genome repository
    this function downloads, filters, and sorts the corresponding sizes
    file to `assembly`.sizes. Chromosomes in `exclude_list` are ignored.

    Args:
        assembly (str, optional): name of genome assembly. Defaults to 'hg38'.
        exclude_list (list, optional): list of chromosome names to exclude in sizes file.
            Defaults to ['EBV', 'M', 'MT'].

    Returns:
        str: path to sizes file, `assembly`.sizes
    """
    assembly = assembly.lower()
    url = f"http://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chrom.sizes"
    output_path = f"{assembly}.sizes"
    print(f'\n{url} --> {output_path}\n')
    download_file(url, output_path)
    with open(output_path, 'r') as file_:
        lines = file_.readlines()
    filtered_lines = [line for line in lines if line.strip() and not any(excluded_chr in line for excluded_chr in exclude_list) and '_' not in line]
    sorted_lines = sorted(filtered_lines, key=lambda x: int(x.split('\t')[0][3:]) if x.split('\t')[0][3:].isnumeric() else ord(x.split('\t')[0][3:]))
    sorted_output_path = f"{assembly}.sizes"
    with open(sorted_output_path, 'w') as file_:
        file_.writelines(sorted_lines)
    return sorted_output_path


def parse_size_file(size_file, exclude_list=['EBV', 'M', 'MT']) -> dict:
    """Parse a size file and return a dictionary {chr: size}.

    Args:
        size_file (str): Path to the size file.
        exclude_list (list): list of chromosomes to exclude. (default: ['EBV', 'M', 'MT'])
    Returns:
        dict: a dictionary where chromosome names are keys and their sizes are values.
    """
    try:
        df = pd.read_csv(size_file, sep='\t', header=None)
    except Exception as ex:
        print(f'\nrocco_aux.parse_size_file(): Could not parse {size_file}. Ensure it is a valid sizes file.\n')
        raise
    size_dict = {key: val
                 for key, val in OrderedDict(zip(df.iloc[:, 0], df.iloc[:, 1])).items()
                    if key not in exclude_list}
    return size_dict


def is_alignment(filepath) -> bool:
    """
    Check if a file is an alignment file.

    Args:
        filepath (str): The path to the file.

    Returns:
        bool: True if the file is an alignment file, False otherwise.
    """
    if filepath.split('.')[-1] in ['bam','sam']:
        try:
            pysam.head("-n 1", filepath)
            return True
        except pysam.SamtoolsError:
            return False
    return False

def run_par(cmd_file, threads=-1, verbose=False):
    """
    Runs shell commands in `cmd_file` in parallel

    Args:
        cmd_file (str): file containing newline-separated commands
        threads (int): number of threads to use. if default '-1',
            then the number of threads is set to 1 + os.cpu_count()//2
        verbose (bool): whether to print subprocess output to stdout
    """
    assert os.path.exists(cmd_file), f'supplied `cmd_file` could not be found'
    if threads == -1:
        threads = 1 + os.cpu_count()//2
    with open(cmd_file, 'r') as file:
        commands = file.readlines()

    commands = [cmd.strip() for cmd in commands][::-1]


    def run_command(command, verbose=verbose):
        if not verbose:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            process.communicate()
            return command, process.returncode
        if verbose:
            process = subprocess.Popen(command, shell=True, stderr=subprocess.DEVNULL)
            process.communicate()
            return command, process.returncode

    with ThreadPoolExecutor(max_workers = threads) as executor:
        futures = [executor.submit(run_command, cmd) for cmd in commands]
        results = [future.result() for future in futures]
    if verbose:
        for command, returncode in results[::-1]:
            print(f"cmd: {command}\nretval: {returncode}\n")

def has_reads(bamfile, min_reads=1, chrom='chrY') -> bool:
    """
    Check if a BAM file contains a minimum number of reads for a specific chromosome.

    Args:
        bamfile (str): path to the BAM file
        min_reads (int): minimum read threshold
        chrom (str, optional): chromosome to check for reads

    Returns:
        bool: True if the BAM file contains at least 'min_reads' for 'chrom', False otherwise.
    """
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    has_min_reads = False

    for elem in bamfile.get_index_statistics():
        if elem[0] == chrom and elem[1] >= min_reads:
            has_min_reads = True
            break
    bamfile.close()
    return has_min_reads

def get_locus_size(wig_path, fixedStep=False) -> int:
    """
    Infers step size from replicates` wig files.

    Args:
        wig_path (str) : path to *directory* containing wig files
        fixedStep (bool): whether the wiggle tracks are in fixedStep format. Defaults to `False`.

    Returns:
        int: an integer locus/step size used in the wig files
    """
    if fixedStep:
        return parse_fixedStep_header([os.path.join(wig_path,x) for x in os.listdir(wig_path) if os.path.splitext(x)[-1] == '.wig'][0])[2]
    i = 0
    pos1 = 0
    pos2 = 0
    wig_file = [os.path.join(wig_path,x) for x in os.listdir(wig_path) if os.path.splitext(x)[-1] == '.wig'][0]
    with open(wig_file, mode='r', encoding='utf-8') as wig:
        for line in wig:
            line = line.strip()
            if line.replace('\t', '').replace(' ', '').isnumeric():
                if i == 0:
                    pos1 = int(line.split('\t')[0])
                if i == 1:
                    pos2 = int(line.split('\t')[0])
                if pos2 - pos1 > 0:
                    break
                i += 1
    return pos2 - pos1

def get_start_end(wig_path, fixedStep=False):
    """
    Infers starting and ending nucleotide position common to the
    multiple signal tracks

    To infer a reasonable genomic region over which to run ROCCO,
    this function returns the median starting point and endpoints
    among the multiple samples' signal files in `wig_path`.

    Args:
        wig_path (str): path to *directory* containing wig files
        fixedStep (bool): whether the wiggle tracks are in fixedStep format. Defaults to `False`.

    Returns:
        start (int): median starting position
        end (int): median ending position
    """
    starts = []
    ends = []
    def get_start(track, fixedStep=fixedStep):
        if fixedStep:
            start = parse_fixedStep_header(track)[1]
        with open(track, mode='r', encoding='utf-8') as wig:
            for line in wig:
                line = line.strip()

                if line.replace('\t', '').replace(' ', '').isnumeric():
                    start = int(line.split('\t')[0])
                    return start
        return None

    def get_end(track, fixedStep=fixedStep):
        if fixedStep:
            chrom,start,step = parse_fixedStep_header(track)
            length = num_nonempty_lines(track) - 1
            return start + step*(length-1)

        with subprocess.Popen(["tail", "-n1", track],
                              stdout=subprocess.PIPE) as proc:
            for line in proc.stdout:
                line = line.decode()
                line = line.strip()
                return int(line.split('\t')[0])
        return None

    wig_files = [os.path.join(wig_path,x) for x in os.listdir(wig_path) if os.path.splitext(x)[-1].lower() == '.wig']
    for wig_file in wig_files:
        starts.append(get_start(wig_file))
        ends.append(get_end(wig_file))
    return sorted(starts)[len(starts) // 2], sorted(ends)[len(ends) // 2]


def num_nonempty_lines(filepath):
    with open(filepath, 'r') as file:
        nonempty_lines = sum(1 for line in file if not line.strip().isspace())
    return nonempty_lines


def parse_fixedStep_header(fs_wig):
    with open(fs_wig, 'r') as file:
        i=0
        for line in file:
            if line.startswith("fixedStep"):
                line = line.strip()
                declaration_line = line.split()
                chrom = declaration_line[1].split('=')[1]
                start = int(declaration_line[2].split('=')[1])
                step = int(declaration_line[3].split('=')[1])
                if i > 0:
                    warnings.warn(f'fixedStep header line is not first line of {fs_wig}')
                return chrom,start,step
            i += 1
        return None

def tmp_fixedStep(fs_wig, delim='\t'):
    loci = []
    signal = []
    chrom_ct = 0
    chrom = None
    with open(fs_wig, 'r') as file:
        for line in file:
            if line.startswith("fixedStep"):
                if chrom_ct > 0:
                    warnings.warn(f'only the first chromosome, {chrom}, was extracted from {fs_wig}')
                    break
                line = line.strip()
                declaration_line = line.split()
                chrom = declaration_line[1].split('=')[1]
                start = int(declaration_line[2].split('=')[1])
                step = int(declaration_line[3].split('=')[1])
                chrom_ct += 1
            elif line.strip():
                loci.append(start)
                signal.append(float(line))
                start += step
    tmpname = f'{fs_wig}.vstep.wig'
    with open(tmpname, 'w') as vstep_file:
        for loc,sig in zip(loci,signal):
            vstep_file.write(f'{loc}{delim}{sig}\n')
    return tmpname


def read_wig(chrom_wig_file, start=0, end=10**10, locus_size=50, delim='\t', fixedStep=False):
    r"""
    Processes a chromosome-specific wig file for inclusion in the signal matrix $\mathbf{S}_{chr}$.
    Assumes fixed step-size between nucleotide positions.

    Args:
        wig_file (str): path to wiggle-formatted signal track
        start (int): inferred or manually-specified starting nucleotide
          position
        end (int): inferred or manually-specified ending nucleotide
          position
        locus_size (int): interval/locus size. Each wig file will have
          signal values at every `start + i*locus_size` nucleotide pos.
          for i = 0,1,2,...
        delim (str): delimiter used to separate nucleotide positions from signal values in
            wiggle tracks. Defaults to tabs `\t`
        fixedStep (bool): whether the wiggle track is in fixedStep format. Defaults to `False`.

    Returns:
        loci: list of starting nucleotide positions for each locus
        signal: enrichment signal value at each locus in `loci`
    """

    print("reading: {}".format(chrom_wig_file))
    tmpfile = None
    if fixedStep:
         tmpfile = tmp_fixedStep(chrom_wig_file, delim=delim)
         chrom_wig_file = tmpfile
    loci = []
    signal = []
    with open(chrom_wig_file, mode='r', encoding='utf-8') as wig:
        for line in wig:
            line = line.strip()
            try:
                line = line.split(delim)
                line[0] = int(line[0])
                line[1] = float(line[1])

                if line[0] < start:
                    continue
                if line[0] > end:
                    break

                # case: negative signal value
                # behavior: set signal value to zero
                if line[1] < 0:
                    line[1] = 0

                # case: wig begins after specified `start`
                # behavior: pad beginning with zero signal values
                if start < line[0] and len(loci) == 0:
                    for loc in range(start, line[0], locus_size):
                        loci.append(loc)
                        signal.append(0)
                    loci.append(line[0])
                    signal.append(line[1])
                    continue



                # case: missing data, gap between two loci
                # behavior: fill in gaps with zero signal vals
                if len(loci) > 0 and line[0] - loci[-1] > locus_size:
                    loci.append(loci[-1] + locus_size)
                    signal.append(0)
                    continue

                loci.append(line[0])
                signal.append(line[1])

            except ValueError:
                continue
            except IndexError:
                continue
    # case: wig ends prematurely
    # behavior: pad with zeros
    if loci[-1] < end:
        for loc in range(loci[-1] + locus_size, end + locus_size, locus_size):
            loci.append(loc)
            signal.append(0)
    if fixedStep and tmpfile is not None:
        try:
            os.remove(tmpfile)
        except:
            warnings.warn(f'temporary file {tmpfile} not deleted')
    return loci, signal

