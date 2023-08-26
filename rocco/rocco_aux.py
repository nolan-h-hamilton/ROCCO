import os
import subprocess
import pandas as pd
import pysam
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor


def trim_path(fname: str) -> str:
    """
    Remove trailing `/`s from directory/file paths

    Args:
        fname (str): name

    Returns:
        str: trimmed name
    """
    while len(fname) > 1 and fname[-1] == '/':
        fname = fname[:-1]
    return fname


def sort_combine_bed(outfile: str, dir_: str = '.', exclude_list: list = ['EBV', 'M', 'MT']):
    """
    Sorts and combines chromosome-specific bed files. Creates a new
    bed file `outfile` to store the combined results.

    Parameters:
        outfile (str): The output file name where the sorted and combined bed file will be written.
        dir_ (str): The directory containing the chromosome-specific bed files.
            Default is the current directory ('.').
        exclude_list (list): list of chromosomes to exclude. (default: ['EBV', 'M', 'MT'])


    Notes:
        - this function assumes the bed files are named with the\
            same convention and delimiter '_'
    """
    dir_ = trim_path(dir_)

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
    filenames = [dir_ + '/' + fname for fname in filenames]
    with open(outfile, mode='w', encoding='utf-8') as outfile_:
        for fname in filenames:
            cat_process = subprocess.Popen(('cat', fname),
                                           stdout=outfile_.fileno())
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
    df = pd.read_csv(size_file, sep='\t', header=None)
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
    if 'bam' == filepath.split('.')[-1]:
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

    def run_command(command):
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
    for command, returncode in results[::-1]:
        print(f"cmd: {command}\nretval: {returncode}\n")
