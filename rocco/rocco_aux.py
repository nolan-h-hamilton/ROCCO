import os
import subprocess
import pandas as pd
import pybedtools
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


def get_size_file(assembly='hg38', exclude_list=['EBV', 'M', 'MT']) -> str:
    """
    If `assembly` is the name of an assembly included in the pybedtools
    genome registry, this function will create a chromosome sizes file
    and return a path to it: `<assembly>.sizes`

    Args:
        assembly (str): genome assembly name
        exclude_list (list): list of chromosome names to exclude\
            Defaults to `['EBV', 'M', 'MT]`

    Returns:
        str: file path to a chromosome sizes file

    Notes: excludes any chromosome names with an underscore
    """
    fname = assembly + '.sizes'
    try:
        assembly_dict = pybedtools.chromsizes(assembly)
    except (AttributeError, OSError) as not_avail_ex:
        print(f'{assembly} is not available in the pybedtools genome registry')
        raise not_avail_ex
    keys = [x for x in assembly_dict.keys()
            if '_' not in x
            and x[3:] not in exclude_list]
    keys = sorted(keys, key=lambda x: int(val)
                    if (val := x[3:]).isnumeric() else ord(val))

    with open(fname, "w", encoding='utf-8') as assembly_file:
        for name in keys:
            assembly_file.write(
                "{}\t{}\n".format(
                    name, assembly_dict[name][1]))
    return fname


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

def run_par(cmd_file, verbose=False):
    """
    Runs shell commands in `cmd_file` in parallel

    Args:
        cmd_file (str): file containing newline-separated commands
        verbose (bool): whether to print subprocess output to stdout
    """
    assert os.path.exists(cmd_file), f'supplied `cmd_file` could not be found'

    with open(cmd_file, 'r') as file:
        commands = file.readlines()

    commands = [cmd.strip() for cmd in commands][::-1]
    print(f'rocco_aux.run_par: running following commmands in parallel\n{commands}\n')

    def run_command(command):
        if not verbose:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            process.communicate()
            return command, process.returncode
        if verbose:
            process = subprocess.Popen(command, shell=True, stderr=subprocess.DEVNULL)
            process.communicate()
            return command, process.returncode

    with ThreadPoolExecutor(max_workers = 1 + (os.cpu_count() // 2)) as executor:
        futures = [executor.submit(run_command, cmd) for cmd in commands]
        results = [future.result() for future in futures]
    for command, returncode in results[::-1]:
        print(f"cmd: {command}\nretval: {returncode}\n")
