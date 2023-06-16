import os
import subprocess
import pandas as pd
import pybedtools
from collections import OrderedDict

def trim_path(fname: str):
    """Remove trailing `/`s from directory paths"""
    while len(fname) > 1 and fname[-1] == '/':
        fname = fname[:-1]
    return fname

def sort_combine_bed(outfile: str, dir_: str = '.', exclude_: list = []) -> bool:
    """
    Sorts and combines chromosome-specific bed files. Creates a new
    bed file `outfile` to store the combined results.

    Parameters:
        outfile (str): The output file name where the sorted and combined bed file will be written.
        dir_ (str, optional): The directory containing the chromosome-specific bed files.
          Default is the current directory ('.').
        exclude_ (list, optional): list of filenames to exclude from result

    Returns:
        `True` if combined file created successfully, `False` otherwise

    Notes:
        - this function assumes the bed files are named with the\
            same convention and delimiter '_'
    """
    dir_ = trim_path(dir_)

    filenames = [f_ for f_ in os.listdir(dir_)
                 if f_.split('.')[-1] == 'bed' and 'ROCCO_out' in f_
                 and f_ not in exclude_]

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
            if cat_process.returncode > 0:
                return False
    return True

def get_size_file(assembly='hg38', exclude_list=['EBV', 'M', 'MT']) -> dict:
    """
    If `assembly` is the name of an assembly included in the pybedtools
    genome collection, this function will return a chromosome sizes file
    `<assembly>.sizes`

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

def parse_size_file(size_file, exclude_list=['EBV', 'M', 'MT']):
    """Parse a size file and return a dictionary {chr: size}.

    Args:
        size_file: Path to the size file.

    Returns:
        A dictionary where chromosome names are keys and their sizes are values.
    """
    df = pd.read_csv(size_file, sep='\t', header=None)
    size_dict = {key: val
                 for key, val in OrderedDict(zip(df.iloc[:, 0], df.iloc[:, 1])).items()
                    if key not in exclude_list}
    return size_dict