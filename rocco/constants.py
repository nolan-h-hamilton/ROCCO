r"""
ROCCO: Constants
====================

Default effective genome sizes and chromosome size files for supported assemblies:
`hg38`, `hg19`, `mm10`, `mm39`, and `dm6`.

Universally optimal parameters/default values cannot be established for all datasets.
When possible, users may find improved results by tuning the CLI options to their specific data.


Variables
---------
GENOME_DICT : dict
    Dictionary of genome assemblies with effective genome size and chromosome sizes file.
    Keys are assembly names, values are dictionaries with keys 'effective_genome_size' and 'sizes_file'.
    Default assemblies: `hg38`, `hg19`, `mm10`, `mm39`, and `dm6`.

"""

import os
import warnings

def get_file_path(filename):
    file_path = os.path.join(os.path.dirname(__file__), filename)
    if not os.path.exists(file_path):
        warnings.warn(f"File not found: {file_path}")
    return file_path

GENOME_DICT = {
    'hg38': {'effective_genome_size': int(2.7e9), 'sizes_file': get_file_path('hg38.sizes')},
    'hg19': {'effective_genome_size': int(2.7e9), 'sizes_file': get_file_path('hg19.sizes')},
    'mm10': {'effective_genome_size': int(1.87e9), 'sizes_file': get_file_path('mm10.sizes')},
    'mm39': {'effective_genome_size': int(1.87e9), 'sizes_file': get_file_path('mm39.sizes')},
    'dm6': {'effective_genome_size': int(1.45e8), 'sizes_file': get_file_path('dm6.sizes')}}
