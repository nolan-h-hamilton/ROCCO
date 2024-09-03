r"""
ROCCO: Constants
====================
"""


import os
import warnings


def get_file_path(filename):
    file_path = os.path.join(os.path.dirname(__file__), filename)
    if not os.path.exists(file_path):
        warnings.warn(f"File not found: {file_path}")
    return file_path

GENOME_DICT = {
    'hg38': {'effective_genome_size': int(2.7e9), 'sizes_file': get_file_path('hg38.sizes'), 'params': get_file_path('hg_params.csv')},
    'hg19': {'effective_genome_size': int(2.7e9), 'sizes_file': get_file_path('hg19.sizes'), 'params': get_file_path('hg_params.csv')},
    'mm10': {'effective_genome_size': int(1.87e9), 'sizes_file': get_file_path('mm10.sizes'), 'params': get_file_path('mm_params.csv')},
    'mm39': {'effective_genome_size': int(1.87e9), 'sizes_file': get_file_path('mm39.sizes'), 'params': get_file_path('mm_params.csv')}
    }
                
