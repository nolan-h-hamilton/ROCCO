import linecache
import os
import sys
import shlex
import subprocess
from pybedtools import BedTool

def wrap_run(cmd):
    print(cmd)
    try:
        cmd_list = shlex.split(cmd)
        proc = subprocess.run(cmd_list,
            shell=False,
            capture_output=True,
            text=True)
        return proc
    except Exception as e:
        print(f"Error executing command: {e}")
        return None

def line_k(fname, k):
    try:
        line = linecache.getline(fname, k)
        line = line.strip()
        return line
    except FileNotFoundError:
        return None

def is_nonempty(fname, min_bytes=10):
    return os.path.getsize(fname) >= min_bytes

def log(txt, script):
    print(f'{script}: {txt}')


wig_min_bytes = 100
min_wig_line_ct = 100
jac_min = .99

log('ensuring cwd is `tests` and parent directory is `ROCCO`...',
    sys.argv[0])
assert os.path.abspath(os.getcwd())[-5:] == 'tests'
assert os.path.abspath(os.path.pardir)[-5:] == 'ROCCO'

log('changing current working directory to parent dir.',
    sys.argv[0])
os.chdir(os.path.pardir)
assert os.path.abspath(os.getcwd())[-5:] == 'ROCCO'

log('running ROCCO_chrom.py for chr21, store output in `tests/output`', sys.argv[0])
wrap_run('rocco chrom --chrom chr21 --wig_path tests/data/tracks_chr21 -b .02 -N 0 --outdir tests/output')
log('running ROCCO_chrom.py for chr22, store output in `tests/output`', sys.argv[0])
wrap_run('rocco chrom --chrom chr22 --wig_path tests/data/tracks_chr22 -b .02 -N 0 --outdir tests/output')

assert os.path.exists('tests/output/ROCCO_out_chr21_0.02_1.0_0.0_1.0_1.0_1.0.bed')
chr21_bed = BedTool('tests/output/ROCCO_out_chr21_0.02_1.0_0.0_1.0_1.0_1.0.bed')
assert os.path.exists('tests/output/ROCCO_out_chr22_0.02_1.0_0.0_1.0_1.0_1.0.bed')
chr22_bed = BedTool('tests/output/ROCCO_out_chr22_0.02_1.0_0.0_1.0_1.0_1.0.bed')
merged_bed = chr21_bed.cat(chr22_bed)

log('running ROCCO_gwide.py using tests/test_params.csv, store output in `tests/output` and combine', sys.argv[0])
wrap_run('rocco gwide -p tests/test_params.csv --outdir tests/output/combined --combine tests/combined.bed -N 0 --multi 2')
assert os.path.exists('tests/combined.bed')
combined_bed = BedTool('tests/combined.bed')

log('Comparing output from the ROCCO_chrom.py jobs and ROCCO.py job to ensure similarity', sys.argv[0])
assert merged_bed.jaccard(combined_bed)['jaccard'] >= jac_min

log('all cases passed', sys.argv[0])
