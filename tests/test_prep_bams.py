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

log('run prep_bams.py to prepare data for ROCCO algorithm...',
    sys.argv[0])
prepbam_proc = wrap_run(f'python3 prep_bams.py --bamdir tests/data --outdir tests/data -s tests/data/test_sizes.sizes')

log('ensure signal tracks have been created and organized properly',
    sys.argv[0])
assert len(os.listdir('tests/data/tracks_chr21')) >= 5
assert len(os.listdir('tests/data/tracks_chr21')) >= 5

for fname in os.listdir('tests/data/tracks_chr21'):
    if 'wig' not in fname.split('.')[-1]:
        continue
    fname = f'tests/data/tracks_chr21/{fname}'
    assert 'chr21' in fname
    assert is_nonempty(fname, wig_min_bytes)
    with open(fname, 'r') as f_:
        line_ct = len(f_.readlines())
    assert line_ct >= min_wig_line_ct
    assert int(line_k(fname, 3).split('\t')[0]) - int(line_k(fname, 2).split('\t')[0]) == 50

for fname in os.listdir('tests/data/tracks_chr22'):
    if 'wig' not in fname.split('.')[-1]:
        continue
    assert 'chr22' in fname
    fname = f'tests/data/tracks_chr22/{fname}'
    assert is_nonempty(fname, wig_min_bytes)
    with open(fname, 'r') as f_:
        line_ct = len(f_.readlines())
    assert line_ct >= min_wig_line_ct
    assert int(line_k(fname, 3).split('\t')[0]) - int(line_k(fname, 2).split('\t')[0]) == 50
log('all cases passed', sys.argv[0])