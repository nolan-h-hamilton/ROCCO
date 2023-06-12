"""
Test ROCCO for some basic use-cases using toy data.

Note: this script will delete any existing bai/bed/tsv files in `tests` dir.
      with the prefix 'test_' or 'ROCCO_'
"""
import os
import subprocess
import pysam
import pandas as pd

def index_bamfiles():
    for file_ in os.listdir('tests'):
        if (
                file_.split('.')[-1] == 'bam' and
                not os.path.exists('tests/' + file_ + '.bai')):
            try:
                pysam.index('tests/' + file_, 'tests/' + file_ + '.bai')
                print('indexing {}'.format('tests/' + file_))
            except pysam.SamtoolsError as ex:
                print('Could not create index file for BAM via pysam.\
                  This will cause errors for the count_matrix.py test')
                print(ex)

def clean():
    for fname in os.listdir('tests'):
        if (fname.split('.')[-1] in ['bed','tsv', '.bai']
            and fname.split('_')[0] in ['test','ROCCO']):
            print('removing tests/{}'.format(fname))
            os.remove('tests/' + fname)


# move to ROCCO main directory
if os.getcwd().split('/')[-1] != 'ROCCO':
    os.chdir('..')
assert os.getcwd().split('/')[-1] == 'ROCCO', 'run script from ROCCO/tests'

# remove ANY bed or tsv files with prefix 'test_'
clean()

# create index files for each BAM file
# these indexes are necessary for III(b)
index_bamfiles()

# CASE 1
print('\n\ncase I: run ROCCO using all samples - no parallel jobs')

cmd = ['python3', 'ROCCO.py',
       '--param_file', 'tests/test_params.csv',
       '--outdir', 'tests',
       '--combine', 'tests/test_out1.bed']

with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as proc_one:
    for line in proc_one.stdout:
        print(line, end='')
    proc_one.wait()

assert proc_one.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_out1.bed'), 'bed file not created'
with open('tests/test_out1.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, "couldn't find chr21 results in bed"
    assert 'chr22' in bed_content, "couldn't find chr22 results in bed"


# CASE 2
print('\n\ncase II: run ROCCO on `group_a` samples - no parallel jobs')

cmd = ['python3', 'ROCCO.py',
       '--param_file', 'tests/test_params.csv',
       '--outdir', 'tests',
       '--combine', 'tests/test_group_a_out.bed',
       '--identifiers', 'tests/test_group_a.txt']

with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as proc_two:
    for line in proc_two.stdout:
        print(line, end='')
    proc_two.wait()

assert proc_two.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_group_a_out.bed'), 'bed file not created'
with open('tests/test_group_a_out.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, "couldn't find chr21 results in bed"
    assert 'chr22' in bed_content, "couldn't find chr22 results in bed"

# CASE 3
print('\n\ncase III: all samples, parallel jobs for each chrom,\
create count matrix')

cmd = ['python3', 'ROCCO.py',
       '--param_file', 'tests/test_params.csv',
       '--outdir', 'tests',
       '--combine', 'tests/test_par_out.bed',
       '--jobs', '2']

with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as proc_three:
    for line in proc_three.stdout:
        print(line, end='')
    proc_three.wait()

assert proc_three.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_par_out.bed'), 'bed file not created'
with open('tests/test_par_out.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, "couldn't find chr21 results in bed"
    assert 'chr22' in bed_content, "couldn't find chr22 results in bed"

print('\ncase III(b): creating count_matrix')

cmd = ['python3', 'count_matrix.py',
       '-i', 'tests/test_par_out.bed',
       '-d', 'tests/test_metadata.txt',
       '--bamdir', 'tests',
       '-o', 'tests/test_countmat.tsv']

with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as proc_four:
    for line in proc_four.stdout:
        print(line, end='')
    proc_four.wait()

cmatrix = pd.read_csv('tests/test_countmat.tsv', sep='\t')
samps = pd.read_csv('tests/test_metadata.txt', sep='\t')
cmatrix_cols = list(cmatrix.columns)
cmatrix_cols.remove('name')
samp_names = list(samps['name'])
assert os.path.exists('tests/test_countmat.tsv')
assert cmatrix_cols == samp_names
#clean()

# if we made it here, all cases passed
print('test_rocco.py: all cases passed.')

