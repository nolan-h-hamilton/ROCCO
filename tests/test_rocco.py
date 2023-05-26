"""test ROCCO for some basic use-cases"""
import os
import subprocess
import pysam
import pandas as pd



# move to ROCCO main directory
if os.getcwd().split('/')[-1] != 'ROCCO':
    os.chdir('..')
assert os.getcwd().split('/')[-1] == 'ROCCO','script needs to be run from ROCCO/tests'
init_files = os.listdir('tests')
# create index files for each BAM file
for file_ in os.listdir('tests'):
    if file_.split('.')[-1] == 'bam' and not os.path.exists('tests/' + file_ + '.bai'):
        try:
            pysam.index('tests/'+file_, 'tests/'+file_ + '.bai')
            print('indexing {}'.format('tests/'+file_))
        except pysam.SamtoolsError as ex:
            print('Could not create index file for BAM via pysam.\
                  This will cause errors for the count_matrix.py test')
            print(ex)

print('\n\ncase 1: run ROCCO using all samples - no parallel jobs')
cmd = ['python3', 'ROCCO.py',
        '--param_file', 'tests/test_params.csv',
          '--outdir', 'tests',
            '--combine', 'tests/test_out1.bed']
with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
    for line in p.stdout:
        print(line, end='')
    p.wait()
assert p.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_out1.bed'), 'bed file not created'
with open('tests/test_out1.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, 'could not find chr21 results in test bed'
    assert 'chr22' in bed_content, 'could not find chr22 results in test bed'


print('\n\ncase II: run ROCCO on `group_a` samples - no parallel jobs')
cmd = ['python3', 'ROCCO.py',
        '--param_file', 'tests/test_params.csv',
        '--outdir', 'tests',
        '--combine', 'tests/test_group_a_out.bed',
        '--identifiers', 'tests/test_group_a.txt']
with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
    for line in p.stdout:
        print(line, end='')
    p.wait()
assert p.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_group_a_out.bed'), 'bed file not created'
with open('tests/test_group_a_out.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, 'could not find chr21 results in test bed'
    assert 'chr22' in bed_content, 'could not find chr22 results in test bed'


print('\n\ncase III: run ROCCO on all samples with parallel jobs for each chrom. create count matrix')
cmd = ['python3', 'ROCCO.py',
        '--param_file', 'tests/test_params.csv',
        '--outdir', 'tests',
        '--combine', 'tests/test_par_out.bed',
        '--jobs', '2']
with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
    for line in p.stdout:
        print(line, end='')
    p.wait()
assert p.returncode == 0, 'job failed: {}'.format(cmd)
assert os.path.exists('tests/test_par_out.bed'), 'bed file not created'
with open('tests/test_par_out.bed', mode='r', encoding="utf-8") as f:
    bed_content = f.read()
    assert 'chr21' in bed_content, 'could not find chr21 results in test bed'
    assert 'chr22' in bed_content, 'could not find chr22 results in test bed'
cmd = ['python3', 'count_matrix.py',
        '-i', 'tests/test_par_out.bed',
        '-d', 'tests/test_metadata.txt',
        '--bamdir', 'tests',
        '-o', 'tests/test_countmat.tsv']
with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
    for line in p.stdout:
        print(line, end='')
    p.wait()
cmatrix = pd.read_csv('tests/test_countmat.tsv', sep='\t')
samps = pd.read_csv('tests/test_metadata.txt', sep='\t')
cols1 = list(cmatrix.columns)
cols1.remove('name')
IDs = list(samps['name'])
assert os.path.exists('tests/test_countmat.tsv')
assert cols1 == IDs

# delete test files
testfiles = [x for x in os.listdir('tests') if x not in init_files]
for test_file in testfiles:
    os.remove('tests/' + test_file)
print('test_rocco.py: all cases passed.')
