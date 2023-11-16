"""
Test 'rocco chrom' utility with pytest
"""

import os
import subprocess
import shutil

solver = 'CLARABEL'
peakfile_minsize = 100
jaccard_min = .90
rr_iter = '100'
chr19_start = '20000000'
chr19_end = '50000000'
chr19_budget= '0.035'
chr20_start = '20000000'
chr20_end = '50000000'
chr20_budget = '0.035'

tracks_chrom19_wig = os.path.relpath('data/tracks_chr19/chr19_test.bam.bw.wig.orig')
tracks_chrom20_wig = os.path.relpath('data/tracks_chr20/chr20_test.bam.bw.wig.orig')

chr19_outfile = f'ROCCO_out_chr19_{chr19_budget}_1.0_0.0_1.0_1.0_1.0.bed'
chr19_ref_file = os.path.relpath('data/ref/ref_chr19.bed')

chr20_outfile = f'ROCCO_out_chr20_{chr20_budget}_1.0_0.0_1.0_1.0_1.0.bed'
chr20_ref_file = os.path.relpath('data/ref/ref_chr20.bed')

gwide_outfile = 'gwide_out.bed'
gwide_outdir = os.path.relpath('output/combined')
gwide_ref_file = os.path.relpath('data/ref/ref_gwide.bed')


def test_sim(var=1, samples=10):
    """
    Add noise to existing wiggle file to generate varying tracks for testing purposes.
    `samples` signal tracks are generated.
    """
    proc1 = subprocess.run(['python', 'sim.py', tracks_chrom19_wig, str(var), str(samples)],
                           capture_output=True, encoding='utf-8')
    proc2 = subprocess.run(['python', 'sim.py', tracks_chrom20_wig, str(var), str(samples)],
                           capture_output=True, encoding='utf-8')
    assert proc1.returncode == 0, f'sim.py failed for {tracks_chrom19_wig}'
    assert proc2.returncode == 0, f'sim.py failed for {tracks_chrom20_wig}'


def test_rocco_exists():
    """
    Ensure command line utility 'rocco' is available
    """
    assert shutil.which('rocco') is not None

def test_chrom_help():
    """
    Ensure 'rocco chrom -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'chrom', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco chrom' in proc.stdout

def test_gwide_help():
    """
    Ensure 'rocco gwide -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'gwide', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco gwide' in proc.stdout

def test_budgets_help():
    """
    Ensure 'rocco budgets -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'budgets', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco budgets' in proc.stdout

def test_prep_help():
    """
    Ensure 'rocco prep -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'prep', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco prep' in proc.stdout

def test_run_chrom_chr19():
    """
    Ensure 'rocco chrom' runs successfully for chr19
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig),
                           '--budget', chr19_budget, '-N', rr_iter, '--start', chr19_start, '--end', chr19_end,
                           '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco chrom returncode != 0'
    # check if output file was created properly
    assert chr19_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(chr19_outfile).st_size > peakfile_minsize


def test_run_chrom_chr20():
    """
    Ensure 'rocco chrom' runs successfully for chr20
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom', '--chrom', 'chr20', '--wig_path', os.path.dirname(tracks_chrom20_wig),
                            '--budget', chr20_budget, '-N', rr_iter, '--start', chr20_start, '--end', chr20_end,  '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0,'rocco chrom returncode != 0'
    # check if output file was created properly
    assert chr20_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(chr20_outfile).st_size > peakfile_minsize

def test_run_gwide():
    # check if command returns exit code '0'
    assert subprocess.run(["rocco", "gwide", "-p", "test_params.csv","--outdir", gwide_outdir, "--combine", gwide_outfile, "-N", rr_iter,  '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco gwide returncode != 0'
    # check if output file was created properly
    assert gwide_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(gwide_outfile).st_size > peakfile_minsize
    chrom_files = [x for x in os.listdir(gwide_outdir) if x.split('.')[-1] == 'bed']
    # ensure outdir contains chromosome-specific BED files
    assert len(chrom_files) == 2

def test_comp_ref_chr19():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", chr19_outfile , "-b", chr19_ref_file ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_chr19 = float(data.split("\t")[2])
    assert jaccard_chr19 >= jaccard_min

def test_comp_ref_chr20():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", chr20_outfile , "-b", chr20_ref_file ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_chr20 = float(data.split("\t")[2])
    assert jaccard_chr20 >= jaccard_min

def test_comp_ref_gwide():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", gwide_outfile, "-b", gwide_ref_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_gwide = float(data.split("\t")[2])
    assert jaccard_gwide >= jaccard_min
