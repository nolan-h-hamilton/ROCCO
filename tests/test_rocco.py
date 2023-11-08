"""
Test 'rocco chrom' utility with pytest
"""

import os
import subprocess
import shutil

peakfile_minsize = 100
tracks_chrom19_wig = os.path.relpath('data/tracks_chr19/chr19_test.bam.bw.wig.orig')
tracks_chrom20_wig = os.path.relpath('data/tracks_chr20/chr20_test.bam.bw.wig.orig')

def test_sim():
    """
    Simulate 5 signal tracks per chromosome for testing
    """
    proc1 = subprocess.run(['python', 'sim.py', tracks_chrom19_wig, '10', '5'],
                           capture_output=True, encoding='utf-8')
    proc2 = subprocess.run(['python', 'sim.py', tracks_chrom20_wig, '10', '5'],
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

def test_chrom_chr19():
    """
    Ensure 'rocco chrom' runs successfully for chr19
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig), '--budget', '0.05',
                           '--start', '20000000', '--end', '30000000'],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco chrom returncode != 0'
    # check if output file was created properly
    assert 'ROCCO_out_chr19_0.05_1.0_0.0_1.0_1.0_1.0.bed' in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat('ROCCO_out_chr19_0.05_1.0_0.0_1.0_1.0_1.0.bed').st_size > peakfile_minsize

def test_chrom_chr20():
    """
    Ensure 'rocco chrom' runs successfully for chr20
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom', '--chrom', 'chr20', '--wig_path', os.path.dirname(tracks_chrom20_wig), '--budget', '0.04',
                           '--start', '20000000', '--end', '30000000'],
                          capture_output=True, encoding='utf-8').returncode == 0,'rocco chrom returncode != 0'
    # check if output file was created properly
    assert 'ROCCO_out_chr20_0.04_1.0_0.0_1.0_1.0_1.0.bed' in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat('ROCCO_out_chr20_0.04_1.0_0.0_1.0_1.0_1.0.bed').st_size > peakfile_minsize
