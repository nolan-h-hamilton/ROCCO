import pybedtools
import pytest
import rocco
import subprocess
import os

MIN_JACCARD = .95
PEAK_SCORE_MIN = 500
BW_FILES = ['data/sample1.bw', 'data/sample2.bw', 'data/sample3.bw']
PR_BED = 'test_pr.bed'

@pytest.mark.regular
def test_rocco_base():
    global BW_FILES
    global MIN_JACCARD
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv', step=100)
    rocco_obj.run() # genome-wide output stored in BED6 file
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD
    os.remove(rocco_obj.outfile)

@pytest.mark.extra
def test_peakscore_filter():
    global BW_FILES
    global PEAK_SCORE_MIN
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv', peak_score_filter = PEAK_SCORE_MIN)
    rocco_obj.run() # genome-wide output stored in BED6 file
    for feature in pybedtools.BedTool(rocco_obj.outfile):
        assert float(feature[4]) >= PEAK_SCORE_MIN, f"{feature}, {PEAK_SCORE_MIN}"
    os.remove(rocco_obj.outfile)

@pytest.mark.extra
def test_constant_params():
    global BW_FILES
    global MIN_JACCARD
    filler_params = {'budget':0.045,
                     'gamma':1.0,
                     'tau':0.0,
                     'c_1':1.0,
                     'c_2':1.0,
                     'c_3':1.0}
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', filler_params=filler_params,
                            outfile='constant_param_test.bed')
    rocco_obj.run()
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD
    os.remove(rocco_obj.outfile)

@pytest.mark.extra
def test_pr_exclusion():
    global PR_BED
    global BW_FILES
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes',
                            outfile='test_pr_exclusion.bed', pr_bed=PR_BED)
    rocco_obj.run()
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b=PR_BED)['jaccard']) == 0.0
    os.remove(rocco_obj.outfile)
