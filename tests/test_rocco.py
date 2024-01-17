import pybedtools
import pytest
import rocco

MIN_JACCARD = 0.95
PEAK_SCORE_MIN = 2000
BW_FILES = ['data/sample1.bw', 'data/sample2.bw', 'data/sample3.bw', 'data/sample4.bw', 'data/sample5.bw']

@pytest.mark.regular
def test_rocco_base():
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv')
    rocco_obj.run() # genome-wide output stored in BED6 file
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD

@pytest.mark.regular
def test_peakscore_filter():
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv', peak_score_filter = PEAK_SCORE_MIN)
    rocco_obj.run() # genome-wide output stored in BED6 file
    for feature in pybedtools.BedTool(rocco_obj.outfile):
        assert float(feature[4]) >= PEAK_SCORE_MIN, f"{feature}, {PEAK_SCORE_MIN}"

@pytest.mark.regular
def test_constant_params():
    filler_params = {'budget':0.045,
                     'gamma':1.0,
                     'tau':0.0,
                     'c_1':1.0,
                     'c_2':1.0,
                     'c_3':1.0}
    rocco_obj = rocco.Rocco(input_files=BW_FILES, genome_file='test_hg38.sizes', filler_params=filler_params)
    rocco_obj.run()
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD