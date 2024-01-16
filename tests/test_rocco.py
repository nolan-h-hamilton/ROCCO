import pybedtools
import pytest
import rocco

MIN_JACCARD = 0.95
PEAK_SCORE_MIN = 200

@pytest.mark.regular
def test_rocco_base():
    bw_files = ['data/sample1.bw', 'data/sample2.bw', 'data/sample3.bw']
    rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv')
    rocco_obj.run() # genome-wide output stored in BED6 file
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD

@pytest.mark.regular
def test_peakscore_filter():
    bw_files = ['data/sample1.bw', 'data/sample2.bw', 'data/sample3.bw']
    rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv', peak_score_filter = PEAK_SCORE_MIN)
    rocco_obj.run() # genome-wide output stored in BED6 file
    for feature in pybedtools.BedTool(rocco_obj.outfile):
        assert float(feature[4]) >= PEAK_SCORE_MIN, f"{feature}, {PEAK_SCORE_MIN}"
