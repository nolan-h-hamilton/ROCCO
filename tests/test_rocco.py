import pybedtools
import pytest
import rocco

MIN_JACCARD = 0.95

@pytest.mark.regular
def test_rocco_base():
    bw_files = ['data/sample1.bw', 'data/sample2.bw', 'data/sample3.bw']
    rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='test_hg38.sizes', chrom_param_file='test_hg38_param_file.csv')
    rocco_obj.run() # genome-wide output stored in BED6 file
    assert float(pybedtools.BedTool(rocco_obj.outfile).jaccard(b='data/ref/test_ref.bed')['jaccard']) > MIN_JACCARD
