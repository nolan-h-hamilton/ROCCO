import os
import pytest
from rocco import Sample
import numpy as np

TEST_BAMFILE_chr21 = 'data/ENCFF231YYD_chr21.bam'
TEST_GENOME_FILE_chr21 = 'test_chr21.sizes'
TEST_STEP = 100
TEST_CHR21_OUTFILE = 'test_out_ENCFF231YYD_chr21.bg'

@pytest.mark.regular
def test_sample_init():
    global TEST_BAMFILE_chr21
    global TEST_GENOME_FILE_chr21
    global TEST_STEP
    global TEST_CHR21_OUTFILE
    test_sample = Sample(TEST_BAMFILE_chr21, TEST_GENOME_FILE_chr21, step=TEST_STEP, output_file=TEST_CHR21_OUTFILE)
    test_sample.write_track()
    
    assert len(np.unique(np.diff(test_sample.coverage_dict['chr21'][0]))) == 1
    assert test_sample.coverage_dict['chr21'][0][1] - test_sample.coverage_dict['chr21'][0][0] == TEST_STEP
    assert test_sample.output_format in ['bg','bedgraph']
    os.remove(TEST_CHR21_OUTFILE)