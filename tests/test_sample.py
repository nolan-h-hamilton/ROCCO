import os
import pytest
from rocco import Sample
import numpy as np

TEST_BAMFILE_chr21 = 'data/ENCFF231YYD_chr21.bam'
TEST_GENOME_FILE_chr21 = 'test_chr21.sizes'
TEST_STEP = 100

@pytest.mark.regular
def test_sample_init():
    global TEST_BAMFILE_chr21
    global TEST_GENOME_FILE_chr21
    global TEST_STEP
    test_sample = Sample(TEST_BAMFILE_chr21, TEST_GENOME_FILE_chr21, step=TEST_STEP)
    
    assert len(np.unique(np.diff(test_sample.get_chrom_data('chr21')[0]))) == 1
    assert TEST_STEP in list(np.unique(np.diff(test_sample.get_chrom_data('chr21')[0])))
    
test_sample_init()