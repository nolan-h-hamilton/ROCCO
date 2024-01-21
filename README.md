# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

<p align="center">
<img width="400" alt="logo" src="docs/logo.png">

ROCCO is a multisample consensus peak caller for open chromatin data, e.g., ATAC-seq.

**Features**

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
2. **Mathematically tractable model** that permits performance and efficiency guarantees.
3. **Efficient for large numbers of samples** with an asymptotic time complexity independent of sample size;
4. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
5. **No required training data** or a heuristically determined set of initial candidate peak regions;


# Paper

If using ROCCO in your research, please cite the [original paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics*


   ```
    Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
    Bioinformatics, Volume 39, Issue 12, December 2023
   ```

**DOI**: ``10.1093/bioinformatics/btad725``

# Documentation

Documentation and example usage are available at https://nolan-h-hamilton.github.io/ROCCO/

# Installation

   ```
   pip install rocco
   ```

# Input
ROCCO accepts samples' **BAM** alignments or **BigWig** coverage tracks as input.

# Output

A **BED** file containing peak regions and scores.

# Minimal Example

Both an API and command-line interface are available to run ROCCO.

## Command-line interface

Run ROCCO on the test data included with this repository (BigWig files). Output will be stored in a BED6 file.

   ```
   rocco -i tests/data/*.bw --genome_file tests/test_hg38.sizes --chrom_param_file tests/test_hg38_param_file.csv
   ```

Default output BED file named as `rocco_peaks_[timestamp].bed`.

## API

Same as the above example, but using the API:

   ```
    >>> import rocco
    >>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw']
    >>> # see Rocco.HG38_PARAMS
    >>> rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='tests/test_hg38_param_file.csv')
    >>> rocco_obj.run() # Default output BED file named as `rocco_peaks_[timestamp].bed`.
   ```

# Testing ROCCO

  ```
  cd tests
  pytest -v -rPA -l -k "regular" test_rocco.py
  ```

# Version History

Previous releases can be found at https://github.com/nolan-h-hamilton/ROCCO/tags
