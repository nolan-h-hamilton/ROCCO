# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

<p align="center">
<img width="400" alt="logo" src="docs/logo.png">

ROCCO is a scalable consensus peak calling algorithm for open chromatin data on multiple samples, e.g., ATAC-seq.

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
ROCCO accepts samples' **BAM** alignments or **BigWig** coverage tracks and a genome sizes file as input.

# Output

A **BED** file containing peak regions and scores.

# Minimal Example

Run ROCCO on the test data included with this repository (BigWig files).

   ```
   rocco -i tests/data/*.bw --genome_file tests/test_hg38.sizes
   ```

Default output BED file named as `rocco_peaks_[timestamp].bed`. Because no `--chrom_param_file` is supplied,
the default genome-wide budget (`0.035`), gamma (`1.0`), etc. are used for all chromosomes. See the [Documentation](https://nolan-h-hamilton.github.io/ROCCO/) for additional examples and details.

# Testing ROCCO

  ```
  cd tests
  pytest -v -rPA -l -k "regular" test_rocco.py
  ```

# Version History

Previous releases can be found at https://github.com/nolan-h-hamilton/ROCCO/tags
