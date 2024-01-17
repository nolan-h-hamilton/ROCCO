# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

<p align="center">
<img width="400" alt="logo" src="docs/logo.png">

Underlying ROCCO is a constrained optimization problem that can be
solved efficiently to predict consensus regions of open chromatin across
multiple samples.

**Features**

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
2. **Mathematically tractable model** that permits performance and efficiency guarantees.
3. **Efficient for large numbers of samples** with an asymptotic time complexity independent of sample size;
4. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
5. **No required training data** or a heuristically determined set of initial candidate peak regions;


# Paper

If using ROCCO in your research, please cite the [original paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics*.

# Installation

   ```
   pip install rocco
   ```

# GitHub (Homepage)

https://github.com/nolan-h-hamilton/ROCCO/

# Documentation

Documentation and example use demos are available at https://nolan-h-hamilton.github.io/ROCCO/

# Input
ROCCO accepts samples' **BAM** alignments or **BigWig** coverage tracks as input.

# Output

A **BED** file containing peak regions and scores.

# Minimal Example

Both an API and command-line interface are available to run ROCCO.

## Command-line interface

Run ROCCO on the test data included with this repository. Output will be stored in a BED6 file.

   ```
   rocco -i tests/data/*.bw --genome_file tests/test_hg38.sizes --chrom_param_file tests/test_hg38_param_file.csv
   ```

## API

Same as the above example, but using the API:

   ```
    >>> import rocco
    >>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw', 'tests/data/sample4.bw', 'tests/data/sample5.bw']
    >>> # see Rocco.HG38_PARAMS
    >>> rocco_obj = Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='tests/test_hg38_param_file.csv')
    >>> rocco_obj.run() # genome-wide output stored in BED6 file
   ```

# Testing ROCCO

  ```
  cd tests
  pytest -v -rPA -l -k "regular" test_rocco.py
  ```

# Notes/Miscellaneous

* If using BedGraph or BigWig input, ensure contiguous intervals within each chromosome (no gaps)

* Users may consider tweaking the default chromosome-specific $b,\gamma,\tau$ parameters or filtering peaks by score with the `--peak_score_filter` argument.

* Peak scores are computed as the average number of reads over the given peak region (w.r.t samples), divided by the length of the region, and then scaled to units of kilobases. If you wish to incorporate the peak score filter, a suitable cutoff can be evaluated by viewing the output histogram of peak scores.


# Version History

Previous releases can be found at https://github.com/nolan-h-hamilton/ROCCO/tags