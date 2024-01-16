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

# Example Usage

ROCCO offers a command-line interface for convenience and also an API for greater programmatic flexibility.

## CLI (Command-line Interface)

See `rocco --help` for a full list of argument descriptions. Wildcards and regular expressions can be used to specify subsets of input files, chromosomes to skip, etc.

### Example 1

* BAM input files
* Default chromosome-specific parameters for hg38 (See code `Rocco.HG38_PARAMS`)

```
rocco --input_files sample1.bam sample2.bam sample3.bam --genome_file genome.sizes --chrom_param_file hg38
```

If BAM input is supplied, [`pysam`](https://pysam.readthedocs.io/en/stable/) is used to compute the coverage tracks marking
the number of reads aligning to each bin/locus at interals of size `--step`. This can require substantial computing time. See
`--proc_num` and `--samtools_threads` for multiprocessing options.


### Example 2

* BigWig input files (Specified with a wildcard)
* Default chromosome-specific parameters for hg38 (See code `Rocco.HG38_PARAMS`)

```
rocco --input_files *.bw --genome_file genome.sizes --chrom_param_file hg38
```

This input format is useful if you have used, e.g., [deepTools
bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html),
for normalization, smoothing, read extension, etc. the samples' initial BAM alignments.

### Example 3

* BedGraph input files (Specified with a wildcard)
* Default chromosome-specific parameters for hg38

```
rocco --input_files *.bg --genome_file genome.sizes --chrom_param_file hg38
```

### Example 4

* BAM input files
* Default chromosome-specific parameters for hg38
* Scale coverage tracks for each sample individually with `--sample_weights`

```
rocco --input_files sample1.bam sample2.bam sample3.bam \
      --genome_file genome.sizes --chrom_param_file hg38 \
      --sample_weights 1.50 1.0 1.0
```

### Example 5

* Use a custom chromosome parameter file (`tests/test_hg38_param_file.csv`)

```
rocco --input_files tests/data/sample1.bw tests/data/sample2.bw \
      tests/data/sample3.bw --genome_file tests/test_hg38.sizes \
      --chrom_param_file tests/test_hg38_param_file.csv
```

## API (Application Programmer Interface)


### Example 6

```
>>> import rocco
>>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw']
>>> rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='tests/test_hg38_param_file.csv')
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```

# Documentation

ROCCO's complete documentation is available at https://nolan-h-hamilton.github.io/ROCCO/


# Testing ROCCO

Run unit tests

  ```
  cd tests
  pytest -v -rPA -l -k "regular" test_rocco.py
  ```

# Notes/Miscellaneous


* If using BedGraph or BigWig input, ensure contiguous intervals within each chromosome (no gaps)
    * Such gaps can be filled with zeros.

* Users may consider tweaking the default chromosome-specific $b,\gamma,\tau$ parameters or filtering peaks by score with
    the `--peak_score_filter` argument.

* Peak scores are computed as the average number of reads over the given peak region (w.r.t samples), divided by the length of the region, and then scaled to units of kilobases. A suitable peak score cutoff can be evaluated by viewing the output histogram of peak scores. In many cases a cutoff around 100.0 is a reasonable starting point.


# Version History

Previous releases can be found at https://github.com/nolan-h-hamilton/ROCCO/tags