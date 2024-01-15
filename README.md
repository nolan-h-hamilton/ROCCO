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

See `rocco --help` for a full list of argument descriptions.

### Example One

Run ROCCO with BAM input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for `hg38` assembly
in `Rocco.HG38_PARAMS`

   ```
   rocco --input_files tests/data/*.bam --genome_file tests/test_hg38.sizes --chrom_param_file hg38
   ```

If BAM input is supplied, [`pysam`](https://pysam.readthedocs.io/en/stable/) is used to compute the coverage tracks marking
the number of reads aligning to each bin/locus at interals of size `--step`.


### Example Two

Run ROCCO with bigwig input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for `hg38` assembly
in `Rocco.HG38_PARAMS`

For instance, if you wish to generate the coverage tracks with
[deepTools
bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
or another utility that produces bigwig signal files with additional
features for normalization, smoothing, read extension, etc. You can
supply the resulting bigwig files as input to ROCCO.


```
rocco --input_files tests/data/*.bw --genome_file tests/test_hg38.sizes --chrom_param_file hg38
```

### Example Three

Run ROCCO with bedgraph input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for the `hg38`
assembly in `Rocco.HG38_PARAMS`

```
rocco --input_files tests/data/*.bg --genome_file tests/test_hg38.sizes --chrom_param_file hg38
```

### Example Four

Scale coverage value of samples before calling peaks

```
rocco --input_files tests/data/*.bw --genome_file tests/test_hg38.sizes --chrom_param_file hg38 --sample_weights 0.50 1.0 1.0
```

### Example Five

Use a custom chromosome parameter file

```
rocco --input_files tests/data/*.bw --genome_file tests/test_hg38.sizes --chrom_param_file tests/test_hg38_param_file.csv
```

## API (Application Programmer Interface)

### Example One

Run ROCCO with BAM input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for `hg38` assembly
in `Rocco.HG38_PARAMS`

```
>>> import rocco
>>> bamfiles = ['tests/data/sample1.bam', 'tests/data/sample2.bam', 'tests/data/sample3.bam']
>>> rocco_obj = rocco.Rocco(input_files=bamfiles, genome_file='tests/test_hg38.sizes', chrom_param_file='hg38') # see Rocco.HG38_PARAMS
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```

### Example Two

Run ROCCO with bigwig input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for the `hg38`
assembly in `Rocco.HG38_PARAMS`

```
>>> import rocco
>>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw']
>>> rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='hg38') # see Rocco.HG38_PARAMS
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```

### Example Three

Run ROCCO with bedgraph input files for each sample using default
chromosome-specific budget, gamma, etc. parameters for the `hg38`
assembly in `Rocco.HG38_PARAMS`

```
>>> import rocco
>>> bedgraph_files = ['tests/data/sample1.bg', 'tests/data/sample2.bg', 'tests/data/sample3.bg']
>>> rocco_obj = rocco.Rocco(input_files=bedgraph_files, genome_file='tests/test_hg38.sizes', chrom_param_file='hg38') # see ROCCO.HG38_PARAMS
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```

### Example Four

Scale coverage value of sample1 before calling peaks

```
>>> import rocco
>>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw']
>>> rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='hg38', sample_weights=[0.50,1.0,1.0])
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```

### Example Five

Use a custom chromosome parameter file

```
>>> import rocco
>>> bw_files = ['tests/data/sample1.bw', 'tests/data/sample2.bw', 'tests/data/sample3.bw']
>>> rocco_obj = rocco.Rocco(input_files=bw_files, genome_file='tests/test_hg38.sizes', chrom_param_file='tests/test_hg38_param_file.csv')
>>> rocco_obj.run() # genome-wide output stored in BED6 file
```
