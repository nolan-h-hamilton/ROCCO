# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)
![PyPI - Version](https://img.shields.io/pypi/v/rocco?logo=Python&logoColor=%23FFFFFF&color=%233776AB&link=https%3A%2F%2Fpypi.org%2Fproject%2Frocco%2F)

## What

ROCCO is an efficient algorithm for detection of "consensus peaks" in large datasets with multiple HTS data samples, where an enrichment in read counts/densities is observed in a nontrivial subset of samples.

### Input/Output

* *Input*: Samples' BAM alignments
* *Output*: BED file of consensus peak regions (Default format is BED3: `chrom,start,end`)

## How

ROCCO models consensus peak calling as a constrained optimization problem with an upper-bound on the total proportion of the genome selected as open/accessible and a fragmentation penalty to promote spatial consistency in active regions and sparsity elsewhere.

## Why

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals
2. **Scaling to large sample sizes (100+)**
3. **Unsupervised** Does not require training data or a heuristically determined set of initial candidate peak regions
4. **No rigid thresholds + less manual tuning** with respect to the minimum number/width of supporting samples/replicates. (A polytime 2-state DP now solves the constrained optimization problem, and budget $b$, fragmentation/TV penalty $\gamma$ are data-driven by default.)
5. **Mathematically tractable model** permitting worst-case analysis of runtime and performance

## Usage

  ```shell
  rocco -i <bam_files> -g <hg38, hg19, mm10, mm39, dm6, ...> -o <output_file.bed> --narrowPeak
  ```

for example:

  ```shell
  rocco -i sample1.bam sample2.bam sample3.bam -g hg38 -o consensus_peaks.bed --narrowPeak
  ```

See `rocco --help` for more options and details.


## Paper/Citation

If using ROCCO in your research, please cite the [original paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics* (DOI: `btad725`)

   ```plaintext
    Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
    Bioinformatics, Volume 39, Issue 12, December 2023
   ```

## Installation

### PyPI (`pip`)

   ```shell
   python -m pip install rocco --upgrade
   ```

If lacking administrative control, you may need to append `--user` to the above.


### Build from Source

If preferred, ROCCO can easily be built from source:

* Clone or download this repository

  ```shell
  git clone https://github.com/nolan-h-hamilton/ROCCO.git
  cd ROCCO
  python setup.py sdist bdist_wheel
  python -m pip install -e .
  ```
