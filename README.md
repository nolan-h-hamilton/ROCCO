# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)
![PyPI - Version](https://img.shields.io/pypi/v/rocco?logo=Python&logoColor=%23FFFFFF&color=%233776AB&link=https%3A%2F%2Fpypi.org%2Fproject%2Frocco%2F)

## What

ROCCO is an efficient algorithm for detection of "consensus peaks" in large datasets with multiple HTS data samples, where an enrichment in read counts/densities is observed in a nontrivial subset of samples.

### Input/Output

* *Input*: Samples' `.bam` alignments (`-i`) and a reference genome assembly (`-g`) for chromosome sizes
* *Output*: BED file of consensus peak regions (Default format is BED3: `chrom,start,end`). narrowPeak and gappedPeak outputs are supported for BAM inputs. Use `--peak_mode both`.

## How

ROCCO models consensus peak calling as an optimization problem: select the most enriched genomic regions as peaks while controlling total genomic coverage and discouraging fragmented calls.

* enrichment is scored: increasing with read density and decreasing with dispersion
* budget proportions calibrate selection penalties that control total genomic coverage
* (excessive) fragmentation is controlled via a total variation penalty, multiplied by a penalty parameter ($\gamma \sum_i |x_{i+1} - x_i|$).


## Why

1. **Consideration of enrichment magnitude and spatial characteristics** of open chromatin signals
2. **Scaling to large sample sizes (100+)**
3. **Unsupervised** Does not require training data or a heuristically determined set of initial candidate peak regions
4. **Less rigid thresholds** with respect to the minimum number/width of supporting samples/replicates.
5. **Mathematically tractable model** permitting worst-case analysis of runtime and performance

## Basic Usage

Run `rocco --help` for the complete command-line reference.

ROCCO accepts either one or more BAM files, separated by spaces.

```shell
rocco -i <BAM> [<BAM> ...] \
  -g <GENOME> \
  -o <OUTPUT.bed>
```
Genome assemblies `hg38`, `hg19`, `mm10`, `mm39`, and `dm6` can be used with `-g`, to supply a built-in chromosome sizes file (canonical chromosomes only) and effective genome size.

If using ROCCO with another assembly, just provide the chromosome sizes file and effective genome size as described explicitly:

* `-s, --chrom_sizes_file <FILE>`: chromosome names and sizes.
* `--effective_genome_size <BASES>`: [effective genome size](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html), only required if using RPGC normalization.

For example:

```shell
rocco \
  -i sample1.bam sample2.bam sample3.bam \
  -s assemblyName.chrom.sizes \
  --effective_genome_size <BASES> \
  -o consensus_peaks.bed
```

When troubleshooting, consider processing a subset of
the genome for faster evaluation:

```shell
rocco \
  -i sample1.bam sample2.bam sample3.bam \
  -g hg38 \
  -o consensus_peaks.bed \
  --chroms chr11 chr20 chr21 chr22 # at least four chroms
```

### Example

Call consensus peaks from three BAM files and generate both narrowPeak and gappedPeak outputs:

```shell
rocco \
  -i sample1.bam sample2.bam sample3.bam \
  -g hg38 \
  -o consensus_peaks.bed \
  --peak_mode both
```

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
