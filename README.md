# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

<p align="center">
<img width="400" alt="logo" src="docs/logo.png">

## What

ROCCO is an algorithm for efficient identification of "consensus peaks" in multiple HTS data samples (namely, ATAC-seq), where read densities are consistently enriched.

## How

ROCCO models consensus peak calling as a constrained optimization problem with an upper-bound on the total proportion of the genome selected as open/accessible and a fragmentation penalty to promote spatially consistency in active regions and sparsity elsewhere.

## Why

ROCCO offers several attractive features:

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals
2. **Scaling to large sample sizes** with an asymptotic time complexity independent of sample size
3. **No required training data** or a heuristically determined set of initial candidate peak regions
4. **No rigid thresholds** on the minimum number/width of supporting samples/replicates
5. **Mathematically tractable model** with worst-case bounds on runtime and performance

### Demo

A brief walkthrough with visualized peak results using publicly available ATAC-seq data:

<https://github.com/nolan-h-hamilton/ROCCO/tree/main/docs/demo/demo.ipynb>

### Paper/Citation

If using ROCCO in your research, please cite the [original paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics*

   ```plaintext
    Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
    Bioinformatics, Volume 39, Issue 12, December 2023
   ```

### Documentation

ROCCO's documentation is available at <https://nolan-h-hamilton.github.io/ROCCO/>

Note that using the module-level functions directly may allow for greater flexibility in applications than using the command-line interface, which is limited in scope.

### Installation

#### PyPI (`pip`)

   ```shell
   pip install rocco
   ```

#### Build from Source

You can also build from source directly if preferred.

* Clone or download this repository

  ```shell
  git clone https://github.com/nolan-h-hamilton/ROCCO.git
  cd ROCCO
  # if you haven't already before, run `pip install setuptools wheel`
  python setup.py sdist bdist_wheel
  pip install -e .
  ```

ROCCO utilizes the popular bioinformatics software [Samtools](http://www.htslib.org) and [bedtools](https://bedtools.readthedocs.io/en/latest/). If not available already, these system dependencies can be installed with standard MacOS or Linux/Unix package managers, e.g., `brew install samtools` (Homebrew), `sudo apt-get install samtools` (APT).

### Input/Output

* *Input*: Samples' BAM alignments or BigWig tracks
* *Output*: BED file of consensus peak regions
