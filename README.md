# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals, the latter of which is an informative but often ignored aspect of ATAC-seq data;
1. Leverages data from multiple samples of varying quality **without imposing arbitrary thresholds** on a minimum number of samples declaring peaks;
1. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample count;
1. **Does not require training data or a heuristically determined set of initial candidate regions**, which are hard to define given the lack of a priori sets of open chromatin regions;
1. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.

## Environment Details
ROCCO has been developed and tested to run on ATAC-seq alignments in a standard unix bioinformatics environment with Python3.7+

A ROCCO-specific conda environment with all dependencies installed can be created using
[rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/CONDA/rocco_conda.yml):
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
The created environment can then be loaded via: `conda activate rocco`.

## Input
- A BAM file for each sample

The signal matrix, $\mathbf{S}_{chr} \in \mathbb{R}^{K\times n}$, as defined in the paper is created from the BAM files by running:
```
python3 prep_bams.py -i <bam_dir, def=.> --cores <num_cores, def=1> -L <step_size, def=50>
```
A [flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png) provides a visual demonstration of this preprocessing step. Run `python3 prep_bams.py --help` for a full list of parameters.

## Getting Started
Clone/download this repository to use ROCCO
  ```
  git clone https://github.com/nolan-h-hamilton/ROCCO.git
  ```
### Jupyter Notebook Demos
  Some lightweight demonstrations using publicly available data. 
  1. **Quick Start Demo**. [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo.ipynb). Includes visualized results in UCSC genoome browser.
  1. **Differential Accessibility**. [heart_da_demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo_files/heart_da_demo.ipynb). A sample ROCCO-->[DESeq2](https://github.com/mikelove/DESeq2) pipeline for differential analysis.

### Auxiliary Scripts
Used for pre/post-processing. Use cases included in the demo notebooks and/or the scripts' module docstrings.
1. `est_budgets.py`: used to determine chromosome-specific budgets as a function of read density
1. `count_matrix.py`: used to create a count matrix for ROCCO peak results
   
### ROCCO API Reference
Source code documentation is available [here](https://nolan-h-hamilton.github.io/ROCCO/).
