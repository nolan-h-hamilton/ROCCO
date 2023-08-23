# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
1. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
1. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample count;
1. **Does not require training data** or initial candidate peak regions which are hard to define given the lack of a priori sets of open chromatin regions;
1. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.

## Getting Started

**Install Dependencies with Conda**

A ROCCO-specific conda environment with all dependencies installed can be created using
[rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/CONDA/rocco_conda.yml):

```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
load via: `conda activate rocco`.

Alternatively, dependencies (standard bioinformatics tools listed in `docs`) can be installed manually.

**Install ROCCO with pip** ([PyPi page](https://pypi.org/project/rocco/))

  ```
  pip install rocco
  ```


#### Quick Start Demo
To see ROCCO in action, refer to the Jupyter notebook: [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo.ipynb).

This demonstration offers an *interactive* overview of the ROCCO pipeline that can be executed by running the commands in each cell. Output from a previous session is included if you do not wish to run the pipeline yourself.

## Citation
```
ROCCO: A Robust Method for Detection of Open Chromatin via Convex Optimization
Nolan H. Hamilton, Terrence S. Furey
bioRxiv 2023.05.24.542132; doi: https://doi.org/10.1101/2023.05.24.542132
```
