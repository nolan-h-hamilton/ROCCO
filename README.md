# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/rocco?label=PyPI%20package&color=blue)](https://pypi.org/project/rocco/)

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
1. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
1. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample count;
1. **Does not require training data** or initial candidate peak regions which are hard to define given the lack of a priori sets of open chromatin regions;
1. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.

## Getting Started

ROCCO is written for compatibility with any unix variant: MacOS, linux, etc.

**Conda**

If using [conda](https://anaconda.org), ROCCO and all dependencies can be installed in a virtual environment specified in [rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/rocco_conda.yml):

```
conda env create -n rocco --file rocco_conda.yml
```

load via: `conda activate rocco`.

**Install ROCCO with pip** ([PyPI](https://pypi.org/project/rocco/))

Alternatively, ROCCO can be obtained via pip/PyPi:
  ```
  pip install rocco --upgrade
  ```

System Dependencies: *samtools, bedtools, [UCSC KentUtils](http://hgdownload.soe.ucsc.edu/admin/exe/): bigWigToWig, wigToBigWig, bigWigCat*

#### Quick Start Demo
For a quick intro, refer to the Jupyter notebook: [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb).

This demonstration offers an interactive overview of the ROCCO pipeline that can be executed by running the commands in each cell. Output from a previous session is included if you do not wish to run the steps yourself.
## Documentation
API Reference: https://nolan-h-hamilton.github.io/ROCCO/index.html
## Citation
```
ROCCO: A Robust Method for Detection of Open Chromatin via Convex Optimization
Nolan H. Hamilton, Terrence S. Furey
bioRxiv 2023.05.24.542132; doi: https://doi.org/10.1101/2023.05.24.542132
```
