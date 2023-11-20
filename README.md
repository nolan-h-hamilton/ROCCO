# [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/4e662eda-4899-46fa-ac9f-998e28f592c1">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/rocco?label=PyPI%20package&color=blue)](https://pypi.org/project/rocco/)

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
2. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
3. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample size;
4. **Does not require training data** or initial candidate peak regions
5. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency

## Installation

ROCCO is written for compatibility with -nix variants: MacOS, linux, etc.


**Install ROCCO with pip** ([PyPI](https://pypi.org/project/rocco/))

ROCCO can be installed with pip:
  ```
  pip install rocco --upgrade
  ```

System Dependencies: *samtools, bedtools, [UCSC KentUtils](http://hgdownload.soe.ucsc.edu/admin/exe/): bigWigToWig, wigToBigWig, bigWigCat*. All of these are available on conda.  Alternatively, for details
on manual installation, see [environment.md](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/environment.md)

**Conda**

If using [conda](https://anaconda.org), ROCCO and all core/optional dependencies can be installed in a virtual environment specified in [rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/rocco_conda.yml):

```
conda env create -n rocco --file rocco_conda.yml
```

load via: `conda activate rocco`.

## Quick Start Demo
Refer to the interactive demonstration: [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb). Includes a walkthrough using publicly available ENCODE ATAC-seq alignments.

## Documentation
See the [Rocco API Reference](https://nolan-h-hamilton.github.io/ROCCO/rocco/rocco.html) for a more in-depth treatment of technical aspects as well as several usage examples, parameter descriptions, etc.

## Testing
[test_rocco.py](https://github.com/nolan-h-hamilton/ROCCO/blob/main/tests/test_rocco.py) is run via GitHub Actions at each release in the [`Tests`](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml) workflow. Simulates test data using  ENCODE ATAC-seq alignments and ensures functionality of fundamental tasks. To run
locally, `cd` in to the `tests` directory and run `pytest test_rocco.py`.

## Citation
```
ROCCO: A Robust Method for Detection of Open Chromatin via Convex Optimization
Nolan H. Hamilton, Terrence S. Furey
bioRxiv 2023.05.24.542132; doi: https://doi.org/10.1101/2023.05.24.542132
```
