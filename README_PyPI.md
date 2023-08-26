# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
1. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
1. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample count;
1. **Does not require training data** or initial candidate peak regions which are hard to define given the lack of a priori sets of open chromatin regions;
1. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.


#### [Quick Start Demo](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb)
The quick start demo is an interactive Jupyter Notebook showcasing ROCCO's functionality.

## Documentation
API Reference: https://nolan-h-hamilton.github.io/ROCCO/index.html

## Citation
```
ROCCO: A Robust Method for Detection of Open Chromatin via Convex Optimization
Nolan H. Hamilton, Terrence S. Furey
bioRxiv 2023.05.24.542132; doi: https://doi.org/10.1101/2023.05.24.542132
```
