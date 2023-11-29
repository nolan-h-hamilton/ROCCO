# [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/4e662eda-4899-46fa-ac9f-998e28f592c1">

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin** across multiple samples.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks;
2. **No arbitrary thresholds** on the minimum number of supporting samples/replicates;
3. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample count;
4. **Does not require training data** or initial candidate peak regions which are hard to define given the lack of a priori sets of open chromatin regions;
5. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.


#### [Quick Start Demo](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo/demo.ipynb)
The quick start demo is an interactive Jupyter Notebook showcasing ROCCO's functionality.

## Documentation
API Reference: https://nolan-h-hamilton.github.io/ROCCO/index.html

## Citation
If using ROCCO in your research, please cite the [corresponding paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics*.

