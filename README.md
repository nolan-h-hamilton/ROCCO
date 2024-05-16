# ROCCO: [R]obust [O]pen [C]hromatin Detection via [C]onvex [O]ptimization

<p align="center">
<img width="400" alt="logo" src="docs/logo.png">

ROCCO is an optimal consensus peak calling algorithm for chromatin accessibility data that is scalable to large sample sizes.

**Features**

1. **Consideration of enrichment and spatial characteristics** of open chromatin signals to capture the full extent of peaks
2. **Mathematically tractable model** that permits performance and efficiency guarantees
3. **Efficient for large numbers of samples** with an asymptotic time complexity independent of sample size
4. **No arbitrary thresholds** on the minimum number of supporting samples/replicates
5. **No required training data** or a heuristically determined set of initial candidate peak regions


# Demo

A brief walkthrough with visualized peak results using publicly available ATAC-seq data:

https://github.com/nolan-h-hamilton/ROCCO/tree/main/docs/demo/demo.ipynb

# Paper/Citation

If using ROCCO in your research, please cite the [original paper](https://doi.org/10.1093/bioinformatics/btad725) in *Bioinformatics*


   ```
    Nolan H Hamilton, Terrence S Furey, ROCCO: a robust method for detection of open chromatin via convex optimization,
    Bioinformatics, Volume 39, Issue 12, December 2023
   ```

**DOI**: [10.1093/bioinformatics/btad725](`10.1093/bioinformatics/btad725`)

# Documentation

ROCCO's documentation is available at https://nolan-h-hamilton.github.io/ROCCO/

# Installation

   ```
   pip install rocco
   ```

ROCCO utilizes the popular bioinformatics software [Samtools](http://www.htslib.org) and [bedtools](https://bedtools.readthedocs.io/en/latest/). If not available already, these system dependencies can be installed with standard MacOS or Linux/Unix package managers, e.g., `brew install samtools` (Homebrew), `sudo apt-get install samtools`) (APT).


# Input

To determine consensus peaks (BED format), ROCCO accepts **BAM** alignments and a genome sizes file as input

   ```
   rocco -i sample1.bam sample2.bam sample3.bam [...] -g hg38.sizes --params hg38
   ```

or with a wildcard:

   ```
   rocco -i *.bam -g hg38.sizes --params hg38
   ```

See `rocco --help` for more details.

