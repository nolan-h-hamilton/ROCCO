# Environment Details
ROCCO uses a standard bioinformatics computing environment plus a few Python packages for optimization.

## Conda
If you have ```conda``` (https://docs.conda.io/en/latest/) installed, you can create a new environment for ROCCO by running:
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
this will create a conda environment `rocco` with dependencies that can be loaded by running `conda activate rocco`.


## Manual Installation

```pip install rocco```

### Dependencies
**samtools**

samtools (https://www.htslib.org/download/) v1.17


***bedtools***

bedtools (https://bedtools.readthedocs.io/en/latest/) v2.30.0 


***ucsctools***

https://hgdownload.soe.ucsc.edu/admin/exe/

```
bigWigToWig
bigWigCat
wigToBigWig
```
