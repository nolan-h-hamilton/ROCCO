# Environment Details
ROCCO uses a standard bioinformatics computing environment plus a few Python packages for optimization.

## Conda
If you have ```conda``` (https://docs.conda.io/en/latest/) installed, you can create a new environment for ROCCO by running:
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
this will create a conda environment `rocco` with dependencies that can be loaded by running `conda activate rocco`.


## Manual Installation

```pip install rocco --upgrade```

### Dependencies
The dependencies can be obtained through [Bioconda](https://bioconda.github.io) or downloaded at the links listed below.

**samtools**

samtools (https://www.htslib.org/download/)

Bioconda: `samtools`

***bedtools***

bedtools (https://bedtools.readthedocs.io/en/latest/)

Bioconda: `bedtools`

***UCSC Binary Utilities***

https://hgdownload.soe.ucsc.edu/admin/exe/

```
bigWigToWig
bigWigCat
wigToBigWig
```

Bioconda:
```
ucsc-bigwigtowig
ucsc-wigtobigwig
ucsc-bigwigcat
```