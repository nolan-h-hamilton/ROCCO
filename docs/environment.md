# Environment Details
ROCCO uses a standard bioinformatics computing environment plus a few Python packages for optimization.

## Conda
If you have ```conda``` (https://docs.conda.io/en/latest/) installed, you can create a new environment for ROCCO by running:
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
this will create a conda environment `rocco` with all necessary dependencies that can be loaded by running `conda activate rocco` at the command line. This environment was assembled and tested using Conda *v23.3.1*

This should install all dependencies for both ROCCO and the ```prep_bams.py``` script that performs the BAM-->WIG pipeline.

## Manual Installation
### Running ROCCO
**Version**

Python 3.7+

**Packages**
```
cvxpy
pandas
numpy
scipy
```

All above packages can be easily installed by running:
```
pip install <package-name>
```

### Running the BAM --> WIG Pipeline
```prep_bams.py``` uses the following packages to convert replicates' BAM files into ROCCO conformable input via PEPATAC tools. Use if you do not already have fixed-step wiggle tracks for each chromosome/replicate.

**Additional Python Packages**

If running this pipeline, you should install these Python packages as well.
```
logmuse
pararead
pysam
pybedtools
```

All above packages can be easily installed by running:
```
pip install <package-name>
```


**samtools**

samtools (https://www.htslib.org/download/) v1.17


***bedtools***

bedtools (https://bedtools.readthedocs.io/en/latest/) v2.30.0 


***ucsctools***

```prep_bams.py``` uses several UCSC binary tools
available at https://hgdownload.soe.ucsc.edu/admin/exe/

```
bigWigToWig
bigWigCat
wigToBigWig
```
