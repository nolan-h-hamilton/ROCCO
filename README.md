# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)

Underlying ROCCO is a constrained optimization problem that can be solved efficiently to **predict consensus regions of open chromatin**.

**Features**

1. Explicitly accounts for both **enrichment and spatial characteristics** of open chromatin signals, the latter of which is an informative but often ignored aspect of ATAC-seq data;
1. Leverages data from multiple samples of varying quality **without imposing arbitrary thresholds** on a minimum number of samples declaring peaks;
1. Is efficient for **large numbers of samples** with an asymptotic time complexity independent of sample/replicate count;
1. **Does not require training data or a heuristically determined set of initial candidate regions**, which are hard to define given the lack of a priori sets of open chromatin regions;
1. Employs a **mathematically tractable model** permitting guarantees of performance and efficiency.


## Environment Details
ROCCO has been developed and tested to run on ATAC-seq alignments in a standard unix bioinformatics environment with Python3.7+

A ROCCO-specific conda environment can be created using 
[rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/CONDA/rocco_conda.yml):
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
will create a conda environment `rocco` with necessary dependencies that can be loaded via `conda activate rocco`.

Alternatively, for *manual installation* of dependencies and other general details, refer to [docs/environment.md](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/environment.md).

## Input 
ROCCO accepts a set of enrichment signal tracks used to generate a sample-by-locus matrix, $\mathbf{S} \in \mathbb{R}^{K\times n}$, for scoring and optimization.
  
These signal tracks are generated and formatted appropriately from samples' BAM files by running
```
python3 prep_bams.py -i <directory_path_containing_bam_files> --cores <num_cores>
```
Run `python3 prep_bams.py --help` for a full list of parameters. 

See a [flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png) for a visual demonstration of the data preprocessing step. Readers may also reference the [tests](https://github.com/nolan-h-hamilton/ROCCO/blob/main/tests) directory which contains toy data used in the testing workflow that is structured conformably for ROCCO.

## Getting Started
Clone/download this repository to use ROCCO
  ```
  git clone https://github.com/nolan-h-hamilton/ROCCO.git
  ```
### Jupyter Notebook Demos
  Some lightweight demonstrations using publicly available data.
  1. **Basic Usage**. [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo.ipynb). Includes visualized results in UCSC genoome browser.
  1. **Differential Accessibility**. [heart_da_demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo_files/heart_da_demo.ipynb). A sample ROCCO-->[DESeq2](https://github.com/mikelove/DESeq2) pipeline for differential analysis.

  
## Main Scripts
### `ROCCO.py`
`ROCCO.py` calls `ROCCO_chrom.py` for each chromosome specified in a CSV file--see [`params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/params.csv)  or [`demo_params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo_files/demo_params.csv) for syntax. *For any `NULL` entries in this file, the corresponding genome-wide parameter value, set with the CLI arguments in the table below, will be used*. The default values for these parameters yield good general performance, but users seeking more/less conservative or topologically simple predictions can modify as necessary. The budget parameter for each chromosome can be inferred as a function
of read density---see `est_budgets.py`.
| **Parameter**                                 	| **Description**                                                                                                                                                                                                                                                                                            	| **Default**                                                                    	|
|-----------------------------------------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|--------------------------------------------------------------------------------	|
| `-p`<br><font size ='1'>`--param_file`</font> 	| File containing chromosome-specific values for $b,\gamma,\tau,c_1,c_2,c_3$ parameters. If any entries are left as `NULL`, ROCCO will replace them with the values specified below. For example, if the budget column is left `NULL`, ROCCO will run using budget `-b` for all chromosomes.                 	| [`params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/params.csv) 	|
| $b$<br>`-b`<br>`--budget`                     	| Budget parameter, i.e., that maximum fraction of chromatin that can be classified as accessible in each chromosome.                                                                                                                                                                                             	| .035                                                                            	|
| $\gamma$<br>`-g`<br>`--gamma`                 	|  Gamma parameter controlling influence of adjacent loci on selections                                                                                                                                                                                                                         	| 1.0                                                                            	|
| $\tau$<br>`-t`<br>`--tau`                     	|  Median enrichment threshold for a locus to receive a positive score.                                                                                                                                                                                                                         	| 0.0                                                                            	|
| `--c1`<br>`--c2`<br>`--c3`                    	| Weights $c_1,c_2,c_3$ for locus score function $\mathcal{S}(i)$, defined in the paper. Note, modifying $c_1 < c_2$ can cause negative locus scores which may cause unexpected results.                                                                                                                                                                                                     	| 1.0,1.0,1.0                                                                    	|
| `-N, --rr_iter`                                 	| number of `RR` solutions to generate when refining the solution to the LP-relaxation 	| `50`      	|
| `--solver`                                    	| Optimization software used to solve the LP underlying ROCCO. ECOS is used by default and is a viable open-source option. `MOSEK` offers significantly greater speed and is free for academic use. Free trial commercial licenses are also available. See https://www.mosek.com/products/academic-licenses/ 	| `ECOS`                                                                         	|
| <font size ='2'>`--bed_format`</font>         	| `--bed_format 6` creates output in BED6 format. `--bed_format 3` creates output in BED3 format. The naming convention for peaks in BED6 format is, for example, `chr22_1` for the first peak called in chromosome 22.                                                                                        	| `6`                                                                            	|
| <font size ='2'>`--verbose`</font>            	| If invoked, detailed output will be generated while ROCCO is running                                                                                                                                                                                                                                       	| `False`                                                                        	|
| <font size ='2'>`--outdir`</font>            	| directory to store output files in                                                                                                                                                                                                                                       	| `.`                                                                        	|
| <font size ='2'>`--identifiers`</font>            	| path to a file with sample identifiers on each line. Allows specifying a subset of samples in `--wig_path` rather than using all `.wig` files. Identifiers must be a substring of the respective file in `--wig_path`                                                                                                                                                                                                                                     	| `None`                                                                        	|
| <font size ='2'>`--combine`</font>            	| If not `None`, concatenates chromosome-specific output bedfiles into the file specified with this parameter. e.g., `--combine bedname.bed` stores output in a single file `bedname.bed`                                                                                                                                                                                                                                     	| `None`     
| <font size ='2'>`--jobs`</font>            	| Number of ROCCO_chrom.py jobs to run in parallel. ROCCO will process each chromosome in-order if left unspecified. If set to `0`, the maximum possible number of jobs will be created. Use judiciously depending on available memory.                               	| `1`     
| <font size ='2'>`--mem`</font>            	| Optional: if running in parallel, only start a ROCCO_chrom.py job if there is at least `-mem` memory avilable. See GNU Parallel's `--memfree` parameter                                    	| `None`     


### `ROCCO_chrom.py`
This script gathers chromosome-specific signal tracks and solves the constrained optimization problem underlying ROCCO to predict open chromatin regions.

*Note*, when computing genome-wide peaks, depending on your computing environment (e.g., SLURM), it may be preferable to create `ROCCO_chrom.py` jobs for each chromosome for simultaneous distribution across multiple nodes rather than running `ROCCO.py`.
| **Parameter**                              	| **Description**                                                                                                                                                                                                                                                                                            	| **Default** 	|
|--------------------------------------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-------------	|
| $b$<br>`-b`<br>`--budget`                  	| Budget parameter, i.e., that maximum fraction of chromatin that can be classified as accessible.                                                                                                                                                                                                           	| .035         	|
| $\gamma$<br>`-g`<br>`--gamma`              	| Gamma parameter controlling influence of adjacent loci on selections                                                                                                                                                                                                                                       	| 0.0         	|
| $\tau$<br>`-t`<br>`--tau`                  	| Median enrichment threshold for a locus to receive a positive score.                                                                                                                                                                                                                                       	| 1.0         	|
| `--c1`<br>`--c2`<br>`--c3`<br>             	| Weights $c_1,c_2,c_3$ for locus score function $\mathcal{S}(i)$ defined in the paper. Note, modifying $c_1 < c_2$ can cause negative locus scores which may cause unexpected results.                                                                                                                                                                                                                     	| 1.0,1.0,1.0 	|
| <font size ='2'>`--wig_path`</font>        	| The directory containing chromosome-specific signal tracks for each replicate.                                                                                                                                                                                                                             	| `cwd`       	|
| `--solver`                                 	| Optimization software used to solve the LP underlying ROCCO. ECOS is used by default and is a viable open-source option. `MOSEK` offers significantly greater speed and is free for academic use. Free trial commercial licenses are also available. See https://www.mosek.com/products/academic-licenses/ 	| `ECOS`      	|
| `-N, --rr_iter`                                 	| number of `RR` solutions to generate when refining the solution to the LP-relaxation 	| `50`      	|
| <font size='2'>`--bed_format`</font>       	| `--bed_format 6` creates output in BED6 format. `--bed_format 3` creates output in BED3 format. The naming convention for peaks in BED6 format is, for example, `chr22_1` for the first peak called in chromosome 22.                                                                                        	| `6`         	|
| <font size ='2'>`--verbose`</font>         	| If invoked, detailed output will be generated while ROCCO is running                                                                                                                                                                                                                                       	| `False`     	|
| <font size ='2'>`--outdir`</font>            	| directory to store output file in                                                                                                                                                                                                                                       	| `.`                                                                        	|
| <font size ='2'>`--identifiers`</font>            	| path to a file with sample identifiers on each line. Allows specifying a subset of samples in `--wig_path` rather than using all `.wig` files. Identifiers must be a substring of the respective file in `--wig_path`                                                                                                                                                                                                                                     	| `None`                                                                        	|
## Auxiliary Scripts
Used for pre/post-processing. Use cases included in the demo notebooks and/or the scripts' module docstrings.
1. `prep_bams.py`: used to generate smooth, fixed-step signal tracks for each chromosome/replicate
1. `est_budgets.py`: used to determine chromosome-specific budgets as a function of read density
1. `count_matrix.py`: used to create a count matrix for ROCCO peak results