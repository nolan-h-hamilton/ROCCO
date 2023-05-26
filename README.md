# [R]obust [O]pen [C]hromatin Dection via [C]onvex [O]ptimization
<p align="center">
<img width="472" alt="logo" src="https://github.com/nolan-h-hamilton/ROCCO/assets/57746191/170478f1-5820-4056-b315-3c8dee3603d9">

[![Tests](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml/badge.svg)](https://github.com/nolan-h-hamilton/ROCCO/actions/workflows/tests.yml)
  
ROCCO determines accessible chromatin by utilizing **multiple samples'** enrichment signals *and* spatial features thereof. The model underlying ROCCO is a constrained optimization problem that can be solved efficiently to predict open chromatin regions. 

**Clone/download this repository to use ROCCO.**
  ```
  git clone https://github.com/nolan-h-hamilton/ROCCO.git
  ```

*For a walkthrough using publicly available data with visualized results, see the step-by-step demonstrations in [demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo.ipynb) or [heart_da_demo.ipynb](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo_files/heart_da_demo.ipynb)*

## Environment Details
ROCCO has been developed and tested to run on ATAC-seq alignments in a standard unix bioinformatics environment with Python3.

A ROCCO-specific conda environment can be created using 
[rocco_conda.yml](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/CONDA/rocco_conda.yml):
```
conda env create -n rocco --file docs/CONDA/rocco_conda.yml
```
will create a conda environment `rocco` with necessary dependencies that can be loaded via `conda activate rocco`.

Alternatively, for *manual installation* of dependencies and other general details, refer to [docs/environment.md](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/environment.md).

## Input 
The core ROCCO algorithm takes a replicate-by-locus matrix, $\mathbf{S} \in \mathbb{R}^{K\times n}$, of enrichment signals as input. This data can be generated from replicates' BAM alignments by running

```
python3 prep_bams.py -i <path_to_bam_files>
```

If not already sorted, uses pysam to call `samtools sort` on the BAM files to sort by coordinates. If BAM files have been sorted but indexes are not available, invoke the `--index` parameter. A full list of parameters can be viewed via `python3 prep_bams.py --help`.

The primary function of `prep_bams.py` is [PEPATAC's](https://github.com/databio/pepatac) `bamSitesToWig.py` tool that creates smooth, fixed-step enrichment signals for each replicate.  Alternative methods (e.g., [deepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)) to generate the necessary signal tracks exist but have not been tested. 

See the [flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png) for a visual demonstration of the data preprocessing step.


## Main Scripts
### `ROCCO.py`
`ROCCO.py` calls `ROCCO_chrom.py` for each chromosome specified in a CSV file--see [`params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/params.csv)  or [`demo_params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/demo_files/demo_params.csv) for syntax. For any `NULL` entries in this file, the corresponding genome-wide default, set with the CLI parameters below, will be used.
| **Parameter**                                 	| **Description**                                                                                                                                                                                                                                                                                            	| **Default**                                                                    	|
|-----------------------------------------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|--------------------------------------------------------------------------------	|
| `-p`<br><font size ='1'>`--param_file`</font> 	| File containing chromosome-specific values for $b,\gamma,\tau,c_1,c_2,c_3$ parameters. If any entries are left as `NULL`, ROCCO will replace them with the values specified below. For example, if the budget column is left `NULL`, ROCCO will run using budget `-b` for all chromosomes.                 	| [`params.csv`](https://github.com/nolan-h-hamilton/ROCCO/blob/main/params.csv) 	|
| $b$<br>`-b`<br>`--budget`                     	| (Genome-wide) Budget parameter, i.e., that maximum fraction of chromatin that can be classified as accessible.                                                                                                                                                                                             	| .035                                                                            	|
| $\gamma$<br>`-g`<br>`--gamma`                 	| (Genome-wide) Gamma parameter controlling influence of adjacent loci on selections                                                                                                                                                                                                                         	| 1.0                                                                            	|
| $\tau$<br>`-t`<br>`--tau`                     	| (Genome-wide) Median enrichment threshold for a locus to receive a positive score.                                                                                                                                                                                                                         	| 0.0                                                                            	|
| `--c1`<br>`--c2`<br>`--c3`                    	| (Genome-wide) Weights $c_1,c_2,c_3$ for locus score function $\mathcal{S}(i)$, defined in the paper.                                                                                                                                                                                                       	| 1.0,1.0,1.0                                                                    	|
| `-N, --rr_iter`                                 	| number of `RR` solutions to generate when refining the solution to the LP-relaxation 	| `50`      	|
| `--solver`                                    	| Optimization software used to solve the LP underlying ROCCO. ECOS is used by default and is a viable open-source option. 'MOSEK' offers significantly greater speed and is free for academic use. Free trial commercial licenses are also available. See https://www.mosek.com/products/academic-licenses/ 	| `ECOS`                                                                         	|
| <font size ='2'>`--bed_format`</font>         	| `--bed_format 6` creates output in BED6 format. `--bed_format 3` creates output in BED3 format. The naming convention for peaks in BED6 format is, for example, `chr22_1` for the first peak called in chromosome 22.                                                                                        	| `6`                                                                            	|
| <font size ='2'>`--verbose`</font>            	| If invoked, detailed output will be generated while ROCCO is running                                                                                                                                                                                                                                       	| `False`                                                                        	|
| <font size ='2'>`--outdir`</font>            	| directory to store output files in                                                                                                                                                                                                                                       	| `.`                                                                        	|
| <font size ='2'>`--identifiers`</font>            	| path to a file with sample identifiers on each line. Allows specifying a subset of samples in `--wig_path` rather than using all `.wig` files. Identifiers must be a substring of the respective filepath in `--wig_path`                                                                                                                                                                                                                                     	| `None`                                                                        	|
| <font size ='2'>`--combine`</font>            	| If not `None`, concatenates chromosome-specific output bedfiles into the file specified with this parameter. e.g., `--combine bedname.bed` stores output in a single file `bedname.bed`                                                                                                                                                                                                                                     	| `None`     
| <font size ='2'>`--jobs`</font>            	| Number of ROCCO_chrom.py jobs to run in parallel. ROCCO will process each chromosome in-order if left unspecified. If set to `0`, the maximum possible number of jobs will be created. Use judiciously depending on available memory.                               	| `1`     
| <font size ='2'>`--mem`</font>            	| Optional: if running in parallel, only start a ROCCO_chrom.py job if there is at least `-mem` memory avilable. See GNU Parallel's `--memfree` parameter                                    	| `None`     

**Example**
Run `ROCCO` genome-wide on hg38 with default parameters 
```
python ROCCO.py -p params.csv
```

### `ROCCO_chrom.py`
This script gathers chromosome-specific signal tracks and solves the constrained optimization problem underlying ROCCO to predict open chromatin regions.

*Note*, when computing genome-wide peaks, depending on your computing environment (e.g., SLURM), it may be preferable to create `ROCCO_chrom.py` jobs for each chromosome for simultaneous distribution across multiple nodes rather than running `ROCCO.py`.
| **Parameter**                              	| **Description**                                                                                                                                                                                                                                                                                            	| **Default** 	|
|--------------------------------------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-------------	|
| $b$<br>`-b`<br>`--budget`                  	| Budget parameter, i.e., that maximum fraction of chromatin that can be classified as accessible.                                                                                                                                                                                                           	| .035         	|
| $\gamma$<br>`-g`<br>`--gamma`              	| Gamma parameter controlling influence of adjacent loci on selections                                                                                                                                                                                                                                       	| 0.0         	|
| $\tau$<br>`-t`<br>`--tau`                  	| Median enrichment threshold for a locus to receive a positive score.                                                                                                                                                                                                                                       	| 1.0         	|
| `--c1`<br>`--c2`<br>`--c3`<br>             	| Weights $c_1,c_2,c_3$ for locus score function $\mathcal{S}(i)$, defined in the paper.                                                                                                                                                                                                                     	| 1.0,1.0,1.0 	|
| <font size ='2'>`--wig_path`</font>        	| The directory containing chromosome-specific signal tracks for each replicate.                                                                                                                                                                                                                             	| `cwd`       	|
| `--solver`                                 	| Optimization software used to solve the LP underlying ROCCO. ECOS is used by default and is a viable open-source option. 'MOSEK' offers significantly greater speed and is free for academic use. Free trial commercial licenses are also available. See https://www.mosek.com/products/academic-licenses/ 	| `ECOS`      	|
| `-N, --rr_iter`                                 	| number of `RR` solutions to generate when refining the solution to the LP-relaxation 	| `50`      	|
| <font size='2'>`--bed_format`</font>       	| `--bed_format 6` creates output in BED6 format. `--bed_format 3` creates output in BED3 format. The naming convention for peaks in BED6 format is, for example, `chr22_1` for the first peak called in chromosome 22.                                                                                        	| `6`         	|
| <font size ='2'>`--verbose`</font>         	| If invoked, detailed output will be generated while ROCCO is running                                                                                                                                                                                                                                       	| `False`     	|
| <font size ='2'>`--outdir`</font>            	| directory to store output file in                                                                                                                                                                                                                                       	| `.`                                                                        	|
| <font size ='2'>`--identifiers`</font>            	| path to a file with sample identifiers on each line. Allows specifying a subset of samples in `--wig_path` rather than using all `.wig` files. Identifiers must be a substring of the respective filepath in `--wig_path`                                                                                                                                                                                                                                     	| `None`                                                                        	|

**Example**
Run ROCCO on a single chromosome, chr22, with default parameters using all samples in `tracks_chr22`. Store results in current working directory.
```
python3 ROCCO_chrom.py --wig_path tracks_chr22 --chrom chr22 
```
## Auxiliary Tools
ROCCO-specific wrappers for common bioinformatics tasks included for convenience.

### `count_matrix.py`
Given a ROCCO peak file in BED6 format, generates a count matrix, with each row corresponding to a peak and each column to an alignment file from which read-counts are derived.

Use the `-i` parameter to designate the peak file. `--bamfiles` specifies a file containing a row for each BAM filepath.

The default output file is `ROCCO_peak_counts.tsv`

**Example**
Assuming BED6 peakfile `ROCCO_combined.bed` and BAM files are listed in `sample.data`
```
python3 count_matrix.py -i ROCCO_combined.bed --bamfiles sample.data
```

### `prep_bams.py`
Processes BAM alignments into ROCCO conformable input. This preprocessing script utilizes [PEPATAC's](https://github.com/databio/pepatac) `bamSitesToWig.py` tool to create smooth, fixed-step enrichment signals for each replicate. A [flowchart](https://github.com/nolan-h-hamilton/ROCCO/blob/main/docs/bamsig_flowchart.png) shows the pipeline carried out by this script.

If not done already, uses pysam to call `samtools sort` on BAM files to coordinate-sort them. If BAM files have been sorted but index files are not available, invoke the `--index` parameter.  *A full list of parameters can be viewed via `python3 prep_bams.py --help`*

**Example**

```
python3 prep_bams.py -i <path_to_bams, def='.'> -L <step-size, def=50> --cores <num. cores, def=1> 
```
