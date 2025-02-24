# Example Use

* For the sake of brevity and to accommodate personal computing environments, we restrict this lightweight demo to $N=3$ alignments over chr16 and chr21.
  * N.B. ROCCO generally offers its greatest performance relative to benchmarks for larger sample sizes, e.g., $N \geq 10$.

## Metadata

Download, index (`samtools index <bamfile>`), and store the following three sequence alignments in a directory `als_bamfiles`

| File (Link) | Sex  | Age | Cell Type     | Condition | Reads (f2/F3840) |
|-------------|------|-----|---------------|-----------|---------------|
| [ENCFF333UUT](https://www.encodeproject.org/files/ENCFF333UUT) | Male | 49  | Motor Neuron | ALS       | 144144508 |
| [ENCFF786BZU](https://www.encodeproject.org/files/ENCFF786BZU) | Male | 47  | Motor Neuron | ALS       | 123872478 |
| [ENCFF742KHO](https://www.encodeproject.org/files/ENCFF742KHO) | Male | 55  | Motor Neuron | ALS       | 216811080 |

```bash
rocco -i als_bamfiles/*.bam -g hg38 -o rocco_peaks_als_n3.bed  --narrowPeak \
--chroms chr16 chr21 --threads 4
```

* Note, the number of threads will affect memory usage during the optimization step, and we have observed diminishing returns with respect to runtime for $\textsf{--threads} > 8$.
* $\textsf{--narrowPeak}$ invokes generation of narrowPeak output, with bootstrapped $p$-values and BH-corrected $q$-values in $-\log_{10}(\cdot)$ scale.
  * This flag also creates a 'raw' peak-by-count matrix (one row for each peak, one column for each sample) that can be used downstream.

### Output Snippets

* `rocco_peaks_als_n3.narrowPeak`

  ```plain
  chr21	31629800	31630000	chr21_31629800_31630000	286	.	37.2162	1.3188	1.2433	31629900
  chr21	31658050	31658350	chr21_31658050_31658350	306	.	57.8229	1.5133	1.2722	31658200
  chr21	31659200	31660700	chr21_31659200_31660700	908	.	668.9143	2.574	1.3943	31659950
  chr21	31665500	31665800	chr21_31665500_31665800	298	.	49.3102	1.3699	1.2566	31665650
  chr21	31730800	31733100	chr21_31730800_31733100	943	.	704.5196	2.8751	1.4623	31731950
  ```

* `rocco_peaks_als_n3.counts.tsv`

    ```plain
    peak_name	ENCFF333UUT	ENCFF742KHO	ENCFF786BZU
    chr21_31629800_31630000	53	106	7
    chr21_31658050_31658350	152	312	11
    chr21_31659200_31660700	3946	6249	4024
    chr21_31665500_31665800	115	92	32
    chr21_31730800_31733100	7009	9298	6408
    ```

### Visualization

50kb region in chromosome 21 centered over $\textsf{SOD1}$ (linked to familial ALS)

![igv_example](als_n3_igv.png)

## Consenrich --> ROCCO

Consensus peak calling is a natural downstream use-case for the multi-sample signal extraction method, [Consenrich](https://github.com/nolan-h-hamilton/Consenrich). See [here](https://github.com/nolan-h-hamilton/Consenrich/blob/main/docs/consensus_peaks.md) for several examples of integration of Consenrich results with ROCCO.
