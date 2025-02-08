import os 
import rocco
# import matplotlib.pyplot as plt
# import seaborn as sns # pip install seaborn

BAM_LIST_FILE = 'bam_files.txt' # text file with paths to bam files (one per line)
PEAK_FILE = 'test_peaks.bed' # ROCCO BED output 
CHROM_SIZES_FILE = 'hg38.chrom.sizes' # sizes file
PEAK_COUNT_MATRIX_FILE = 'peak_count_matrix.tsv' # path to count matrix--if it doesn't exist already, it'll be generated
OUTPUT_SCORED_PEAKS_FILE = 'rocco_peaks_scored.bed' # output scored peaks

# I've not yet tested this with a large peak file. `peak_count_matrix.sh' may be faster
print('\nGenerating raw count matrix...')
if not os.path.exists(PEAK_COUNT_MATRIX_FILE):
    PEAK_COUNT_MATRIX_FILE = rocco.readtracks.raw_count_matrix(bam_list_file=BAM_LIST_FILE,
                                                                        peak_file=PEAK_FILE,
                                                                        output_file=PEAK_COUNT_MATRIX_FILE)
    print('done.\n')

# just checking the above worked
if not os.path.exists(PEAK_COUNT_MATRIX_FILE):
    raise ValueError(f"Failed to find or create count matrix: {PEAK_COUNT_MATRIX_FILE}")


# score peaks and store output in narrowPeak-like format
# note that p-values and q-values are, per MACS convention, in -log10 scale
# likewise, the bed6 scores are scaled in [0,1000] for browser compliance
print(f'\nScoring peaks in {PEAK_COUNT_MATRIX_FILE}...')
peak_scores, bed6_scores, pvals = rocco.readtracks.score_peaks(bam_files=BAM_LIST_FILE,
                                           chrom_sizes_file=CHROM_SIZES_FILE,
                                           peak_file=PEAK_FILE,
                                           count_matrix_file=PEAK_COUNT_MATRIX_FILE,
                                           output_file=OUTPUT_SCORED_PEAKS_FILE)
print(f'Done. Output --> {OUTPUT_SCORED_PEAKS_FILE} \n')


# for reference, you can plot narrowPeak-like scores, signal values, pvals
# uncomment `import seaborn as sns` and the code below

# sns.histplot(peak_scores, kde=True, color='blue', label='Regular Scale')
# plt.savefig('peak_scores_regular.png',dpi=300)
# plt.close()

# sns.histplot(bed6_scores, kde=True, color='red', label='Bed6 Scale')
# plt.savefig('bed6_scores_regular.png',dpi=300)
# plt.close()


# optionally, if you want to ensure the p-values are calibrated, supply a 'null' peak file 
# above in a separate run and checking they follow a uniform distribution: 
# ...e.g., multiple calls`bedtools random -n 100 -l [50,100,250,500,750,1000,...] >> null.bed`
# For 'real' peaks, the distribution will be left-skewed

# sns.histplot(pvals, kde=True, color='green', label='Pvals')
# plt.savefig('pvals_regular.png',dpi=300)
# plt.close()