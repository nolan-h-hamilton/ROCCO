#!/bin/bash
#
#   <bam_list_file>   A file listing one BAM filepath per line
#   <bed_file>        rocco peak BED output
#   <output_file>     count matrix output

if [ $# -ne 3 ]; then
  echo "format: $0 <bam_list_file> <bed_file> <output_file>" >&2
  exit 1
fi

BAM_LIST_FILE=$1
BED_FILE=$2
OUTPUT_FILE=$3

BAMS=""
HEADER="peak_name"

while IFS= read -r bam
do
    sample_name=$(basename "$bam")
    sample_name=$(echo -e "$sample_name" | sed 's/\.bam$//')
    BAMS="$BAMS $bam"
    HEADER="$HEADER\t$sample_name"
done < "$BAM_LIST_FILE"


bedtools multicov -bams $BAMS -bed "$BED_FILE" | \
awk -v header="$HEADER" '
  BEGIN {
    print header
  }
  {
    # joins bed entries into a single column `peak_name` -- hopefully more convenient for DESeq2
    peak_name = $1"_"$2"_"$3
    printf "%s", peak_name
    for(i=4; i<=NF; i++){
      printf "\t%s", $i
    }
    printf "\n"
  }
' > "$OUTPUT_FILE"