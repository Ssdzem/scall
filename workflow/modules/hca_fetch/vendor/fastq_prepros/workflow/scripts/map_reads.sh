#!/usr/bin/env bash

set -euo pipefail

# Usage: map_reads.sh <FASTQ_R1> <FASTQ_R2> <GENOME_DIR> <WHITELIST> <OUTPUT_BAM> <THREADS>
if [[ $# -ne 8 ]]; then
  echo "Usage: $0 FASTQ_R1 FASTQ_R2 GENOME_DIR WHITELIST OUTPUT_BAM THREADS" >&2
  exit 1
fi

FASTQ_R1="$1"
FASTQ_R2="$2"
GENOME_DIR="$3"
WHITELIST="$4"
OUTPUT_BAM="$5"
CB_LEN="$6"
UMI_LEN="$7"
THREADS="$8"

# Derive prefix for STAR outputs
OUTPUT_PREFIX="${OUTPUT_BAM%.bam}"

# Detect R1 length and decide barcode check
# (disable pipefail around this pipeline to avoid SIGPIPE abort)
set +o pipefail
read1_seq=$(zcat "$FASTQ_R1" | sed -n '2{p;q;}')
set -o pipefail

R1_LEN=${#read1_seq}
if [[ $((CB_LEN+UMI_LEN)) -eq $R1_LEN ]]; then
  BARCODE_LEN=$((CB_LEN+UMI_LEN))
else
  echo "Note: R1 length ($R1_LEN) != CB+UMI ($((CB_LEN+UMI_LEN))); disabling length check" >&2
  BARCODE_LEN=0
fi

# Ensure output folder exists
mkdir -p "$(dirname "$OUTPUT_BAM")"

echo "Running STAR mapping..." >&2
echo "[$(date)] map_reads.sh starting" >&2
echo "STAR -> $(command -v STAR)" >&2
STAR --version >&2
echo "Running STAR mapping…" >&2

# ---- STAR Solo mapping ----
STAR --runThreadN "$THREADS" \
     --readFilesIn "$FASTQ_R2" "$FASTQ_R1" \
     --readFilesCommand zcat \
     --genomeDir "$GENOME_DIR" \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist "$WHITELIST" \
     --soloCBstart 1 --soloCBlen "$CB_LEN" \
     --soloUMIstart $((CB_LEN+1)) --soloUMIlen "$UMI_LEN" \
     --soloBarcodeReadLength "$BARCODE_LEN"

# Move the BAM into place
mv "${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam" "$OUTPUT_BAM"
