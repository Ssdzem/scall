#!/usr/bin/env bash
set -euo pipefail

# Usage: map_reads.sh FASTQ_R1 FASTQ_R2 GENOME_DIR WHITELIST OUTPUT_BAM CB_LEN UMI_LEN THREADS
if [[ $# -ne 8 ]]; then
  echo "Usage: $0 FASTQ_R1 FASTQ_R2 GENOME_DIR WHITELIST OUTPUT_BAM CB_LEN UMI_LEN THREADS" >&2
  exit 2
fi

FASTQ_R1="$1"
FASTQ_R2="$2"
GENOME_DIR="$3"
WHITELIST="$4"
OUTPUT_BAM="$5"
CB_LEN="$6"
UMI_LEN="$7"
THREADS="$8"

# Output prefix: same folder as OUTPUT_BAM, file stem + underscore (STAR appends filenames)
outdir="$(dirname "$OUTPUT_BAM")"
prefix="${outdir}/$(basename "${OUTPUT_BAM%.bam}")_"
mkdir -p "$outdir"

echo "[$(date)] map_reads.sh starting" >&2
echo "STAR at: $(command -v STAR)" >&2
STAR --version >&2

# Important details:
# - For 10x, STARsolo expects cDNA read FIRST, barcode read SECOND -> R2 then R1.  (docs below)
# - Disable strict barcode-read-length check to tolerate 26/28 bp mixes: --soloBarcodeReadLength 0
# - CB/UMI positions: CB starts at 1 (1-based), UMI starts right after CB.
STAR \
  --runThreadN "$THREADS" \
  --readFilesIn "$FASTQ_R2" "$FASTQ_R1" \
  --readFilesCommand zcat \
  --genomeDir "$GENOME_DIR" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$prefix" \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WHITELIST" \
  --soloCBstart 1 --soloCBlen "$CB_LEN" \
  --soloUMIstart "$((CB_LEN+1))" --soloUMIlen "$UMI_LEN" \
  --soloBarcodeReadLength 0

# Move coordinate-sorted BAM to the requested path
mv "${prefix}Aligned.sortedByCoord.out.bam" "$OUTPUT_BAM"
