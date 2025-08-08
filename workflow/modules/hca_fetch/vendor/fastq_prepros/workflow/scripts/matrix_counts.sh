#!/usr/bin/env bash
set -euo pipefail

FASTQ_R1=$1            # always R1
FASTQ_R2=$2            # always R2
GENOME_DIR=$3
CB_WHITELIST=$4
OUTPUT_PREFIX=$5       # e.g. output/counts/.../{sample}_
OUTPUT_MATRIX=$6       # e.g. output/counts/.../{sample}_matrix.txt
CB_LEN=$7
UMI_LEN=$8
THREADS=$9

if [[ $# -ne 9 ]]; then
  echo "Usage: $0 FASTQ_R1 FASTQ_R2 GENOME_DIR CB_WHITELIST OUTPUT_PREFIX OUTPUT_MATRIX CB_LEN UMI_LEN THREADS" >&2
  exit 1
fi

# 1) Ensure output dir
mkdir -p "$(dirname "$OUTPUT_MATRIX")"

# 2) Detect R1 length and decide barcode check
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

# 3) Run STARsolo—R1 then R2, guaranteed
STAR --runThreadN "$THREADS" \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$FASTQ_R2" "$FASTQ_R1" \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist "$CB_WHITELIST" \
     --soloFeatures Gene \
     --soloCBlen "$CB_LEN" \
     --soloUMIlen "$UMI_LEN" \
     --soloBarcodeReadLength "$BARCODE_LEN" \
     --soloUMIfiltering MultiGeneUMI \
     --outFileNamePrefix "${OUTPUT_PREFIX}"

# 4) Find and move the matrix file
SOLO_DIR="${OUTPUT_PREFIX}Solo.out"
if [[ -f "${SOLO_DIR}/Gene/matrix.txt" ]]; then
  mv "${SOLO_DIR}/Gene/matrix.txt" "$OUTPUT_MATRIX"
elif [[ -f "${SOLO_DIR}/Gene/raw/matrix.mtx" ]]; then
  mv "${SOLO_DIR}/Gene/raw/matrix.mtx" "$OUTPUT_MATRIX"
else
  echo "Error: could not find STAR Solo matrix under ${SOLO_DIR}/Gene" >&2
  exit 1
fi
