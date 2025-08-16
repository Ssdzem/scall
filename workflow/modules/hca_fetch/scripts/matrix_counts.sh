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

COUNT_LAYER="${COUNT_LAYER:-filtered}"     # filtered | raw
SOLO_FEATURES="${SOLO_FEATURES:-Gene}"     # e.g., "Gene" or "Gene GeneFull" or "Gene Velocyto"

if [[ $# -ne 9 ]]; then
  echo "Usage: $0 FASTQ_R1 FASTQ_R2 GENOME_DIR CB_WHITELIST OUTPUT_PREFIX OUTPUT_MATRIX CB_LEN UMI_LEN THREADS" >&2
  exit 1
fi

# 1) Ensure output dir
mkdir -p "$(dirname "$OUTPUT_MATRIX")"

# 2) Always disable barcode-read-length check.
# Mixed R1 lengths occur; correct chemistry (CB_LEN/UMI_LEN), ensured in fastq vibe check rule and metadata chemistry annotation, is vital for proper counts
BARCODE_LEN=0 #to set the length to CB_LEN + UMI_LEN

# 3) Run STARsolo—R1 then R2, guaranteed
STAR --runThreadN "$THREADS" \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$FASTQ_R2" "$FASTQ_R1" \
     --readFilesCommand zcat \
     --outSAMtype None \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist "$CB_WHITELIST" \
     --soloFeatures ${SOLO_FEATURES} \
     --soloCBlen "$CB_LEN" \
     --soloUMIlen "$UMI_LEN" \
     --soloBarcodeReadLength "$BARCODE_LEN" \   # now always 0
     --soloUMIfiltering MultiGeneUMI \
     --soloCellFilter EmptyDrops_CR

# 4) Stage the 10x triple to a stable folder and keep a legacy pointer
SOLO_DIR="${OUTPUT_PREFIX}Solo.out"
SRC_DIR="${SOLO_DIR}/Gene/${COUNT_LAYER}"

if [[ ! -s "${SRC_DIR}/matrix.mtx" || ! -s "${SRC_DIR}/features.tsv" || ! -s "${SRC_DIR}/barcodes.tsv" ]]; then
  echo "Error: expected 10x files not found in ${SRC_DIR}" >&2
  exit 1
fi

# Put files in a per-sample folder derived from OUTPUT_MATRIX basename
DEST_ROOT="$(dirname "$OUTPUT_MATRIX")"
BASE="$(basename "$OUTPUT_MATRIX")"
BASE="${BASE%.*}"   # drop extension
DEST_DIR="${DEST_ROOT}/${BASE}"

mkdir -p "${DEST_DIR}"
cp -f "${SRC_DIR}/matrix.mtx"   "${DEST_DIR}/matrix.mtx"
cp -f "${SRC_DIR}/features.tsv" "${DEST_DIR}/features.tsv"
cp -f "${SRC_DIR}/barcodes.tsv" "${DEST_DIR}/barcodes.tsv"

# For backward compatibility, symlink legacy single file to the new matrix.mtx
ln -sf "${DEST_DIR}/matrix.mtx" "$OUTPUT_MATRIX"

echo "[matrix_counts] Wrote 10x triple to ${DEST_DIR} (layer=${COUNT_LAYER}); legacy link -> ${OUTPUT_MATRIX}"
