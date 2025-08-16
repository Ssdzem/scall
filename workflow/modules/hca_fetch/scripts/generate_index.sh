#!/bin/bash

set -euo pipefail
# Arguments
GENOME_FASTA=$1
GTF_FILE=$2
WHITELIST=$3
OUTPUT_DIR=$4
THREADS=$5

# (Opcional) Verificar que exista la whitelist
if [[ ! -f "$WHITELIST" ]]; then
  echo "Error: whitelist file $WHITELIST not found!" >&2
  exit 1
fi

# STAR index generation usando THREADS
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$OUTPUT_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 12 \
     --soloCBwhitelist "$WHITELIST"
