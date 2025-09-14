# inputs

PROJ="rs"
SAMPLE="pdl57"
BAM="../../../results/te/"$PROJ"/"$SAMPLE"/star_alignment/Aligned.sortedByCoord.out.bam"
WL="../../../results/te/"$PROJ"/"$SAMPLE"/star_alignment/Solo.out/Gene/filtered/barcodes.tsv"

# output
OUT="../../../results/te/"$PROJ"/"$SAMPLE"/cellsort/Aligned.sortedByCB.bam"

# go
mkdir -p "$(dirname "$OUT")"
stellarscope cellsort \
  --nproc 12 \
  --tempdir /datos/sensence/emilio/tmp \
  --outfile "$OUT" \
  "$BAM" \
  "$WL"
