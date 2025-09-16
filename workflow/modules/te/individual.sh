PROJ="control"
SAMPLE="PDL24"
# create a directory for stellarscope individual results
mkdir -p ../../../results/te/"$PROJ"/"$SAMPLE"/individual

# stellarscope analysis (pooling mode individual)
stellarscope assign \
  --exp_tag "$SAMPLE" \
  --outdir ../../../results/te/"$PROJ"/"$SAMPLE"/individual \
  --nproc 18 \
  --stranded_mode F \
  --whitelist ../../../results/te/"$PROJ"/"$SAMPLE"/star_alignment/Solo.out/Gene/filtered/barcodes.tsv \
  --pooling_mode individual \
  --reassign_mode best_exclude \
  --max_iter 500 \
  --updated_sam \
  ../../../results/te/"$PROJ"/"$SAMPLE"/cellsort/Aligned.sortedByCB.bam \
  resources/retro.hg38.v1.gtf
