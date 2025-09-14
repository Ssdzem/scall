PROJ="etoposide"
SAMPLE="50m_10d"
# run STAR alignment

STAR \
  --runThreadN 18 \
  --genomeDir /datos/sensence/emilio/scall/workflow/modules/te/resources/STAR_GRCh38.d1.vd1_gencode.v38 \
  --readFilesIn /datos/sensence/emilio/liver_sc/fastq_prepros/input/fastq/"$PROJ"/"$SAMPLE"/R2.fastq.gz /datos/sensence/emilio/liver_sc/fastq_prepros/input/fastq/"$PROJ"/"$SAMPLE"/R1.fastq.gz \
  --readFilesCommand gunzip -c \
  --soloCBwhitelist /datos/sensence/emilio/scall/workflow/modules/te/resources/whitelist_10x/3M-february-2018.txt \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
  --outSAMtype BAM SortedByCoordinate \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --limitOutSJcollapsed 5000000 \
  --outFilterMultimapNmax 500 \
  --outFilterMultimapScoreRange 5 \
  --outFileNamePrefix /datos/sensence/emilio/scall/results/te/"$PROJ"/"$SAMPLE"/star_alignment/
