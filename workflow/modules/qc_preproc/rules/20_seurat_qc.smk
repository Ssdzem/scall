# Include this from your qc_preproc module Snakefile when backend == "seurat"

from pathlib import Path

QC = config["qc_preproc"]
COUNTS_LAYER = QC.get("counts_layer", "filtered")

def starsolo_dir(wildcards):
    # STARsolo 10x-like output folder (matrix.mtx, features.tsv[.gz], barcodes.tsv[.gz])
    return f"results/hca_fetch/counts/{wildcards.sample}/Solo.out/Gene/{COUNTS_LAYER}"

rule qc_preproc__seurat_write:
    input:
        mtx = lambda w: f"{starsolo_dir(w)}/matrix.mtx",
        features = lambda w: f"{starsolo_dir(w)}/features.tsv.gz",
        barcodes = lambda w: f"{starsolo_dir(w)}/barcodes.tsv.gz",
    output:
        h5seurat = "results/qc_preproc/{sample}.h5seurat",
    params:
        mito_regex = QC.get("mito_regex", "^MT-"),
        min_genes  = QC.get("min_genes", 200),
        max_mt_pct = QC.get("max_mt_pct", 20),
    threads: 4
    conda: "envs/r-seurat.yaml"
    log:
        "logs/qc_preproc/{sample}.seurat_qc.log"
    script:
        "../scripts/qc_seurat.R"

rule qc_preproc__h5seurat_to_h5ad:
    input:
        h5seurat = rules.qc_preproc__seurat_write.output.h5seurat
    output:
        h5ad = "results/qc_preproc/{sample}.h5ad"
    threads: 1
    conda: "envs/r-seurat.yaml"
    log:
        "logs/qc_preproc/{sample}.seurat_convert.log"
    shell:
        r"""
        Rscript -e "suppressPackageStartupMessages(library(SeuratDisk));
                    SeuratDisk::Convert('{input.h5seurat}', dest='h5ad', overwrite=TRUE)"
        """

# Optional convenience target
rule qc_preproc__all:
    input:
        expand("results/qc_preproc/{{sample}}.{ext}", ext=["h5seurat"] + (["h5ad"] if QC.get("convert_to_h5ad", True) else []))
