# workflow/modules/qc_preproc/rules/10_write_h5ad.smk

rule write_h5ad:
    """
    Make {H5AD_BASENAME} from STARsolo 10x MTX dir.

    Uses Scanpy’s documented APIs:
      • scanpy.read_10x_mtx(path, var_names='gene_symbols', make_unique=True)
      • AnnData.write_h5ad(..., compression='gzip')
    """
    input:
        mtx_dir = star_gene_layer_dir
    output:
        h5ad = f"{QC_ROOT}/{{project}}/{{sample}}/{H5AD_BASENAME}"
    conda:
        CONDA_ENV
    threads:
        THREADS_DEFAULT
    log:
        f"{LOGS_DIR}/{{project}}/{{sample}}/00_load.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.h5ad})" "$(dirname {log})"
        python "{SCRIPT_LOADER}" \
          --input "{input.mtx_dir}" \
          --output "{output.h5ad}" \
          --min-genes {MIN_GENES} \
          --min-cells-gene {MIN_CELLS_GENE} \
          &> "{log}"
        """
