# workflow/modules/qc_preproc/rules/01_load_h5ad.smk

rule load_h5ad:
    input:
        h5ad=lambda w: os.path.join(QC_ROOT, w.project, w.sample, H5AD_BASENAME)
    output:
        h5ad = f"{QC_ROOT}/{{project}}/{{sample}}/{LOAD_TAG}"
    log:
        f"{QC_ROOT}/{{project}}/{{sample}}/logs/01_nb_load_write_h5ad.log"
    conda:
        f"{config.get("envs_dir", "../envs")}/scanpy_pm.yaml"
    notebook:
        config.get("notebooks", {}).get(
    "load_h5ad",
    "workflow/modules/qc_preproc/vendor/hcc_notebooks/pipeline/01_load_h5ad.py.ipynb",
)
