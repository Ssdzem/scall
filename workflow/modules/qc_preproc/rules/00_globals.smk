# workflow/modules/qc_preproc/rules/00_globals.smk
import os
from glob import glob

ENV_DIR         = config.get("envs_dir", "../envs")
CONDA_ENV       = f"{ENV_DIR}/scanpy_pm.yaml"

COUNTS_DIR      = config.get("counts_dir", "results/hca_fetch/counts")
QC_ROOT         = config.get("qc_root", "results/qc_preproc")
LOGS_DIR        = config.get("logs_dir", "logs/qc_preproc")
LAYER           = config.get("layer", "filtered")
THREADS_DEFAULT = int(config.get("threads_default", 2))
H5AD_BASENAME   = config.get("h5ad_basename", "00_raw.h5ad")
LOAD_TAG        = config.get("load_tag", "01_loaded.h5ad")

SCRIPT_LOADER   = config.get("scripts", {}).get(
    "write_h5ad",
    "workflow/modules/qc_preproc/scripts/00_write_h5ad.py",
)

MIN_GENES       = int(config.get("min_genes", 200))
MIN_CELLS_GENE  = int(config.get("min_cells_gene", 10))
PROJECT         = config.get("project")  # None → all projects

def star_gene_layer_dir(w):
    # function is allowed in *input*; receives `wildcards`
    return f"{COUNTS_DIR}/{w.project}/{w.sample}/Gene/{LAYER}"

# Discover concrete final targets on disk (so master can target the module)
_final = []
if os.path.isdir(COUNTS_DIR):
    projects = [PROJECT] if PROJECT and PROJECT not in ("", "null") else \
               sorted([d for d in os.listdir(COUNTS_DIR)
                       if os.path.isdir(os.path.join(COUNTS_DIR, d))])
    for pr in projects:
        for p in glob(os.path.join(COUNTS_DIR, pr, "*", "Gene", LAYER)):
            sample = p.split(os.sep)[-3]
            _final.append(f"{QC_ROOT}/{pr}/{sample}/{H5AD_BASENAME}")
final_targets = sorted(set(_final))

# >>> aggregator in globals, like hca_fetch <<<
rule all:
    input:
        final_targets
