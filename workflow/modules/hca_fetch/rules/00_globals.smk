# workflow/rules/00_globals.smk

import os
import pandas as pd
from snakemake.io import directory, temp

# ---- Configuration ----

raw_dir        = config["raw_dir"]
fastq_dir      = config["fastq_dir"]
index_dir      = config["index_dir"]
bam_dir        = config["bam_dir"]
counts_dir     = config["counts_dir"]
threads        = config["star_threads"]
sjdb_overhang  = config["sjdb_overhang"]
sa_nbases      = config["genomeSAindexNbases"]
wl_map         = config.get("whitelists", {})
qc_dir         = config["qc_dir"]
scripts        = config["scripts"]

# Phase switch (QC-only vs full pipeline)
_first_raw   = str(config.get("first_fastq", "false")).strip().lower()
first_fastq  = _first_raw in ("1", "true", "yes", "y")

# ---- Metadata ----
df_curls  = pd.read_csv(config["project_curls_path"])
df_sample = pd.read_csv(config["metadata_path"])

df_curls.columns  = [c.lower() for c in df_curls.columns]    # 'project','curls'
df_sample.columns = [c.lower() for c in df_sample.columns]   # 'bundle_uuid','proyect','ident_sample','chemistry'

project_list = df_curls["project"].unique().tolist()

sample_map = (
    df_sample
    .groupby("proyect")["ident_sample"]
    .apply(lambda x: sorted(x.unique().tolist()))
    .to_dict()
)

bundles_per_sample = (
    df_sample
    .groupby(["proyect", "ident_sample"])["bundle_uuid"]
    .apply(list)
    .to_dict()
)

curl_dict = df_curls.set_index("project")["curls"].to_dict()

# sample -> chemistry (assumes consistent chemistry per sample)
chem_map = (
    df_sample
    .groupby("ident_sample")["chemistry"]
    .agg(lambda x: x.unique().tolist()[0])
    .to_dict()
)

def wl_for_sample(sample: str) -> str:
    """Return whitelist path for a given sample based on its chemistry."""
    chem = chem_map[sample]
    try:
        return wl_map[chem]
    except KeyError:
        raise ValueError(f"No whitelist configured for chemistry '{chem}'")

# 10x CB/UMI defaults (keep aligned with your earlier choices)
param_map = {
    "10x3'v2": {"cb": 16, "umi": 10},
    "10x3'v3": {"cb": 16, "umi": 12},
    "10x5'v1": {"cb": 16, "umi": 10},
}

project_cfg = [config["project"]] if config.get("project") else project_list

# Build targets (full vs QC-only)
all_targets, vibe_targets = [], []
for proj in project_cfg:
    for sample in sample_map[proj]:
        all_targets += [
            os.path.join(counts_dir, proj, f"{sample}_matrix.txt"),
            os.path.join(bam_dir,    proj, f"{sample}.bam"),
        ]
        vibe_targets += [
            os.path.join(qc_dir, proj, f"{sample}_vibe.json"),
            os.path.join(qc_dir, proj, f"{sample}_chemistry.txt"),
        ]

final_targets = vibe_targets if first_fastq else all_targets

# Top-level aggregation target
rule all:
    input:
        final_targets
