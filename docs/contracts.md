# Data Contracts (initial draft)
## QC output (example)
file: results/qc_preproc/adata_merged.h5ad
obs required: ["sample","batch","cell_type","cluster"]
## TE output (example)
file: results/te/te_counts.parquet
index: cell barcode == adata.obs_names
