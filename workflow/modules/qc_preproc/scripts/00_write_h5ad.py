#!/usr/bin/env python3
# SCALL qc_preproc: turn a 10x MTX dir (or .loom) into a standardized .h5ad

import argparse, sys
from pathlib import Path
import scanpy as sc

def read_any(path: Path):
    """Read either a 10x MTX directory or a .loom file into AnnData."""
    if path.is_dir():
        # 10x MTX directory. Keep gene symbols as var_names and ensure uniqueness.
        # Docs: scanpy.read_10x_mtx(path, var_names='gene_symbols', make_unique=True)
        ad = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)  # :contentReference[oaicite:1]{index=1}
        # Drop feature_types if present (not needed downstream)
        if "feature_types" in ad.var.columns:
            ad.var = ad.var.drop(columns=["feature_types"])
        return ad
    if path.suffix == ".loom":
        # Loom file. Default var_names in docs is "Gene"; some files store "gene_symbols".
        # Try symbols first, then fallback. Docs: scanpy.read_loom(...) :contentReference[oaicite:2]{index=2}
        try:
            ad = sc.read_loom(str(path), var_names="gene_symbols")  # may fail if that key doesn't exist
        except Exception:
            ad = sc.read_loom(str(path))  # default var_names='Gene'
        ad.var_names_make_unique()
        return ad
    raise ValueError(f"Unsupported input: {path} (expect 10x MTX folder or .loom)")

def main(argv=None):
    p = argparse.ArgumentParser(description="Load single-sample matrix → .h5ad")
    p.add_argument("--input", required=True, help="Path to 10x MTX directory or a .loom file")
    p.add_argument("--output", required=True, help="Output .h5ad path")
    p.add_argument("--min-genes", type=int, default=200, help="Filter cells with <min_genes")
    p.add_argument("--min-cells-gene", type=int, default=10, help="Filter genes detected in <min_cells")
    args = p.parse_args(argv)

    inp = Path(args.input)
    ad = read_any(inp)

    # basic QC filters (mirrors the notebook’s intent)
    sc.pp.filter_genes(ad, min_cells=args.min_cells_gene)  # keep genes seen in ≥N cells
    sc.pp.filter_cells(ad, min_genes=args.min_genes)       # keep cells with ≥N genes

    # write compressed h5ad (official API allows compression='gzip') :contentReference[oaicite:3]{index=3}
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.write_h5ad(str(out), compression="gzip")
    return 0

if __name__ == "__main__":
    sys.exit(main())
