#!/usr/bin/env python3
import argparse, os, gzip, pandas as pd, numpy as np
from scipy import io
from scipy.sparse import coo_matrix, vstack, csc_matrix

def mmread_any(p):
    if p.endswith(".gz"):
        with gzip.open(p, "rb") as fh: return io.mmread(fh)
    return io.mmread(p)

def mmwrite_any(p, M):
    if p.endswith(".gz"):
        with gzip.open(p, "wb") as fh: io.mmwrite(fh, M)
    else:
        io.mmwrite(p, M)

ap = argparse.ArgumentParser(description="Merge gene counts with TE-only (Stellarscope-style correction).")
ap.add_argument("--umi-labels", required=True)  # from label_umis_with_genes.py
ap.add_argument("--te-features", required=True)
ap.add_argument("--te-barcodes", required=True)
ap.add_argument("--gene-features", required=True)
ap.add_argument("--gene-barcodes", required=True)
ap.add_argument("--gene-mtx", required=True)      # features x cells (STARsolo matrix.mtx[.gz])
ap.add_argument("--outdir", required=True)
ap.add_argument("--out-gz", action="store_true")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# 1) Load UMI labels and keep TE-only
lab = pd.read_csv(args.umi_labels, sep="\t")
lab = lab[(lab["label"] == "TE-only") & lab["TE_id"].astype(str).ne("")][["CB","UB","TE_id"]]
# one count per unique molecule (CB,UB,TE_id)
lab = lab.drop_duplicates()

# 2) TE index (features and barcodes from Stellarscope)
# --- TE features (accept 1/2/3 columns; drop __no_feature)
te_feat_raw = pd.read_csv(args.te_features, sep="\t", header=None, dtype=str)
# drop special counters like __no_feature if present
te_feat_raw = te_feat_raw[~te_feat_raw.iloc[:, 0].astype(str).str.startswith("__")]

if te_feat_raw.shape[1] == 1:
    te_feat = pd.DataFrame({
        "te_id": te_feat_raw.iloc[:,0].astype(str),
        "te_name": te_feat_raw.iloc[:,0].astype(str),
        "type": "TE"
    })
elif te_feat_raw.shape[1] == 2:
    te_feat = te_feat_raw.copy()
    te_feat.columns = ["te_id","te_name"]
    te_feat["type"] = "TE"
elif te_feat_raw.shape[1] >= 3:
    te_feat = te_feat_raw.iloc[:, :3].copy()
    te_feat.columns = ["te_id","te_name","type"]
else:
    raise SystemExit("Unexpected TE features.tsv format (no columns).")

te_bcs = pd.read_csv(args.te_barcodes, header=None)[0].astype(str).tolist()

te_row = {t:i for i,t in enumerate(te_feat["te_id"].astype(str))}
te_col = {b:i for i,b in enumerate(te_bcs)}

# Count TE-only UMIs per (TE_id, CB)
grp = lab.groupby(["TE_id","CB"])["UB"].nunique().reset_index(name="n")
# Keep only cells present in TE barcodes
grp = grp[grp["CB"].isin(te_col)]
# Build sparse TE-only matrix (features x cells), fill zeros for missing rows/cols
rows = [te_row[t] for t in grp["TE_id"] if t in te_row]
cols = [te_col[c] for t,c in zip(grp["TE_id"], grp["CB"]) if t in te_row]
data = grp.loc[[t in te_row for t in grp["TE_id"]], "n"].astype(int).to_numpy()
TE_only = coo_matrix((data, (rows, cols)), shape=(len(te_feat), len(te_bcs))).tocsr()

# 3) Gene matrix + indices
gene_feat = pd.read_csv(args.gene_features, sep="\t", header=None)
if gene_feat.shape[1] == 3: gene_feat.columns = ["gene_id","gene_name","type"]
elif gene_feat.shape[1] == 2: gene_feat.columns = ["gene_id","gene_name"]; gene_feat["type"]="Gene"
else: raise SystemExit("Unexpected gene features.tsv format")
gene_bcs = pd.read_csv(args.gene_barcodes, header=None)[0].astype(str).tolist()
G = mmread_any(args.gene_mtx).tocsr()  # features x cells

# 4) Align barcodes (use intersection; keep order stable by gene barcodes)
bcs_inter = [b for b in gene_bcs if b in set(te_bcs)]
if not bcs_inter: raise SystemExit("No overlapping barcodes between gene and TE matrices.")
# column reindex maps
g_col = {b:i for i,b in enumerate(gene_bcs)}
t_col = te_col
gc_idx = [g_col[b] for b in bcs_inter]
tc_idx = [t_col[b] for b in bcs_inter]
G_i = G[:, gc_idx]
TE_i = TE_only[:, tc_idx]

# 5) Stack features: [genes; TE-only]
combined = vstack([G_i, TE_i], format="csr")

# 6) Write outputs
# Barcodes
out_bar = os.path.join(args.outdir, "barcodes.tsv" + (".gz" if args.out_gz else ""))
if args.out_gz:
    with gzip.open(out_bar, "wt") as f:
        pd.Series(bcs_inter).to_csv(f, index=False, header=False)
else:
    pd.Series(bcs_inter).to_csv(out_bar, index=False, header=False)

# Features
gene_feat2 = gene_feat.copy(); gene_feat2["type"] = gene_feat2.get("type", "Gene")
te_feat2   = te_feat.copy();   te_feat2["type"]   = te_feat2.get("type", "TE")
comb_feat  = pd.concat([gene_feat2, te_feat2], ignore_index=True)

out_feat = os.path.join(args.outdir, "features.tsv" + (".gz" if args.out_gz else ""))
if args.out_gz:
    with gzip.open(out_feat, "wt") as f:
        comb_feat.to_csv(f, sep="\t", index=False, header=False)
else:
    comb_feat.to_csv(out_feat, sep="\t", index=False, header=False)

# Matrix
out_mtx = os.path.join(args.outdir, "matrix.mtx" + (".gz" if args.out_gz else ""))
mmwrite_any(out_mtx, combined)

print("[ok] combined matrix written:")
print("  features:", out_feat)
print("  barcodes:", out_bar)
print("  matrix  :", out_mtx)
