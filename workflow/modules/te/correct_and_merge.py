#!/usr/bin/env python3
import argparse, os, re, gzip
from collections import defaultdict
import pandas as pd
import numpy as np
import pysam
from scipy import io
from scipy.sparse import coo_matrix, vstack

def mmread_any(p):
    if p.endswith(".gz"):
        with gzip.open(p, "rb") as fh: return io.mmread(fh).tocsr()
    return io.mmread(p).tocsr()

def mmwrite_any(p, M):
    if p.endswith(".gz"):
        with gzip.open(p, "wb") as fh: io.mmwrite(fh, M)
    else:
        io.mmwrite(p, M)

def load_te_features(path):
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    mask = ~df.iloc[:,0].astype(str).str.startswith("__")
    df = df[mask]
    if df.shape[1] == 1:
        df = pd.DataFrame({"te_id": df.iloc[:,0], "name": df.iloc[:,0], "type": "TE"})
    elif df.shape[1] >= 2:
        df = df.iloc[:,:2]; df.columns = ["te_id","name"]; df["type"] = "TE"
    return df.reset_index(drop=True), mask

def load_gene_features(path):
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    # star/10x often 3 cols: gene_id, gene_name, type
    if df.shape[1] >= 2:
        gene_id = df.iloc[:,0]
        gene_name = df.iloc[:,1].fillna("")
        name = np.where(gene_name.str.len()>0, gene_name, gene_id)
        return pd.DataFrame({"gene_id": gene_id, "name": name, "type": "Gene"})
    else:
        # fallback: single col
        return pd.DataFrame({"gene_id": df.iloc[:,0], "name": df.iloc[:,0], "type": "Gene"})

def parse_attrs(s):
    return {k:v for k,v in re.findall(r'(\S+)\s+"([^"]*)"', s)}

def build_te_index(gtf):
    idx = {}
    with open(gtf) as fh:
        for ln in fh:
            if not ln or ln[0] == "#": continue
            chrom,_,_,beg,end,_,_,_,attrs = ln.rstrip("\n").split("\t")
            a = parse_attrs(attrs)
            te_id = a.get("transcript_id") or a.get("gene_id") or a.get("gene_name") or a.get("ID") or a.get("Name")
            if not te_id: continue
            idx.setdefault(chrom, []).append((int(beg)-1, int(end), te_id))
    return idx

def te_hits(idx, chrom, s, e):
    for a,b,t in idx.get(chrom, ()):
        if a < e and s < b:  # overlap
            yield t

def main():
    ap = argparse.ArgumentParser(description="Stellarscope-style TE correction and merge (no intermediate UMI labels).")
    ap.add_argument("--updated-bam", required=True)
    ap.add_argument("--te-gtf", required=True)
    ap.add_argument("--te-features", required=True)
    ap.add_argument("--te-barcodes", required=True)
    ap.add_argument("--te-mtx", required=True)      # Stellarscope *-TE_counts.mtx
    ap.add_argument("--gene-features", required=True)
    ap.add_argument("--gene-barcodes", required=True)
    ap.add_argument("--gene-mtx", required=True)    # STARsolo matrix.mtx(.gz or not)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-gz", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # load indices
    te_feat, mask = load_te_features(args.te_features)
    te_bcs  = pd.read_csv(args.te_barcodes, header=None, dtype=str)[0].tolist()
    gene_feat = load_gene_features(args.gene_features)
    gene_bcs  = pd.read_csv(args.gene_barcodes, header=None, dtype=str)[0].tolist()

    # matrixes
    TE_raw = mmread_any(args.te_mtx)     # features x cells (TE)
    G      = mmread_any(args.gene_mtx)   # features x cells (genes)

    # remove the __no_features
    TE_raw = TE_raw[mask.values, :]

    # align barcodes: keep intersection in gene order (Seurat-friendly)
    te_set = set(te_bcs)
    bcs = [b for b in gene_bcs if b in te_set]
    if not bcs:
        raise SystemExit("No overlapping barcodes between TE and Gene matrices.")
    g_col = {b:i for i,b in enumerate(gene_bcs)}
    t_col = {b:i for i,b in enumerate(te_bcs)}
    gc_idx = [g_col[b] for b in bcs]
    tc_idx = [t_col[b] for b in bcs]
    G = G[:, gc_idx]
    TE_raw = TE_raw[:, tc_idx]

    # build TE interval index
    idx = build_te_index(args.te_gtf)

    # scan BAM once, collect for each UMI whether it hits any TE and whether it has GX
    bam = pysam.AlignmentFile(args.updated_bam, "rb")
    umi_has_gene = set()                   # (CB,UB)
    umi_te = defaultdict(set)              # (CB,UB) -> set(TE_id)

    def gx_present(rec):
        # treat "-" or "." as missing; check both GX and gx
        for tag in ("GX","gx"):
            if rec.has_tag(tag):
                val = str(rec.get_tag(tag)).strip()
                if val and val not in ("-","."): return True
        return False

    for a in bam:  # sequential; no index required
        if a.is_unmapped or a.is_secondary or a.is_supplementary: continue
        try:
            cb = a.get_tag("CB"); ub = a.get_tag("UB")
        except KeyError:
            continue
        if gx_present(a):
            umi_has_gene.add((cb,ub))
        chrom = bam.get_reference_name(a.reference_id)
        for s,e in a.get_blocks():
            for te in te_hits(idx, chrom, s, e):
                umi_te[(cb,ub)].add(te)

    # build correction counts (TE x cell): number of unique UMIs that are TE & Gene
    te_row = {t:i for i,t in enumerate(te_feat["te_id"].astype(str))}
    te_col = {b:i for i,b in enumerate(bcs)}
    rows, cols, data = [], [], []
    corr = defaultdict(int)
    for (cb,ub), tes in umi_te.items():
        if (cb,ub) in umi_has_gene and cb in te_col:
            for te in tes:
                if te in te_row:
                    corr[(te_row[te], te_col[cb])] += 1
    if corr:
        for (r,c), v in corr.items():
            rows.append(r); cols.append(c); data.append(v)
        TE_corr = coo_matrix((data,(rows,cols)), shape=TE_raw.shape, dtype=np.int32).tocsr()
    else:
        TE_corr = TE_raw.copy(); TE_corr.data[:] = 0

    # corrected TE matrix
    TE_corrd = TE_raw - TE_corr
    TE_corrd.data = np.maximum(TE_corrd.data, 0)

    # make 2-col features: name, type
    gene_feats_2 = gene_feat[["name","type"]].copy()
    te_feats_2   = te_feat[["name","type"]].copy()
    comb_feat = pd.concat([gene_feats_2, te_feats_2], ignore_index=True)
    comb_feat = pd.concat([comb_feat.iloc[:,0], comb_feat], axis=1)
    comb_feat.columns = ["id", "name", "type"]

    # combine matrices: [genes; TE_corrected]
    combined = vstack([G, TE_corrd], format="csr")

    # write outputs
    out_bar = os.path.join(args.outdir, "barcodes.tsv.gz")
    pd.Series(bcs).to_csv(out_bar, index=False, header=False, compression="gzip")

    out_feat = os.path.join(args.outdir, "features.tsv.gz")
    comb_feat.to_csv(out_feat, sep="\t", index=False, header=False, compression="gzip")

    out_mtx = os.path.join(args.outdir, "matrix.mtx" + (".gz" if args.out_gz else ""))
    mmwrite_any(out_mtx, combined)

    print("[ok] wrote:")
    print("  barcodes:", out_bar)
    print("  features:", out_feat, "(2 columns: name, type)")
    print("  matrix  :", out_mtx)

if __name__ == "__main__":
    main()
