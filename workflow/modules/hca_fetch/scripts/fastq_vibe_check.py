#!/usr/bin/env python3
"""
fastq_vibe_check.py

Purpose
-------
A fast “vibe check” for 10x-style scRNA-seq FASTQs (R1/R2 and optional I2).
It samples reads to answer:
  • What are the observed read-length distributions in R1 and R2?
  • Which 10x barcode inclusion list (“whitelist”) fits best (v2 vs v3/3.1 vs v4)?
  • Is there a 3′ vs 5′ hint (polyT presence near start of R2 for 3′ assays)?
  • Is there an I2 file present (dual index, typical of v3.1/v4)?
  • If R1 contains mixed 26/28 lengths (merged runs), warn about splitting.

Outputs
-------
  1) JSON report (machine-readable): recommended STARsolo parameters, counts, warnings
  2) chemistry.txt: single-line “chemistry” string you can consume in downstream rules

Background
----------
• 10x 3′ GEX: R1 carries CB+UMI (16+10 in v2; 16+12 in v3/v4).
• STARsolo expects --readFilesIn R2 R1 (cDNA first, then barcode read).
• 3′ assays often show polyT early in R2 (heuristic signal; not a hard rule).
• Mixed R1 lengths 26/28 often come from mixing v2 and v3 runs; better to separate.
"""

import argparse
import gzip
import json
import os
import re
import sys
from collections import Counter
from typing import Dict, Iterable, Optional, Set, Tuple

# ---------- Utility helpers ----------

def open_auto(path: Optional[str]):
    """
    Open a text file transparently whether it’s gzipped or not.
    Returns a file handle or None if the path is falsy or missing.

    We avoid external deps to keep the tool light and easy to run on HPC nodes.
    """
    if not path or not os.path.exists(path):
        return None
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fastq_sequences(path: str, max_reads: Optional[int] = None) -> Iterable[str]:
    """
    Stream the SEQUENCE lines from a FASTQ (line 2, 6, 10, ...).
    Stops after `max_reads` sequences for speed.

    FASTQ layout reminder:
      1: @header
      2: sequence
      3: + (optional header)
      4: qualities
    """
    fh = open_auto(path)
    if fh is None:
        return
    try:
        i = 0
        taken = 0
        for line in fh:
            # Zero-based index: the sequence is every 4th line with (i % 4) == 1
            if (i % 4) == 1:
                yield line.rstrip()
                taken += 1
                if max_reads is not None and taken >= max_reads:
                    break
            i += 1
    finally:
        fh.close()


def load_whitelist_as_set(path: Optional[str]) -> Set[str]:
    """
    Load a barcode inclusion list (one barcode per line) into a Python set.
    Skips empty lines and comments. Missing/empty path -> empty set.

    10x calls these “barcode inclusion lists” (formerly "whitelists").
    """
    if not path:
        return set()
    fh = open_auto(path)
    if fh is None:
        return set()
    try:
        s = set()
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            s.add(t)
        return s
    finally:
        fh.close()


# ---------- Sampling/statistics helpers ----------

def fraction_polyT_first_k(seq_iter: Iterable[str], k: int = 20, min_run: int = 10, limit: int = 200_000) -> Tuple[float, int]:
    """
    Heuristic: among the first `k` bases of R2 (cDNA read), how many reads have
    a run of >= `min_run` T’s? 3′ GEX libraries often show polyT early in R2.

    Returns:
      (fraction_with_polyT, n_checked)
    """
    pat = re.compile(r"T{%d,}" % min_run)
    n = 0
    hits = 0
    for seq in seq_iter:
        n += 1
        if pat.search(seq[:k]):
            hits += 1
        if n >= limit:
            break
    return (hits / n) if n else 0.0, n


def length_histogram(seq_iter: Iterable[str], limit: int = 200_000) -> Tuple[Dict[int, int], int]:
    """
    Build a histogram of sequence lengths from a stream of sequences.
    Returns:
      ({length: count, ...}, n_sampled)
    """
    c = Counter()
    n = 0
    for s in seq_iter:
        c[len(s)] += 1
        n += 1
        if n >= limit:
            break
    return dict(c), n


def first16_unique_set(seq_iter: Iterable[str], limit: int = 1_000_000) -> Tuple[Set[str], int]:
    """
    Collect unique 16-mers from the start of R1 sequences (barcode span).
    Returns the set and how many reads were considered.
    """
    s = set()
    n = 0
    for seq in seq_iter:
        if len(seq) >= 16:
            s.add(seq[:16])
        n += 1
        if n >= limit:
            break
    return s, n


def vote_from_whitelists(v2_hits: int, v3_hits: int, v4_hits: int) -> Tuple[str, float]:
    """
    Pick the chemistry (3' v2, 3' v3, 3' v4) with the largest number of 16-mer overlaps.
    Returns:
      (chemistry_string or "unknown", confidence_fraction)
    """
    candidates = [("3p_v2", v2_hits), ("3p_v3", v3_hits), ("3p_v4", v4_hits)]
    name, val = max(candidates, key=lambda x: x[1])
    total = max(1, v2_hits + v3_hits + v4_hits)
    conf = val / total
    return (name if val > 0 else "unknown", conf)


# ---------- Main ----------

def main():
    p = argparse.ArgumentParser(description="FASTQ vibe check for 10x scRNA-seq")
    p.add_argument("--r1", required=True, help="Path to R1 (barcode+UMI) FASTQ[.gz]")
    p.add_argument("--r2", required=True, help="Path to R2 (cDNA) FASTQ[.gz]")
    p.add_argument("--i2", default=None, help="Optional path to I2 FASTQ[.gz] (dual index)")
    p.add_argument("--wl-v2", required=True, help="Path to 3' v2 inclusion list (737K-august-2016.txt)")
    p.add_argument("--wl-v3", required=True, help="Path to 3' v3/v3.1 inclusion list (3M-february-2018*.txt.gz)")
    p.add_argument("--wl-v4", default=None, help="Optional path to 3' v4 inclusion list (3M-3pgex-may-2023*.txt.gz)")
    p.add_argument("--metadata-chem", default="unknown", help="Chemistry label from your metadata (for cross-check)")
    p.add_argument("--max-reads", type=int, default=1_000_000, help="Max reads to sample per measurement")
    p.add_argument("--json-out", required=True, help="Output JSON path")
    p.add_argument("--chem-out", required=True, help="Output single-line chemistry file path")
    args = p.parse_args()

    # Sampling budgets per signal (tuned for fast turn-around on big files).
    LEN_N   = min(args.max_reads, 200_000)   # reads sampled for length histograms
    POLYT_N = min(args.max_reads, 200_000)   # reads checked for polyT heuristic
    CBSET_N = min(args.max_reads, 1_000_000) # reads for 16-mer overlap with lists

    # 1) R1/R2 length histograms — catches obviously mixed inputs early
    r1_hist, r1n = length_histogram(iter_fastq_sequences(args.r1, max_reads=LEN_N), limit=LEN_N)
    r2_hist, r2n = length_histogram(iter_fastq_sequences(args.r2, max_reads=LEN_N), limit=LEN_N)

    # 2) Unique first-16 collection from R1 — we’ll intersect with the inclusion lists
    first16, seen_n = first16_unique_set(iter_fastq_sequences(args.r1, max_reads=CBSET_N), limit=CBSET_N)

    # 3) Load inclusion lists (a.k.a. “whitelists”)
    wl_v2 = load_whitelist_as_set(args.wl_v2)  # 737K-august-2016.txt (3' v2)
    wl_v3 = load_whitelist_as_set(args.wl_v3)  # 3M-february-2018*.txt.gz (3' v3/3.1)
    wl_v4 = load_whitelist_as_set(args.wl_v4) if args.wl_v4 else set()  # 3' v4 (optional)

    # Count overlaps: how many observed 16-mers exist in each inclusion list?
    v2_hits = len(first16 & wl_v2)
    v3_hits = len(first16 & wl_v3)
    v4_hits = len(first16 & wl_v4) if wl_v4 else 0

    inferred_short, conf = vote_from_whitelists(v2_hits, v3_hits, v4_hits)

    # 4) 3′ vs 5′ soft signal: polyT near the start of R2 (many 3′ libraries show this)
    polyT_frac, polyT_n = fraction_polyT_first_k(
        iter_fastq_sequences(args.r2, max_reads=POLYT_N),
        k=20,       # look at the first 20 bases of R2
        min_run=10, # consider >=10 consecutive T’s as a “hit”
        limit=POLYT_N
    )

    # 5) Dual-index presence: if I2 exists on disk, we record it (v3.1/v4 often dual-index)
    has_i2 = bool(args.i2 and os.path.exists(args.i2))

    # 6) Quick red flag for mixed R1 lengths (26 and 28 both present in non-trivial numbers)
    mixed_r1_flag = (26 in r1_hist and 28 in r1_hist) and (r1_hist[26] > 0 and r1_hist[28] > 0)

    # 7) Translate the “short” label into a human-friendly chemistry string
    #    and assemble a STARsolo suggestion block (not executing STAR here).
    recommendation = {
        # STARsolo expects cDNA read first and barcode read second:
        #   --readFilesIn <R2> <R1>
        "readFilesIn_order": "R2_then_R1",
        "CBlen": 16,  # 10x 3′ CB length
    }

    if inferred_short == "3p_v2":
        # 3′ v2: 16bp CB + 10bp UMI, R1=26 used for CB+UMI
        recommendation.update({
            "chemistry": "10x_3prime_v2",
            "UMIlen": 10,
            "soloBarcodeReadLength": 26,  # let STARsolo know R1 effective length
            "whitelist": os.path.abspath(args.wl_v2),
        })
    elif inferred_short == "3p_v3":
        # 3′ v3/3.1: 16bp CB + 12bp UMI, R1=28 used for CB+UMI
        recommendation.update({
            "chemistry": "10x_3prime_v3_or_v3.1",
            "UMIlen": 12,
            "soloBarcodeReadLength": 28,
            "whitelist": os.path.abspath(args.wl_v3),
        })
    elif inferred_short == "3p_v4":
        # 3′ v4 (GEM-X): still 16+12; inclusion list differs from v3/3.1
        recommendation.update({
            "chemistry": "10x_3prime_v4",
            "UMIlen": 12,
            "soloBarcodeReadLength": 28,
            "whitelist": os.path.abspath(args.wl_v4) if args.wl_v4 else None,
        })
    else:
        recommendation.update({
            "chemistry": "unknown",
            "UMIlen": None,
            "soloBarcodeReadLength": None,
            "whitelist": None,
        })

    # STAR options that emulate Cell Ranger >= v4 behavior in STARsolo
    # (we only *suggest* them here; you’ll apply them in your mapping rule)
    recommendation["star_like_cellranger4"] = {
        "clipAdapterType": "CellRanger4",
        "outFilterScoreMin": 30,
        "soloCBmatchWLtype": "1MM_multi_Nbase_pseudocounts",
        "soloUMIfiltering": "MultiGeneUMI_CR",
        "soloUMIdedup": "1MM_CR"
    }

    # Cross-checks and user-facing warnings
    warnings = []
    if mixed_r1_flag:
        warnings.append(
            "Mixed R1 effective lengths (26 and 28) detected."
        )
    if inferred_short == "unknown":
        warnings.append(
            "Whitelist vote is inconclusive; consider a small dual-run STARsolo smoke test (v2 vs v3) and pick the better Summary."
        )
    if has_i2 and recommendation["chemistry"] not in ("10x_3prime_v3_or_v3.1", "10x_3prime_v4", "unknown"):
        warnings.append(
            "I2 detected (dual-index) but inferred chemistry is not v3.1/v4. Re-check metadata."
        )
    if args.metadata_chem and args.metadata_chem.lower() != "unknown" and recommendation["chemistry"] != "unknown":
        if args.metadata_chem.lower() not in recommendation["chemistry"].lower():
            warnings.append(
                f"Metadata chemistry '{args.metadata_chem}' disagrees with inferred '{recommendation['chemistry']}'."
            )

    # Assemble the JSON report
    report = {
        "inputs": {"R1": args.r1, "R2": args.r2, "I2": args.i2},
        "sample_counts": {
            "r1_len_hist": r1_hist, "r1_sampled": r1n,
            "r2_len_hist": r2_hist, "r2_sampled": r2n,
            "first16_seen": len(first16), "first16_sampled": seen_n
        },
        "whitelist_hits": {"v2": v2_hits, "v3": v3_hits, "v4": v4_hits},
        "polyT_first20_frac": polyT_frac,
        "polyT_reads_checked": polyT_n,
        "inferred": {"chemistry": recommendation["chemistry"], "confidence": conf},
        "has_dual_index_I2": has_i2,
        "recommendation": recommendation,
        "warnings": warnings,
        "notes": [
            "STARsolo expects --readFilesIn <R2> <R1> (cDNA first, barcode/UMI second).",
            "Set --soloCBlen=16 and --soloUMIlen per chemistry; optionally set --soloBarcodeReadLength to 26 (v2) or 28 (v3+/v4)."
        ]
    }

    # Write outputs
    os.makedirs(os.path.dirname(args.json_out), exist_ok=True)
    with open(args.json_out, "w") as f:
        json.dump(report, f, indent=2)

    os.makedirs(os.path.dirname(args.chem_out), exist_ok=True)
    with open(args.chem_out, "w") as f:
        f.write(recommendation["chemistry"] + "\n")


if __name__ == "__main__":
    main()
