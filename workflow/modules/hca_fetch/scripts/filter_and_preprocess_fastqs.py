#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
filter_and_preprocess_fastqs.py

Concatenate R1/R2 across all PRUNED bundle_uuids for a given sample,
and also preserve the raw PRUNED bundle_uuids per sample in:
  <keep_dir>/<bundle_uuid>/

Supports common naming schemes:
- Illumina/mkfastq: ..._L001_R1_001.fastq.gz (and R2; I1/I2 are index reads)
- SRA/fasterq-dump: SRRxxxxxxx_1.fastq(.gz) / SRRxxxxxxx_2.fastq(.gz)
- Aliases: read1/read2, pair1/pair2, fwd/rev, forward/reverse

Inputs may be .fastq(.gz) or .fq(.gz). Uncompressed inputs are gzipped
on-the-fly and appended as separate gzip members; pre-gzipped inputs
are byte-copied (valid multi-member .gz).
"""

import argparse
import gzip
import os
import re
import shutil
import sys
from glob import glob
from typing import List, Dict

import pandas as pd


# ---------------------- helpers: metadata ---------------------- #
def get_bundles_for_sample(metadata_csv: str, project: str, sample_id: str) -> List[str]:
    """Return the list of bundle_uuids for a (project, ident_sample)."""
    df = pd.read_csv(metadata_csv)
    df.columns = df.columns.str.lower()
    subset = df[(df["proyect"] == project) & (df["ident_sample"] == sample_id)]
    if subset.empty:
        raise ValueError(
            f"[!] No rows for project='{project}', ident_sample='{sample_id}'. "
            f"Check metadata file: {metadata_csv}"
        )
    # Keep ordering as in metadata
    return subset["bundle_uuid"].astype(str).tolist()


# ---------------------- helpers: file discovery ---------------------- #
def discover_fastqs(lane_dir: str) -> List[str]:
    """Return all plausible FASTQ files in a lane directory."""
    exts = ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq")
    files: List[str] = []
    for pat in exts:
        files.extend(glob(os.path.join(lane_dir, pat)))
    return sorted(files)


# ---------------------- helpers: R1/R2/I1/I2 classification ---------------------- #
# Boundaries: start or one of _-. ; trailing: followed by _-. or end.
BOUND = r"(?:^|[_\.-])"
TRAIL = r"(?=[_\.-]|$)"

R1_PATTERNS = [
    rf"{BOUND}R1{TRAIL}",
    rf"{BOUND}read1{TRAIL}",
    rf"{BOUND}pair1{TRAIL}",
    rf"{BOUND}fwd{TRAIL}",
    rf"{BOUND}forward{TRAIL}",
    rf"{BOUND}1{TRAIL}",  # SRA style: *_1.fastq(.gz)
]
R2_PATTERNS = [
    rf"{BOUND}R2{TRAIL}",
    rf"{BOUND}read2{TRAIL}",
    rf"{BOUND}pair2{TRAIL}",
    rf"{BOUND}rev{TRAIL}",
    rf"{BOUND}reverse{TRAIL}",
    rf"{BOUND}2{TRAIL}",  # SRA style: *_2.fastq(.gz)
]
I1_PATTERNS = [rf"{BOUND}I1{TRAIL}", rf"{BOUND}index1{TRAIL}"]
I2_PATTERNS = [rf"{BOUND}I2{TRAIL}", rf"{BOUND}index2{TRAIL}"]

R1_RE = [re.compile(p, re.IGNORECASE) for p in R1_PATTERNS]
R2_RE = [re.compile(p, re.IGNORECASE) for p in R2_PATTERNS]
I1_RE = [re.compile(p, re.IGNORECASE) for p in I1_PATTERNS]
I2_RE = [re.compile(p, re.IGNORECASE) for p in I2_PATTERNS]


def _match_any(regex_list: List[re.Pattern], s: str) -> bool:
    return any(r.search(s) for r in regex_list)


def classify_fastq_role(fname: str) -> str:
    """
    Return one of: 'R1', 'R2', 'I1', 'I2', or 'UNKNOWN'
    based on filename tokens (case-insensitive).
    """
    base = os.path.basename(fname)
    if _match_any(I1_RE, base):
        return "I1"
    if _match_any(I2_RE, base):
        return "I2"
    if _match_any(R1_RE, base):
        return "R1"
    if _match_any(R2_RE, base):
        return "R2"
    return "UNKNOWN"


# ---------------------- helpers: IO ---------------------- #
def safe_link_or_copy(src: str, dst: str) -> None:
    """
    Try to hard-link to save space; if that fails (cross-device, perms),
    fall back to a regular copy that preserves metadata.
    """
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        if os.path.exists(dst):
            os.remove(dst)
        os.link(src, dst)
    except Exception:
        shutil.copy2(src, dst)


def _is_gz(path: str) -> bool:
    return path.endswith(".gz") or path.endswith(".gzip")


def append_into_gz_aware_output(src: str, out_handle) -> None:
    """
    Append src FASTQ into output stream `out_handle`, which is an already-open
    binary file for the final .fastq.gz.

    - If src is gzipped: byte-copy (becomes another gzip member).
    - If src is plain text: compress on-the-fly into a new gzip member.
    """
    if _is_gz(src):
        with open(src, "rb") as fh:
            shutil.copyfileobj(fh, out_handle)
    else:
        # Write a new gzip member at current end of file
        with open(src, "rb") as fin, gzip.GzipFile(fileobj=out_handle, mode="wb") as gzout:
            shutil.copyfileobj(fin, gzout)


# ---------------------- main routine ---------------------- #
def filter_and_preprocess(
    raw_dir: str,
    sample_id: str,
    project: str,
    metadata_csv: str,
    output_r1: str,
    output_r2: str,
    keep_dir: str,
) -> None:
    """
    - Look up PRUNED bundle_uuids for this (project, sample)
    - Concatenate R1 across bundles into output_r1 (.fastq.gz)
    - Concatenate R2 across bundles into output_r2 (.fastq.gz)
    - Preserve PRUNED raw R1/R2 files under keep_dir/<bundle_uuid>/
    """
    bundles = get_bundles_for_sample(metadata_csv, project, sample_id)
    os.makedirs(os.path.dirname(output_r1), exist_ok=True)
    os.makedirs(keep_dir, exist_ok=True)

    # Open outputs once; we always produce gzipped outputs.
    with open(output_r1, "wb") as r1_out, open(output_r2, "wb") as r2_out:
        for bundle in bundles:
            lane_dir = os.path.join(raw_dir, bundle)
            if not os.path.isdir(lane_dir):
                raise FileNotFoundError(f"[!] Missing directory: {lane_dir}")

            raw_fastqs = discover_fastqs(lane_dir)
            if not raw_fastqs:
                raise ValueError(f"[!] No FASTQs found in {lane_dir}")

            # classify files
            roles: Dict[str, str] = {f: classify_fastq_role(f) for f in raw_fastqs}
            r1_files = [f for f, role in roles.items() if role == "R1"]
            r2_files = [f for f, role in roles.items() if role == "R2"]
            unknown  = [f for f, role in roles.items() if role == "UNKNOWN"]

            if unknown:
                # Don't fail immediately; we just ignore UNKNOWN (e.g., md5s)
                sys.stderr.write(
                    f"[warn] Ignoring non-R1/R2 files in {lane_dir}: "
                    f"{', '.join(os.path.basename(x) for x in unknown)}\n"
                )

            if not r1_files or not r2_files:
                found = ", ".join(os.path.basename(x) for x in raw_fastqs) or "NONE"
                raise ValueError(
                    f"[!] No R1/R2 in {lane_dir} (after classification). Found: {found}"
                )

            # 1) Append to concatenated outputs
            for f in r1_files:
                append_into_gz_aware_output(f, r1_out)
            for f in r2_files:
                append_into_gz_aware_output(f, r2_out)

            # 2) Preserve PRUNED raw R1/R2s under keep_dir/<bundle_uuid>/
            dest_bundle_dir = os.path.join(keep_dir, bundle)
            for f in (r1_files + r2_files):
                dest = os.path.join(dest_bundle_dir, os.path.basename(f))
                safe_link_or_copy(f, dest)

    print(
        f"[OK] Sample '{sample_id}' ({project}): "
        f"concatenated {len(bundles)} bundles -> "
        f"{output_r1} & {output_r2}; "
        f"preserved raw PRUNED bundles under {keep_dir}"
    )


# ---------------------- CLI ---------------------- #
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="filter_and_preprocess_fastqs.py",
        description=("Concatenate R1/R2 across PRUNED bundles for a sample and "
                     "preserve per-bundle raw R1/R2 files.")
    )
    p.add_argument("--raw-dir", required=True,
                   help="Path to raw/<project> directory containing bundle UUID subdirs.")
    p.add_argument("--sample", "--ident-sample", dest="sample_id", required=True,
                   help="ident_sample for this run.")
    p.add_argument("--project", required=True, help="Project name.")
    p.add_argument("--metadata-csv", required=True,
                   help="CSV with columns: proyect, ident_sample, bundle_uuid, chemistry, ...")
    p.add_argument("--out-r1", required=True, help="Output concatenated R1.fastq.gz.")
    p.add_argument("--out-r2", required=True, help="Output concatenated R2.fastq.gz.")
    p.add_argument("--keep-dir", required=True,
                   help="Directory to store preserved per-bundle R1/R2 files.")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    filter_and_preprocess(
        raw_dir=args.raw_dir,
        sample_id=args.sample_id,
        project=args.project,
        metadata_csv=args.metadata_csv,
        output_r1=args.out_r1,
        output_r2=args.out_r2,
        keep_dir=args.keep_dir,
    )
