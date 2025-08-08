#!/usr/bin/env python3
import sys
import os
import subprocess

def download_project(curl_cmd, output_dir):
    # 1. Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # 2. Run your curl pipeline inside it
    subprocess.run(curl_cmd, shell=True, cwd=output_dir, check=True)

    # 3. Inspect each file in the download folder
    for fname in os.listdir(output_dir):
        fpath = os.path.join(output_dir, fname)
        if not os.path.isfile(fpath):
            continue

        if fname.endswith(".fastq.gz"):
            # already compressed; leave it
            print(f"{fname} is already compressed")

        elif fname.endswith(".fastq"):
            # compress uncompressed FASTQ
            print(f"Compressing {fname} → {fname}.gz")
            subprocess.run(["gzip", fpath], check=True)

        else:
            # not a FASTQ: report its extension
            ext = fname.rsplit(".", 1)[-1] if "." in fname else ""
            print(f"not a fastq, is {ext}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: download_project.py <curl_cmd> <output_dir>", file=sys.stderr)
        sys.exit(1)

    curl_cmd, output_dir = sys.argv[1], sys.argv[2]
    download_project(curl_cmd, output_dir)
