#!/usr/bin/env python3
"""
Robust downloader for HCA/Azul-style curl pipelines.

- Wraps all `curl` invocations in safe defaults:
  * -L (follow redirects)           [curl manpage]
  * --fail (nonzero exit on HTTP>=400 so Snakemake retries) [curl manpage]
  * --retry-all-errors (retry on 429/5xx and more)          [everything.curl.dev]
  * --retry 8 + --retry-max-time 900 (overall backoff cap)  [curl manpage]
  * -C - (resume partial downloads)                         [curl manpage]
  * optional --limit-rate (be polite to servers)            [Stenberg blog]
- Adds small random jitter to avoid “thundering herd”.
"""
import argparse
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path

def build_shell_script(curl_cmd: str, rate_limit: str | None) -> str:
    """
    Creates a bash script that:
      1) enables aliases,
      2) defines a hardened curl alias,
      3) executes the provided curl pipeline verbatim.
    """
    alias = (
        "alias curl='command curl -L --fail "
        "--retry-all-errors --retry 8 --retry-max-time 900 -C -'"
    )
    if rate_limit:
        # e.g., 5M or 5000K
        alias = alias[:-1] + f" --limit-rate {rate_limit}'"

    script = f"""#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s expand_aliases
{alias}

# Execute the provided pipeline exactly as given
{curl_cmd}
"""
    tf = tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False)
    tf.write(script)
    tf.flush()
    tf.close()
    os.chmod(tf.name, 0o755)
    return tf.name

def maybe_sleep_jitter(jitter_range: str | None):
    """
    Sleep a small random amount (min-max seconds), e.g., "5-15".
    """
    rng = jitter_range or os.environ.get("DOWNLOAD_JITTER_RANGE", "3-12")
    try:
        lo_s, hi_s = [int(x) for x in rng.split("-")]
    except Exception:
        lo_s, hi_s = 3, 12
    time.sleep(random.uniform(lo_s, hi_s))

def download_project(curl_cmd: str, output_dir: Path, rate_limit: str | None, jitter_range: str | None):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Jitter to desynchronize many simultaneous jobs
    maybe_sleep_jitter(jitter_range)

    # Build & run hardened script in the output directory
    script_path = build_shell_script(curl_cmd, rate_limit)
    try:
        subprocess.run(["bash", script_path], cwd=str(output_dir), check=True)
    finally:
        # Best-effort cleanup
        try:
            os.unlink(script_path)
        except OSError:
            pass

    # Post-process: gzip any .fastq that arrived uncompressed
    for p in output_dir.iterdir():
        if not p.is_file():
            continue
        if p.name.endswith(".fastq.gz"):
            print(f"{p.name} is already compressed")
        elif p.name.endswith(".fastq"):
            print(f"Compressing {p.name} -> {p.name}.gz")
            subprocess.run(["gzip", str(p)], check=True)
        else:
            # Not a FASTQ: report extension (harmless for manifolds/aux files)
            ext = p.suffix[1:] if p.suffix else ""
            print(f"not a fastq, is {ext}")

def main():
    ap = argparse.ArgumentParser(description="Download a project via curl pipeline with safe defaults.")
    ap.add_argument("--curl-cmd", required=True, help="The exact curl pipeline to execute (quote it).")
    ap.add_argument("--output-dir", required=True, help="Directory to place downloaded files.")
    ap.add_argument("--rate-limit", default=os.environ.get("DOWNLOAD_LIMIT"),
                    help="Optional curl --limit-rate value (e.g., 6M).")
    ap.add_argument("--jitter-range", default=None,
                    help='Optional jitter range "min-max" seconds (default from env DOWNLOAD_JITTER_RANGE or 3-12).')
    args = ap.parse_args()

    download_project(
        curl_cmd=args.curl_cmd,
        output_dir=Path(args.output_dir),
        rate_limit=args.rate_limit,
        jitter_range=args.jitter_range,
    )

if __name__ == "__main__":
    main()
