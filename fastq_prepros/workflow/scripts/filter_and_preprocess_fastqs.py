import sys
import os
import shutil
from glob import glob
import pandas as pd

def get_bundles_for_sample(metadata_csv, project, sample_id):
    df = pd.read_csv(metadata_csv)
    df.columns = df.columns.str.lower()
    # filter rows for this project + ident_sample
    subset = df[
        (df["proyect"] == project) &
        (df["ident_sample"] == sample_id)
    ]
    if subset.empty:
        raise ValueError(
            f"[!] No rows for project='{project}', ident_sample='{sample_id}'."
        )
    return subset["bundle_uuid"].tolist()

def filter_and_preprocess(raw_dir, sample_id, project, metadata_csv, output_r1, output_r2):
    # get all bundle_UUIDs for this ident_sample
    bundles = get_bundles_for_sample(metadata_csv, project, sample_id)

    os.makedirs(os.path.dirname(output_r1), exist_ok=True)

    # open outputs once, then append each lane's reads
    with open(output_r1, "wb") as r1_out, open(output_r2, "wb") as r2_out:
        for bundle in bundles:
            lane_dir = os.path.join(raw_dir, bundle)
            if not os.path.isdir(lane_dir):
                raise FileNotFoundError(f"[!] Missing directory: {lane_dir}")

            raw_fastqs = sorted(glob(os.path.join(lane_dir, "*.fastq.gz")))
            r1_files = [f for f in raw_fastqs if "R1" in f]
            r2_files = [f for f in raw_fastqs if "R2" in f]

            if not r1_files or not r2_files:
                raise ValueError(f"[!] No R1/R2 in {lane_dir}")

            for f in r1_files:
                with open(f, "rb") as fh: shutil.copyfileobj(fh, r1_out)
            for f in r2_files:
                with open(f, "rb") as fh: shutil.copyfileobj(fh, r2_out)

if __name__ == "__main__":
    import sys
    raw_dir     = sys.argv[1]
    sample_id   = sys.argv[2]     # now ident_sample
    project     = sys.argv[3]
    metadata    = sys.argv[4]
    output_r1   = sys.argv[5]
    output_r2   = sys.argv[6]

    # Step 1: Look up all bundle_uuids for this sample
    # Step 2: Concatenate across all lanes
    filter_and_preprocess(raw_dir, sample_id, project, metadata, output_r1, output_r2)
