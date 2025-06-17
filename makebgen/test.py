import pandas as pd
import subprocess
import json
from pathlib import Path

# User inputs
input_file = Path("/Users/alish.palmos/PycharmProjects/mrcepid-makebgen/duplicates_only.txt") # Change this to your actual file path
columns_with_dx_ids = ["output_bcf", "output_bcf_idx", "output_vep", "output_vep_idx"]
archive_dir = "filtered_vcfs/chr21_old"  # Change this if you want a different archive folder

# Read the input file
df = pd.read_csv(input_file, sep="\t")

df.head()

# Identify duplicates based on chrom, start, end, vcf_prefix
dup_keys = ["chrom", "start", "end", "vcf_prefix"]
duplicates = df[df.duplicated(subset=dup_keys, keep=False)]

# Group duplicates by identifying keys
grouped = duplicates.groupby(dup_keys)

# Function to get creation timestamp from dx
def get_dx_creation_time(file_id):
    try:
        desc = subprocess.check_output(["dx", "describe", file_id, "--json"])
        return json.loads(desc)["created"]
    except Exception as e:
        print(f"Error getting creation time for {file_id}: {e}")
        return None

# Loop through each group of duplicates (should be 2 each)
for group_key, group_df in grouped:
    if len(group_df) != 2:
        print(f"Skipping non-pair duplicate group: {group_key}")
        continue

    row1, row2 = group_df.iloc[0], group_df.iloc[1]

    for col in columns_with_dx_ids:
        id1, id2 = row1[col], row2[col]

        if pd.isna(id1) or pd.isna(id2):
            continue

        time1 = get_dx_creation_time(id1)
        time2 = get_dx_creation_time(id2)

        if time1 is None or time2 is None:
            continue

        # Identify the older file
        if time1 < time2:
            old_file = id1
        else:
            old_file = id2

        # Move the old file
        try:
            subprocess.run(["dx", "mv", old_file, f"{archive_dir}/"], check=True)
            print(f"Moved {old_file} to {archive_dir}/")
        except subprocess.CalledProcessError as e:
            print(f"Failed to move {old_file}: {e}")
