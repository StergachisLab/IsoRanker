import os
import pandas as pd
import gzip


# def update_files_with_haplotype_info(sample_info_with_haplotype_location, read_stats_path, output_dir):
#     """
#     Updates read statistics and sample metadata with haplotype information.

#     This function performs the following steps:
#     - Loads `sample_info_with_haplotype_location`, which contains paths to haplotype assignment files.
#     - Loads a read-level summary table (`read_stats.txt`) with IDs for each read.
#     - For each sample, appends haplotype assignments (H0/H1/H2) to the read ID in `read_stats`.
#     - Generates a new version of `sample_info` where each sample is duplicated per haplotype.
#     - Writes two outputs:
#         - `updated_read_stats.txt.gz`: updated read stats with haplotype-annotated IDs.
#         - `updated_sample_info.tsv.gz`: expanded sample metadata for each haplotype.

#     Parameters:
#     - sample_info_with_haplotype_location (str): Path to a TSV file containing columns `sample` and `haplotype`, where `haplotype` points to a file with read-level haplotype assignments.
#     - read_stats_path (str): Path to the original read stats file (`read_stats.txt`), with a column `id`.
#     - output_dir (str): Directory where updated output files will be written.

#     Returns:
#     - None: Writes updated read stats and sample info files to the specified output directory.
#     """
    

#     os.makedirs(output_dir, exist_ok=True)

#     # Load sample_info_with_haplotype_location.tsv
#     sample_info = pd.read_csv(sample_info_with_haplotype_location, sep="\t")

#     # Load read_stats.txt. Allow the file to be zipped or not
#     read_stats = pd.read_csv(read_stats_path, sep="\t", dtype={0: str}, compression='infer') # Ensure read id column (first column) is string

#     # Ensure the first column is named "id"
#     read_stats.rename(columns={read_stats.columns[0]: 'id'}, inplace=True)

#     # Processing sample_info to extract haplotype assignment files
#     all_haplotypes = []  # Store haplotype assignments in a list for efficient merging

#     for _, row in sample_info.iterrows():
#         sample = row["sample"]
#         haplotype_file = row["haplotype"]

#         print(f"Processing haplotype assignment file for {sample}", flush=True)


#         if pd.notna(haplotype_file) and haplotype_file.strip():
#             try:
#                 # Check if file is empty or contains only headers
#                 with open(haplotype_file, "r") as f:
#                     lines = f.readlines()
#                     if len(lines) <= 1:  # Only header row present
#                         print(f"Warning: {haplotype_file} contains only headers or is empty. Skipping.", flush=True)
#                         continue  # Skip this file

#                 # Read haplotype file in chunks if it is large
#                 for chunk in pd.read_csv(haplotype_file, sep="\t", chunksize=500000, usecols=["#readname", "haplotype"], dtype=str):
#                     if chunk.empty:
#                         continue  # Skip empty chunks
#                     chunk.rename(columns={"#readname": "id"}, inplace=True)
#                     all_haplotypes.append(chunk)

#             except FileNotFoundError:
#                 print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
#             except pd.errors.EmptyDataError:
#                 print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)

#     # Combine haplotype data into a single DataFrame
#     if all_haplotypes:
#         haplotype_df = pd.concat(all_haplotypes, ignore_index=True)
#     else:
#         haplotype_df = pd.DataFrame(columns=["id", "haplotype"])  # Empty DataFrame to avoid errors

#     # Merge haplotype info with read_stats
#     read_stats = read_stats.merge(haplotype_df, on="id", how="left")

#     # Replace NaN or "none" haplotypes with "H0"
#     read_stats["haplotype"] = read_stats["haplotype"].fillna("H0").replace("none", "H0")

#     # Modify `id` column efficiently using vectorized string operations
#     # read_stats["id"] = read_stats["id"].str.split("_", n=1).str[0] + read_stats["haplotype"] + "_" + read_stats["id"].str.split("_", n=1).str[1]
#     read_stats["id"] = read_stats["id"].str.split("_m", n=1).str[0] + read_stats["haplotype"] + "_m" + read_stats["id"].str.split("_m", n=1).str[1]

#     # Drop the `haplotype` column as it's no longer needed
#     read_stats.drop(columns=["haplotype"], inplace=True)

#     # Save updated read_stats
#     updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
#     read_stats.to_csv(updated_read_stats_path, sep="\t", index=False, compression = "gzip")
#     print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)


#     # Update sample_info with haplotype assignments
#     # updated_sample_info = []

#     # for _, row in sample_info.iterrows():
#     #     sample, individual, condition, haplotype = row["sample"], row["individual"], row["condition"], row["haplotype"]

#     #     if pd.notna(haplotype) and haplotype.strip():
#     #         updated_sample_info.append([f"{sample}H0", individual, condition, "H0"])
#     #         updated_sample_info.append([f"{sample}H1", individual, condition, "H1"])
#     #         updated_sample_info.append([f"{sample}H2", individual, condition, "H2"])
#     #     else:
#     #         updated_sample_info.append([f"{sample}H0", individual, condition, "H0"])

#     # Retain all columns from the original sample_info
#     updated_sample_info = []

#     for _, row in sample_info.iterrows():
#         row_dict = row.to_dict()  # Convert row to dictionary to preserve all columns
#         sample = row_dict["sample"]
#         haplotype = row_dict["haplotype"]

#         if pd.notna(haplotype) and haplotype.strip():
#             for hap in ["H0", "H1", "H2"]:
#                 new_row = row_dict.copy()
#                 new_row["sample"] = f"{sample}{hap}"
#                 new_row["haplotype"] = hap
#                 updated_sample_info.append(new_row)
#         else:
#             new_row = row_dict.copy()
#             new_row["sample"] = f"{sample}H0"
#             new_row["haplotype"] = "H0"
#             updated_sample_info.append(new_row)

#     # Convert to DataFrame and save
#     updated_sample_info_df = pd.DataFrame(updated_sample_info)
#     updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")
#     updated_sample_info_df.to_csv(updated_sample_info_path, index=False, compression = "gzip", sep="\t")

#     print(f"Updated sample info saved to {updated_sample_info_path}", flush=True)



def update_files_with_haplotype_info(
    sample_info_with_haplotype_location,
    read_stats_path,
    output_dir,
    read_stats_chunksize=500_000,
):
    os.makedirs(output_dir, exist_ok=True)

    sample_info = pd.read_csv(
        sample_info_with_haplotype_location,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # Build lightweight read_id -> haplotype dictionary
    haplotype_map = {}

    for row in sample_info.itertuples(index=False):
        sample = getattr(row, "sample")
        haplotype_file = getattr(row, "haplotype")

        if not haplotype_file or not haplotype_file.strip():
            continue

        print(f"Processing haplotype assignment file for {sample}", flush=True)

        if not os.path.exists(haplotype_file):
            print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
            continue

        try:
            for chunk in pd.read_csv(
                haplotype_file,
                sep="\t",
                usecols=["#readname", "haplotype"],
                dtype=str,
                chunksize=500_000,
            ):
                if chunk.empty:
                    continue

                chunk["haplotype"] = chunk["haplotype"].replace({"none": "H0"}).fillna("H0")

                haplotype_map.update(
                    zip(chunk["#readname"].to_numpy(), chunk["haplotype"].to_numpy())
                )

        except pd.errors.EmptyDataError:
            print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)
        except ValueError as e:
            print(f"Warning: Could not read {haplotype_file}: {e}. Skipping.", flush=True)

    updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")

    first_chunk = True

    with gzip.open(updated_read_stats_path, "wt") as out_f:
        for chunk in pd.read_csv(
            read_stats_path,
            sep="\t",
            dtype=str,
            compression="infer",
            chunksize=read_stats_chunksize,
            keep_default_na=False,
        ):
            # Rename first column to id
            if chunk.columns[0] != "id":
                chunk = chunk.rename(columns={chunk.columns[0]: "id"})

            hap = chunk["id"].map(haplotype_map).fillna("H0").replace("none", "H0")

            # Safer and faster than splitting twice
            parts = chunk["id"].str.split("_m", n=1, expand=True)

            has_m = parts.shape[1] == 2
            if has_m:
                chunk["id"] = parts[0] + hap + "_m" + parts[1]
            else:
                # Fallback if some IDs do not contain "_m"
                chunk["id"] = chunk["id"] + hap

            chunk.to_csv(
                out_f,
                sep="\t",
                index=False,
                header=first_chunk,
            )
            first_chunk = False

    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

    # Expand sample_info by haplotype without iterrows
    expanded_rows = []

    for row in sample_info.to_dict(orient="records"):
        sample = row["sample"]
        haplotype_file = row.get("haplotype", "")

        if haplotype_file and haplotype_file.strip():
            for hap in ["H0", "H1", "H2"]:
                new_row = row.copy()
                new_row["sample"] = f"{sample}{hap}"
                new_row["haplotype"] = hap
                expanded_rows.append(new_row)
        else:
            new_row = row.copy()
            new_row["sample"] = f"{sample}H0"
            new_row["haplotype"] = "H0"
            expanded_rows.append(new_row)

    updated_sample_info_df = pd.DataFrame(expanded_rows)

    updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")
    updated_sample_info_df.to_csv(
        updated_sample_info_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(f"Updated sample info saved to {updated_sample_info_path}", flush=True)



def filter_based_on_counts(df, count_threshold=10, group_col='Isoform'):
    """
    Filter isoforms based on count thresholds.
    
    Parameters:
    - df (pd.DataFrame): The input long-format DataFrame.
    - count_threshold (int): The threshold for filtering counts.
    - group_col (str): The column to group by (e.g., 'Isoform', 'Isoform_PBid').
    
    Returns:
    - pd.DataFrame: Filtered DataFrame with only the groups meeting the count threshold.
    """
    # Determine isoforms/groups to keep based on the threshold
    isoforms_to_keep = df.groupby(group_col).apply(
        lambda group: any(group['cyclo_count'] >= count_threshold) or any(group['noncyclo_count'] >= count_threshold),
        include_groups=False
    )
    isoforms_to_keep = isoforms_to_keep[isoforms_to_keep].index.tolist()

    # Filter the DataFrame to include only the isoforms/groups to keep
    return df[df[group_col].isin(isoforms_to_keep)]
