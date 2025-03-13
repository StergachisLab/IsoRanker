import os
import pandas as pd

def update_files_with_haplotype_info(sample_info_with_haplotype_location, read_stats_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Load sample_info_with_haplotype_location.csv
    sample_info = pd.read_csv(sample_info_with_haplotype_location)

    # Load read_stats.txt
    read_stats = pd.read_csv(read_stats_path, sep="\t", dtype={"id": str})  # Ensure id column is string

    # Processing sample_info to extract haplotype assignment files
    all_haplotypes = []  # Store haplotype assignments in a list for efficient merging

    for _, row in sample_info.iterrows():
        sample = row["sample"]
        haplotype_file = row["haplotype"]

        print(f"Processing haplotype assignment file for {sample}", flush=True)

        if pd.notna(haplotype_file) and haplotype_file.strip():
            try:
                # Read haplotype file in chunks if it is large
                for chunk in pd.read_csv(haplotype_file, sep="\t", chunksize=500000, usecols=["#readname", "haplotype"], dtype=str):
                    chunk.rename(columns={"#readname": "id"}, inplace=True)
                    all_haplotypes.append(chunk)
            except FileNotFoundError:
                print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
            except pd.errors.EmptyDataError:
                print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)

    # Combine haplotype data into a single DataFrame
    if all_haplotypes:
        haplotype_df = pd.concat(all_haplotypes, ignore_index=True)
    else:
        haplotype_df = pd.DataFrame(columns=["id", "haplotype"])  # Empty DataFrame to avoid errors

    # Merge haplotype info with read_stats
    read_stats = read_stats.merge(haplotype_df, on="id", how="left")

    # Replace NaN or "none" haplotypes with "H0"
    read_stats["haplotype"] = read_stats["haplotype"].fillna("H0").replace("none", "H0")

    # Modify `id` column efficiently using vectorized string operations
    # read_stats["id"] = read_stats["id"].str.split("_", n=1).str[0] + read_stats["haplotype"] + "_" + read_stats["id"].str.split("_", n=1).str[1]
    read_stats["id"] = read_stats["id"].str.split("_m", n=1).str[0] + read_stats["haplotype"] + "_m" + read_stats["id"].str.split("_m", n=1).str[1]

    # Drop the `haplotype` column as it's no longer needed
    read_stats.drop(columns=["haplotype"], inplace=True)

    # Save updated read_stats
    updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
    read_stats.to_csv(updated_read_stats_path, sep="\t", index=False, compression = "gzip")
    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

    # Update sample_info with haplotype assignments
    updated_sample_info = []

    for _, row in sample_info.iterrows():
        identifier, individual, condition, haplotype = row["sample"], row["individual"], row["condition"], row["haplotype"]

        if pd.notna(haplotype) and haplotype.strip():
            updated_sample_info.append([f"{sample}H0", individual, condition, "H0"])
            updated_sample_info.append([f"{sample}H1", individual, condition, "H1"])
            updated_sample_info.append([f"{sample}H2", individual, condition, "H2"])
        else:
            updated_sample_info.append([f"{sample}H0", individual, condition, "H0"])

    # Convert to DataFrame and save
    updated_sample_info_df = pd.DataFrame(updated_sample_info, columns=["sample", "individual", "condition", "haplotype"])
    updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.csv.gz")
    updated_sample_info_df.to_csv(updated_sample_info_path, index=False, compression = "gzip")

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
        lambda group: any(group['cyclo_count'] >= count_threshold) or any(group['noncyclo_count'] >= count_threshold)
    )
    isoforms_to_keep = isoforms_to_keep[isoforms_to_keep].index.tolist()

    # Filter the DataFrame to include only the isoforms/groups to keep
    return df[df[group_col].isin(isoforms_to_keep)]
