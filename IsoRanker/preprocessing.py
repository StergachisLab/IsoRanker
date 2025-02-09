import os
import pandas as pd

def update_files_with_haplotype_info(sample_info_with_haplotype_location, read_stats_path, output_dir):
    """
    Update read statistics and sample information files with haplotype data.
    
    This function performs the following tasks:
    1. Reads `sample_info_with_haplotype_location.csv` to extract sample and haplotype assignment file information.
    2. Reads `read_stats.txt` and updates read identifiers based on haplotype assignments.
    3. Generates an updated `sample_info.csv` with explicit haplotype assignments (`H0`, `H1`, `H2`).
    4. Saves the updated `read_stats.txt` and `updated_sample_info.csv` to the specified output directory.

    Parameters:
    - sample_info_with_haplotype_location (str): Path to the sample information CSV file (`sample_info_with_haplotype_location.csv`).
      - **Input Format:** This file must contain the columns:
        - `sample`: Sample identifier.
        - `patient`: Patient identifier.
        - `cyclo`: Cycloheximide treatment status.
        - `haplotype`: **File path to the haplotype assignment file** for each sample.
    - read_stats_path (str): Path to the read statistics TSV file (`read_stats.txt`).
      - **Input Format:** Must contain at least the column `id`, representing read identifiers.
    - output_dir (str): Path to the directory where updated files will be saved.

    Returns:
    - None: Saves `updated_read_stats.txt` and `updated_sample_info.csv` in the specified output directory.

    Output Files:
    - `updated_read_stats.txt` (TSV): Read statistics with modified `id` values incorporating haplotype assignments.
    - `updated_sample_info.csv` (CSV): Sample metadata with **haplotype values updated to `H0`, `H1`, and `H2`.**
      - **Output Format:** The `haplotype` column now contains discrete values (`H0`, `H1`, `H2`), rather than file paths.

    Raises:
    - FileNotFoundError: If any required file (`sample_info_with_haplotype_location.csv`, `read_stats.txt`, or haplotype files) is missing.
    - pd.errors.EmptyDataError: If any required file is empty.
    - KeyError: If required columns are missing from input files.
    """

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load sample_info_with_haplotype_location.csv
    sample_info = pd.read_csv(sample_info_with_haplotype_location)
    
    # Load read_stats.txt
    read_stats = pd.read_csv(read_stats_path, sep="\t")
    
    # Dictionary to store haplotype info for updating read_stats
    haplotype_dict = {}

    # Process sample_info to extract haplotype assignment files and update read_stats
    for _, row in sample_info.iterrows():
        sample = row["sample"]
        haplotype_file = row["haplotype"]
    
        print(f"Processing haplotype assignment file for {sample}", flush=True)
    
        if pd.notna(haplotype_file) and haplotype_file.strip():
            try:
                # Load haplotype assignment file (list.txt)
                list_df = pd.read_csv(haplotype_file, sep="\t")
    
                # Check if the DataFrame is empty
                if list_df.empty:
                    print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)
                    continue  # Skip to the next iteration
    
                # Ensure required columns exist before processing
                required_columns = {"#readname", "haplotype"}
                if not required_columns.issubset(list_df.columns):
                    print(f"Warning: {haplotype_file} is missing required columns. Skipping.", flush=True)
                    continue  # Skip processing this file
    
                # Map readname to haplotype
                for _, l_row in list_df.iterrows():
                    readname = l_row["#readname"]
                    haplotype = l_row["haplotype"]
                    haplotype_dict[readname] = haplotype
    
            except FileNotFoundError:
                print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
            except pd.errors.EmptyDataError:
                print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)

    def modify_id(row):
        read_id = row["id"]
        haplotype = haplotype_dict.get(read_id, "H0")  # Default to "H0" if not found
        return f"{read_id.split('_')[0]}{haplotype}_{'_'.join(read_id.split('_')[1:])}"

    # Update read_stats based on haplotype_dict
    read_stats["id"] = read_stats.apply(modify_id, axis=1)

    # Save updated read_stats
    updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt")
    read_stats.to_csv(updated_read_stats_path, sep="\t", index=False)
    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

    # Update sample_info with new sample names and haplotype assignments
    updated_sample_info = []

    for _, row in sample_info.iterrows():
        sample, patient, cyclo, haplotype = row["sample"], row["patient"], row["cyclo"], row["haplotype"]
        
        if pd.notna(haplotype) and haplotype.strip():
            updated_sample_info.append([f"{sample}H0", patient, cyclo, "H0"])
            updated_sample_info.append([f"{sample}H1", patient, cyclo, "H1"])
            updated_sample_info.append([f"{sample}H2", patient, cyclo, "H2"])
        else:
            updated_sample_info.append([f"{sample}H0", patient, cyclo, "H0"])

    # Convert to DataFrame and save
    updated_sample_info_df = pd.DataFrame(updated_sample_info, columns=["sample", "patient", "cyclo", "haplotype"])
    updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.csv")
    updated_sample_info_df.to_csv(updated_sample_info_path, index=False)

    print(f"Updated sample info saved to {updated_sample_info_path}", flush=True)


    # Update read_stats based on haplotype_dict
    def modify_id(row):
        read_id = row["id"]
        haplotype = haplotype_dict.get(read_id, "none")
        # If haplotype is "none", update it to "HP0"
        if haplotype == "none":
            haplotype = "H0"
        return f"{read_id.split('_')[0]}{haplotype}_{'_'.join(read_id.split('_')[1:])}"


    read_stats["id"] = read_stats.apply(modify_id, axis=1)

    # Save updated read_stats
    updated_read_stats_path = "updated_read_stats.txt"
    read_stats.to_csv(updated_read_stats_path, sep="\t", index=False)
    
    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

    # Update sample_info with new sample names
    updated_sample_info = []

    for _, row in sample_info.iterrows():
        sample, patient, cyclo, haplotype = row["sample"], row["patient"], row["cyclo"], row["haplotype"]
        
        if pd.notna(haplotype) and haplotype.strip():
            updated_sample_info.append([f"{sample}H0", patient, cyclo, "H0"])
            updated_sample_info.append([f"{sample}H1", patient, cyclo, "H1"])
            updated_sample_info.append([f"{sample}H2", patient, cyclo, "H2"])
        else:
            updated_sample_info.append([f"{sample}H0", patient, cyclo, "H0"])

    # Convert to DataFrame and save
    updated_sample_info_df = pd.DataFrame(updated_sample_info, columns=["sample", "patient", "cyclo", "haplotype"])
    updated_sample_info_path = "updated_sample_info.csv"
    updated_sample_info_df.to_csv(updated_sample_info_path, index=False)

    print(f"Updated sample_info saved to {updated_sample_info_path}", flush=True)



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
