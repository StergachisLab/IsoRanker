import pandas as pd
from collections import defaultdict

def parse_read_stats(file_path):
    """
    Parse the read_stats.txt file to extract sample and PB identifiers.

    Parameters:
    - file_path: Path to the read_stats.txt file.

    Returns:
    - A nested dictionary with PB identifiers as keys and sample counts as values.
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    with open(file_path, 'r') as f:
        for line in f:
            # Split the line into components
            read, pb_id = line.strip().split()
            
            # Extract the sample name (before the first "_")
            sample = read.split('_')[0]
            
            # Increment the count for the PB identifier in the sample
            counts[pb_id][sample] += 1
    
    return counts

def create_expression_matrix(file_path, output_file=None):
    """
    Create an expression matrix from a read_stats.txt file.

    Parameters:
    - file_path: Path to the read_stats.txt file.
    - output_file: Path to save the resulting expression matrix as a CSV (optional).

    Returns:
    - A pandas DataFrame representing the expression matrix.
    """
    # Parse the file and aggregate counts
    counts = parse_read_stats(file_path)
    
    # Create a DataFrame from the nested dictionary
    df = pd.DataFrame.from_dict(counts, orient='index').fillna(0).astype(int)
    
    # Sort rows and columns for better readability
    df.sort_index(inplace=True)
    df.sort_index(axis=1, inplace=True)
    
    # Save to a CSV file if specified
    if output_file:
        df.to_csv(output_file)
    
    return df


def create_long_format(expression_matrix, sample_info=None):
    """
    Create a long-format DataFrame where each isoform-Sample combination has only one row, including Cyclo_TPM and Noncyclo_TPM.
    
    Parameters:
    - expression_matrix (pd.DataFrame): Expression matrix with PB IDs as the index and sample names as columns.
    - sample_info (pd.DataFrame, optional): Sample info DataFrame with columns `sample`, `patient`, `cyclo`, and `haplotype`.
    
    Returns:
    - pd.DataFrame: Aggregated DataFrame with columns: Isoform, Sample, cyclo_count, noncyclo_count,
                    HP1_cyclo_count, HP2_cyclo_count, HP1_noncyclo_count, HP2_noncyclo_count, Cyclo_TPM, Noncyclo_TPM.
    """
    # Step 1: Handle the case where no sample_info is provided
    if sample_info is None:
        sample_info = pd.DataFrame({
            "sample": expression_matrix.columns,
            "patient": expression_matrix.columns,
            "cyclo": "noncyclo",  # Default to noncyclo if no info is provided
            "haplotype": None     # No haplotype information
        })

    # If the patient column is not provided, set it to the same values as the sample column
    if "patient" not in sample_info.columns:
        sample_info["patient"] = sample_info["sample"]

    # Step 2: Filter sample columns to keep only those found in both the expression matrix and sample info
    valid_samples = expression_matrix.columns.intersection(sample_info["sample"])
    sample_info = sample_info[sample_info["sample"].isin(valid_samples)]

    # Step 3: Melt the expression matrix into long format
    expression_data = expression_matrix[valid_samples].stack().reset_index()
    expression_data.columns = ["Isoform", "Sample", "count"]

    # Step 4: Merge with sample_info to enrich with metadata
    expression_data = expression_data.merge(sample_info, left_on="Sample", right_on="sample", how="left")

    # Step 5: Precompute haplotype-specific and cyclo/noncyclo counts
    expression_data["HP1_cyclo_count"] = (
        (expression_data["haplotype"] == "HP1") & (expression_data["cyclo"] == "cyclo")
    ) * expression_data["count"]
    expression_data["HP2_cyclo_count"] = (
        (expression_data["haplotype"] == "HP2") & (expression_data["cyclo"] == "cyclo")
    ) * expression_data["count"]
    expression_data["HP1_noncyclo_count"] = (
        (expression_data["haplotype"] == "HP1") & (expression_data["cyclo"] == "noncyclo")
    ) * expression_data["count"]
    expression_data["HP2_noncyclo_count"] = (
        (expression_data["haplotype"] == "HP2") & (expression_data["cyclo"] == "noncyclo")
    ) * expression_data["count"]

    # Precompute raw cyclo and noncyclo counts
    expression_data["cyclo_count_raw"] = (expression_data["cyclo"] == "cyclo") * expression_data["count"]
    expression_data["noncyclo_count_raw"] = (expression_data["cyclo"] == "noncyclo") * expression_data["count"]

    # Step 6: Aggregate counts by Isoform-Sample
    aggregated_data = expression_data.groupby(["Isoform", "patient"]).agg(
        HP1_cyclo_count=("HP1_cyclo_count", "sum"),
        HP2_cyclo_count=("HP2_cyclo_count", "sum"),
        HP1_noncyclo_count=("HP1_noncyclo_count", "sum"),
        HP2_noncyclo_count=("HP2_noncyclo_count", "sum"),
        cyclo_count=("cyclo_count_raw", "sum"),
        noncyclo_count=("noncyclo_count_raw", "sum")
    ).reset_index()

    # Rename the 'patient' column to 'Sample'
    aggregated_data.rename(columns={"patient": "Sample"}, inplace=True)
    
    # Step 7: Calculate total reads within each Sample
    sample_totals = aggregated_data.groupby("Sample")[["cyclo_count", "noncyclo_count"]].sum()
    sample_totals = sample_totals.rename(columns={"cyclo_count": "total_cyclo", "noncyclo_count": "total_noncyclo"})

    # Merge totals back to the aggregated data
    aggregated_data = aggregated_data.merge(sample_totals, on="Sample", how="left")

    # Avoid division by zero
    # aggregated_data["total_cyclo"] = aggregated_data["total_cyclo"].replace(0, 1)
    # aggregated_data["total_noncyclo"] = aggregated_data["total_noncyclo"].replace(0, 1)

    # Calculate Cyclo_TPM and Noncyclo_TPM
    aggregated_data["Cyclo_TPM"] = (aggregated_data["cyclo_count"] / aggregated_data["total_cyclo"]) * 1e6
    aggregated_data["Noncyclo_TPM"] = (aggregated_data["noncyclo_count"] / aggregated_data["total_noncyclo"]) * 1e6

    
    # Calculate Cyclo_TPM_rank and Noncyclo_TPM_rank with average ranking for ties. Should go from 1 to number of patients. The lower the rank, the larger the TPM.
    aggregated_data["Cyclo_TPM_Rank"] = aggregated_data.groupby("Isoform")["Cyclo_TPM"].rank(ascending=False, method="average")
    aggregated_data["Noncyclo_TPM_Rank"] = aggregated_data.groupby("Isoform")["Noncyclo_TPM"].rank(ascending=False, method="average")


    # Step 8: Drop unnecessary columns.
    # If the haplotype column are all empty or NaN, then we are not evalauting for haplotype separated information so these columns can be dropped.
    if sample_info['haplotype'].replace('', float('NaN'), inplace=False).isna().all():
        aggregated_data = aggregated_data.drop(columns=["HP1_cyclo_count", "HP2_cyclo_count", "HP1_noncyclo_count", "HP2_noncyclo_count"])

    # Step 9: Return the aggregated DataFrame
    return aggregated_data
