import pandas as pd
from collections import defaultdict
import gzip
import numpy as np



def parse_read_stats(file_path, samples_to_keep=None):
    """
    Parse the read_stats.txt or read_stats.txt.gz file to extract sample and PB identifiers.

    Parameters:
    - file_path (str): Path to the read_stats.txt or read_stats.txt.gz file.

    Returns:
    - dict: A nested dictionary with PB identifiers as keys and sample counts as values.
    """
    counts = defaultdict(lambda: defaultdict(int))

    open_func = gzip.open if file_path.endswith(".gz") else open

    with open_func(file_path, "rt") as f:
        for line in f:
            read, pb_id = line.strip().split()

            sample = read.split("_m")[0]

            if samples_to_keep is not None and sample not in samples_to_keep:
                continue

            counts[pb_id][sample] += 1

    return counts


def create_expression_matrix(file_path, output_file=None, sample_info=None):
    """
    Create an expression matrix from a read_stats.txt/read_stats.txt.gz file.

    If sample_info is provided, only samples listed in sample_info["sample"]
    are included. If sample_info is not provided, all samples in read_stats
    are included.


    Parameters:
    - file_path: Path to the read_stats.txt file.
    - output_file: Path to save the resulting expression matrix as a CSV (optional).

    Returns:
    - A pandas DataFrame representing the expression matrix.
    """

    samples_to_keep = None

    if sample_info is not None:
        if isinstance(sample_info, str):
            sample_info = pd.read_csv(sample_info, sep="\t")

        if "sample" not in sample_info.columns:
            raise ValueError("sample_info must contain a 'sample' column")

        samples_to_keep = set(sample_info["sample"].astype(str))

    counts = parse_read_stats(file_path, samples_to_keep=samples_to_keep)

    df = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)

    df.sort_index(inplace=True)
    df.sort_index(axis=1, inplace=True)

    if output_file:
        df.to_csv(output_file, index=True, compression="gzip", sep="\t")

    return df




def create_long_format(expression_matrix, sample_info=None):
    """
    Create a long-format DataFrame where each isoform-Sample combination has only one row, including Cyclo_TPM and Noncyclo_TPM.
    
    Parameters:
    - expression_matrix (pd.DataFrame): Expression matrix with PB IDs as the index and sample names as columns.
    - sample_info (pd.DataFrame, optional): Sample info DataFrame with columns `sample`, `patient`, `cyclo`, and `haplotype`.
    
    Returns:
    - pd.DataFrame: Aggregated DataFrame with columns: Isoform, Sample, cyclo_count, noncyclo_count,
                    H1_cyclo_count, H2_cyclo_count, H1_noncyclo_count, H2_noncyclo_count, Cyclo_TPM, Noncyclo_TPM.
    """

    expr = expression_matrix

    if sample_info is None:
        sample_info = pd.DataFrame({
            "sample": expr.columns,
            "individual": expr.columns,
            "condition": "noncyclo",
            "haplotype": "H0",
        }).copy()
    else:
        sample_info = sample_info.copy()

    if "individual" not in sample_info.columns:
        sample_info["individual"] = sample_info["sample"]

    if "haplotype" not in sample_info.columns:
        sample_info["haplotype"] = "H0"

    if "condition" not in sample_info.columns:
        sample_info["condition"] = "noncyclo"

    valid_samples = expr.columns.intersection(sample_info["sample"])
    sample_info = sample_info[sample_info["sample"].isin(valid_samples)].copy()

    sample_info["condition"] = sample_info["condition"].astype(str).str.lower()
    sample_info["haplotype"] = sample_info["haplotype"].fillna("H0").astype(str)

    individuals = pd.Index(sample_info["individual"].drop_duplicates())
    isoforms = expr.index

    n_iso = len(isoforms)
    n_ind = len(individuals)

    ind_to_col = {ind: i for i, ind in enumerate(individuals)}

    dtype = np.float32 if np.issubdtype(expr.dtypes.iloc[0], np.floating) else np.int64

    H0_cyclo = np.zeros((n_iso, n_ind), dtype=dtype)
    H1_cyclo = np.zeros((n_iso, n_ind), dtype=dtype)
    H2_cyclo = np.zeros((n_iso, n_ind), dtype=dtype)

    H0_noncyclo = np.zeros((n_iso, n_ind), dtype=dtype)
    H1_noncyclo = np.zeros((n_iso, n_ind), dtype=dtype)
    H2_noncyclo = np.zeros((n_iso, n_ind), dtype=dtype)

    cyclo = np.zeros((n_iso, n_ind), dtype=dtype)
    noncyclo = np.zeros((n_iso, n_ind), dtype=dtype)

    for row in sample_info.itertuples(index=False):
        sample = row.sample
        individual = row.individual
        condition = row.condition
        haplotype = row.haplotype

        j = ind_to_col[individual]
        values = expr[sample].to_numpy(dtype=dtype, copy=False)

        if condition == "cyclo":
            cyclo[:, j] += values

            if haplotype == "H0":
                H0_cyclo[:, j] += values
            elif haplotype == "H1":
                H1_cyclo[:, j] += values
            elif haplotype == "H2":
                H2_cyclo[:, j] += values

        elif condition == "noncyclo":
            noncyclo[:, j] += values

            if haplotype == "H0":
                H0_noncyclo[:, j] += values
            elif haplotype == "H1":
                H1_noncyclo[:, j] += values
            elif haplotype == "H2":
                H2_noncyclo[:, j] += values

    index = pd.MultiIndex.from_product(
        [isoforms, individuals],
        names=["Isoform", "Sample"]
    )

    aggregated_data = pd.DataFrame({
        "cyclo_count": cyclo.ravel(),
        "noncyclo_count": noncyclo.ravel(),
        "H0_cyclo_count": H0_cyclo.ravel(),
        "H1_cyclo_count": H1_cyclo.ravel(),
        "H2_cyclo_count": H2_cyclo.ravel(),
        "H0_noncyclo_count": H0_noncyclo.ravel(),
        "H1_noncyclo_count": H1_noncyclo.ravel(),
        "H2_noncyclo_count": H2_noncyclo.ravel(),
    }, index=index).reset_index()

    total_cyclo = aggregated_data.groupby("Sample", sort=False)["cyclo_count"].transform("sum")
    total_noncyclo = aggregated_data.groupby("Sample", sort=False)["noncyclo_count"].transform("sum")

    aggregated_data["total_cyclo"] = total_cyclo
    aggregated_data["total_noncyclo"] = total_noncyclo

    aggregated_data["Cyclo_TPM"] = np.where(
        total_cyclo > 0,
        aggregated_data["cyclo_count"] / total_cyclo * 1e6,
        0,
    )

    aggregated_data["Noncyclo_TPM"] = np.where(
        total_noncyclo > 0,
        aggregated_data["noncyclo_count"] / total_noncyclo * 1e6,
        0,
    )

    aggregated_data["Cyclo_TPM_Rank"] = (
        aggregated_data.groupby("Isoform", sort=False)["Cyclo_TPM"]
        .rank(ascending=False, method="average")
    )

    aggregated_data["Noncyclo_TPM_Rank"] = (
        aggregated_data.groupby("Isoform", sort=False)["Noncyclo_TPM"]
        .rank(ascending=False, method="average")
    )

    grouped = aggregated_data.groupby("Isoform", sort=False)

    aggregated_data["Cyclo_TPM_Median"] = grouped["Cyclo_TPM"].transform("median")
    aggregated_data["Noncyclo_TPM_Median"] = grouped["Noncyclo_TPM"].transform("median")
    aggregated_data["Cyclo_TPM_Mean"] = grouped["Cyclo_TPM"].transform("mean")
    aggregated_data["Noncyclo_TPM_Mean"] = grouped["Noncyclo_TPM"].transform("mean")

    hap = sample_info["haplotype"]
    no_haplotype_info = (hap.isna() | hap.isin(["", "none", "H0"])).all()

    if no_haplotype_info:
        aggregated_data = aggregated_data.drop(columns=[
            "H0_cyclo_count", "H1_cyclo_count", "H2_cyclo_count",
            "H0_noncyclo_count", "H1_noncyclo_count", "H2_noncyclo_count",
        ])

    return aggregated_data



# def create_long_format(expression_matrix, sample_info=None):
#     """
#     Create a long-format DataFrame where each isoform-Sample combination has only one row, including Cyclo_TPM and Noncyclo_TPM.
    
#     Parameters:
#     - expression_matrix (pd.DataFrame): Expression matrix with PB IDs as the index and sample names as columns.
#     - sample_info (pd.DataFrame, optional): Sample info DataFrame with columns `sample`, `patient`, `cyclo`, and `haplotype`.
    
#     Returns:
#     - pd.DataFrame: Aggregated DataFrame with columns: Isoform, Sample, cyclo_count, noncyclo_count,
#                     H1_cyclo_count, H2_cyclo_count, H1_noncyclo_count, H2_noncyclo_count, Cyclo_TPM, Noncyclo_TPM.
#     """
#     # Step 1: Handle the case where no sample_info is provided
#     if sample_info is None:
#         sample_info = pd.DataFrame({
#             "sample": expression_matrix.columns,
#             "individual": expression_matrix.columns,
#             "condition": "noncyclo",  # Default to noncyclo if no info is provided
#             "haplotype": "H0"     # No haplotype information
#         })

#     # If the patient column is not provided, set it to the same values as the sample column
#     if "individual" not in sample_info.columns:
#         sample_info["individual"] = sample_info["sample"]

#     # If the haplotype column is not provided, set it to none. 
#     if "haplotype" not in sample_info.columns:
#         sample_info["haplotype"] = "H0"

#     # Step 2: Filter sample columns to keep only those found in both the expression matrix and sample info
#     valid_samples = expression_matrix.columns.intersection(sample_info["sample"])
#     sample_info = sample_info[sample_info["sample"].isin(valid_samples)]

#     # Step 3: Melt the expression matrix into long format
#     expression_data = expression_matrix[valid_samples].stack().reset_index()
#     expression_data.columns = ["Isoform", "Sample", "count"]

#     # Step 4: Merge with sample_info to enrich with metadata
#     expression_data = expression_data.merge(sample_info, left_on="Sample", right_on="sample", how="left")

#     # Step 5: Precompute haplotype-specific and cyclo/noncyclo counts
#     expression_data["H0_cyclo_count"] = (
#         (expression_data["haplotype"] == "H0") & (expression_data["condition"] == "cyclo")
#     ) * expression_data["count"]
#     expression_data["H1_cyclo_count"] = (
#         (expression_data["haplotype"] == "H1") & (expression_data["condition"] == "cyclo")
#     ) * expression_data["count"]
#     expression_data["H2_cyclo_count"] = (
#         (expression_data["haplotype"] == "H2") & (expression_data["condition"] == "cyclo")
#     ) * expression_data["count"]
#     expression_data["H0_noncyclo_count"] = (
#         (expression_data["haplotype"] == "H0") & (expression_data["condition"] == "noncyclo")
#     ) * expression_data["count"]
#     expression_data["H1_noncyclo_count"] = (
#         (expression_data["haplotype"] == "H1") & (expression_data["condition"] == "noncyclo")
#     ) * expression_data["count"]
#     expression_data["H2_noncyclo_count"] = (
#         (expression_data["haplotype"] == "H2") & (expression_data["condition"] == "noncyclo")
#     ) * expression_data["count"]

#     # Precompute raw cyclo and noncyclo counts
#     expression_data["cyclo_count_raw"] = (expression_data["condition"] == "cyclo") * expression_data["count"]
#     expression_data["noncyclo_count_raw"] = (expression_data["condition"] == "noncyclo") * expression_data["count"]

#     # Step 6: Aggregate counts by Isoform-Sample
#     aggregated_data = expression_data.groupby(["Isoform", "individual"]).agg(
#         H0_cyclo_count=("H0_cyclo_count", "sum"),
#         H1_cyclo_count=("H1_cyclo_count", "sum"),
#         H2_cyclo_count=("H2_cyclo_count", "sum"),
#         H0_noncyclo_count=("H0_noncyclo_count", "sum"),
#         H1_noncyclo_count=("H1_noncyclo_count", "sum"),
#         H2_noncyclo_count=("H2_noncyclo_count", "sum"),
#         cyclo_count=("cyclo_count_raw", "sum"),
#         noncyclo_count=("noncyclo_count_raw", "sum")
#     ).reset_index()

#     # Rename the 'patient' column to 'Sample'
#     aggregated_data.rename(columns={"individual": "Sample"}, inplace=True)
    
#     # Step 7: Calculate total reads within each Sample
#     sample_totals = aggregated_data.groupby("Sample")[["cyclo_count", "noncyclo_count"]].sum()
#     sample_totals = sample_totals.rename(columns={"cyclo_count": "total_cyclo", "noncyclo_count": "total_noncyclo"})

#     # Merge totals back to the aggregated data
#     aggregated_data = aggregated_data.merge(sample_totals, on="Sample", how="left")

#     # Avoid division by zero
#     # aggregated_data["total_cyclo"] = aggregated_data["total_cyclo"].replace(0, 1)
#     # aggregated_data["total_noncyclo"] = aggregated_data["total_noncyclo"].replace(0, 1)

#     # Calculate Cyclo_TPM and Noncyclo_TPM
#     aggregated_data["Cyclo_TPM"] = (aggregated_data["cyclo_count"] / aggregated_data["total_cyclo"]) * 1e6
#     aggregated_data["Noncyclo_TPM"] = (aggregated_data["noncyclo_count"] / aggregated_data["total_noncyclo"]) * 1e6

    
#     # Calculate Cyclo_TPM_rank and Noncyclo_TPM_rank with average ranking for ties. Should go from 1 to number of patients. The lower the rank, the larger the TPM.
#     aggregated_data["Cyclo_TPM_Rank"] = aggregated_data.groupby("Isoform")["Cyclo_TPM"].rank(ascending=False, method="average")
#     aggregated_data["Noncyclo_TPM_Rank"] = aggregated_data.groupby("Isoform")["Noncyclo_TPM"].rank(ascending=False, method="average")

#     # Calculate mean and median TPM for each Isoform
#     aggregated_data["Cyclo_TPM_Median"] = aggregated_data.groupby("Isoform")["Cyclo_TPM"].transform("median")
#     aggregated_data["Noncyclo_TPM_Median"] = aggregated_data.groupby("Isoform")["Noncyclo_TPM"].transform("median")
#     aggregated_data["Cyclo_TPM_Mean"] = aggregated_data.groupby("Isoform")["Cyclo_TPM"].transform("mean")
#     aggregated_data["Noncyclo_TPM_Mean"] = aggregated_data.groupby("Isoform")["Noncyclo_TPM"].transform("mean")


#     # Step 8: Drop unnecessary columns.
#     # If the haplotype column are all empty or NaN, or "none", or "H0", then we are not evalauting for haplotype separated information so these columns can be dropped.

#     hap = sample_info['haplotype']
#     no_haplotype_info = (hap.isna() | hap.isin(['', 'none', 'H0'])).all()

#     if no_haplotype_info:
#         aggregated_data = aggregated_data.drop(
#             columns=[
#                 "H0_cyclo_count", "H1_cyclo_count", "H2_cyclo_count",
#                 "H0_noncyclo_count", "H1_noncyclo_count", "H2_noncyclo_count",
#             ]
#         )

#     # Step 9: Return the aggregated DataFrame
#     return aggregated_data
