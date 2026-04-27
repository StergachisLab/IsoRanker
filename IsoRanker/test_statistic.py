import numpy as np

from .preprocessing import filter_based_on_counts

from .calculations import apply_hypothesis_test

from .z_score import calculate_z_score

from .z_score import calculate_z_score_MAD

from .ranking import calculate_ranks_for_sample

def _leave_one_out_min(s):
    vals = s.to_numpy()
    n = len(vals)
    if n <= 1:
        return np.full(n, np.nan)

    min1 = np.nanmin(vals)
    min_count = np.sum(vals == min1)

    if min_count > 1:
        return np.full(n, min1)

    min2 = np.nanmin(vals[vals != min1])
    return np.where(vals == min1, min2, min1)


def _leave_one_out_max(s):
    vals = s.to_numpy()
    n = len(vals)
    if n <= 1:
        return np.full(n, np.nan)

    max1 = np.nanmax(vals)
    max_count = np.sum(vals == max1)

    if max_count > 1:
        return np.full(n, max1)

    max2 = np.nanmax(vals[vals != max1])
    return np.where(vals == max1, max2, max1)


def _leave_one_out_median(s):
    vals = s.to_numpy()
    n = len(vals)
    if n <= 1:
        return np.full(n, np.nan)

    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(np.delete(vals, i))
    return out

# def NMD_test_statistic(group):
#     """
#     Calculate test statistic for NMD.
#     This function calculates a unique test statistic for each sample in the group.
#     """
#     results = []
#     for _, row in group.iterrows():
#         # Calculate the test statistic for this specific sample
#         ratio = (row['Cyclo_TPM'] + 1) / (row['Noncyclo_TPM'] + 1)
#         test_statistic = np.log2(ratio) * np.log2(row['Cyclo_TPM'] + 2)
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def NMD_test_statistic(group):
    """
    Calculate test statistic for NMD.
    This function calculates a unique test statistic for each sample in the group.
    """
    group = group.copy()

    ratio = (group["Cyclo_TPM"] + 1) / (group["Noncyclo_TPM"] + 1)
    group["test_statistic"] = np.log2(ratio) * np.log2(group["Cyclo_TPM"] + 2)

    return group

def NMD_hap1_test_statistic(group):
    """
    Calculate test statistic for NMD.
    This function calculates a unique test statistic for each sample in the group.
    """
    results = []
    for _, row in group.iterrows():
        # Calculate the test statistic for this specific sample

        # Adequate phasing for noncyclo
        noncyclo_reads_phased = (row['H1_noncyclo_count'] + row['H2_noncyclo_count'])
        noncyclo_proportion_phased = (noncyclo_reads_phased + 1) / (noncyclo_reads_phased + row['H0_noncyclo_count'] + 1)

        # Adequate phasing for cyclo
        cyclo_reads_phased = (row['H1_cyclo_count'] + row['H2_cyclo_count'])
        cyclo_proportion_phased = (cyclo_reads_phased + 1) / (cyclo_reads_phased + row['H0_cyclo_count'] + 1)

        # Convert to TPM
        row_Cyclo_hap1_TPM = (row['H1_cyclo_count'] / row["total_cyclo"]) * 1e6
        row_Noncyclo_hap1_TPM = (row['H1_noncyclo_count'] / row["total_noncyclo"]) * 1e6

        #If the proportion phased is too small, then the test stat can't be calculated.
        if (noncyclo_proportion_phased < 0.1) | (noncyclo_reads_phased < 10) | (cyclo_proportion_phased < 0.1) | (cyclo_reads_phased < 10):
            test_statistic = 0
        else:
            # Calculate the test statistic for this specific sample
            ratio = (row_Cyclo_hap1_TPM + 1) / (row_Noncyclo_hap1_TPM + 1)
            test_statistic = np.log2(ratio) * np.log2(row_Cyclo_hap1_TPM + 2)

        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group

def NMD_hap2_test_statistic(group):
    """
    Calculate test statistic for NMD.
    This function calculates a unique test statistic for each sample in the group.
    """
    results = []
    for _, row in group.iterrows():
        # Calculate the test statistic for this specific sample

        # Adequate phasing for noncyclo
        noncyclo_reads_phased = (row['H1_noncyclo_count'] + row['H2_noncyclo_count'])
        noncyclo_proportion_phased = (noncyclo_reads_phased + 1) / (noncyclo_reads_phased + row['H0_noncyclo_count'] + 1)

        # Adequate phasing for cyclo
        cyclo_reads_phased = (row['H1_cyclo_count'] + row['H2_cyclo_count'])
        cyclo_proportion_phased = (cyclo_reads_phased + 1) / (cyclo_reads_phased + row['H0_cyclo_count'] + 1)

        # Convert to TPM
        row_Cyclo_hap2_TPM = (row['H2_cyclo_count'] / row["total_cyclo"]) * 1e6
        row_Noncyclo_hap2_TPM = (row['H2_noncyclo_count'] / row["total_noncyclo"]) * 1e6

        #If the proportion phased is too small, then the test stat can't be calculated.
        if (noncyclo_proportion_phased < 0.1) | (noncyclo_reads_phased < 10) | (cyclo_proportion_phased < 0.1) | (cyclo_reads_phased < 10):
            test_statistic = 0
        else:
            # Calculate the test statistic for this specific sample
            ratio = (row_Cyclo_hap2_TPM + 1) / (row_Noncyclo_hap2_TPM + 1)
            test_statistic = np.log2(ratio) * np.log2(row_Cyclo_hap2_TPM + 2)

        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group

# def Noncyclo_Expression_Outlier_LOE(group):
#     """Calculate test statistic for Noncyclo Expression Outlier - Loss of Expression (LOE)."""
#     results = []
#     for _, row in group.iterrows():
#         # Exclude the current sample and calculate metrics from other samples
#         other_samples = group.drop(row.name)
#         min_noncyclo_TPM = other_samples['Noncyclo_TPM'].min()
#         median_noncyclo_TPM = other_samples['Noncyclo_TPM'].median()
        
#         # Calculate the test statistic
#         ratio = (min_noncyclo_TPM + 1) / (row['Noncyclo_TPM'] + 1)
#         test_statistic = np.log2(ratio) * np.log2(median_noncyclo_TPM + 2)
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def Noncyclo_Expression_Outlier_LOE(group):
    """Calculate test statistic for Noncyclo Expression Outlier - Loss of Expression (LOE)."""
    group = group.copy()

    loo_min = _leave_one_out_min(group["Noncyclo_TPM"])
    loo_median = _leave_one_out_median(group["Noncyclo_TPM"])

    ratio = (loo_min + 1) / (group["Noncyclo_TPM"] + 1)
    group["test_statistic"] = np.log2(ratio) * np.log2(loo_median + 2)

    return group

# def Noncyclo_Expression_Outlier_GOE(group):
#     """Calculate test statistic for Noncyclo Expression Outlier - Gain of Expression (GOE)."""
#     results = []
#     for _, row in group.iterrows():
#         # Exclude the current sample and calculate metrics from other samples
#         other_samples = group.drop(row.name)
#         max_noncyclo_TPM = other_samples['Noncyclo_TPM'].max()
#         median_noncyclo_TPM = other_samples['Noncyclo_TPM'].median()
        
#         # Calculate the test statistic
#         ratio = (row['Noncyclo_TPM'] + 1) / (max_noncyclo_TPM + 1)
#         test_statistic = np.log2(ratio) * np.log2(median_noncyclo_TPM + 2)
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def Noncyclo_Expression_Outlier_GOE(group):
    """Calculate test statistic for Noncyclo Expression Outlier - Gain of Expression (GOE)."""
    group = group.copy()

    loo_max = _leave_one_out_max(group["Noncyclo_TPM"])
    loo_median = _leave_one_out_median(group["Noncyclo_TPM"])

    ratio = (group["Noncyclo_TPM"] + 1) / (loo_max + 1)
    group["test_statistic"] = np.log2(ratio) * np.log2(loo_median + 2)

    return group

# def Cyclo_Expression_Outlier_LOE(group):
#     """Calculate test statistic for Cyclo Expression Outlier - Loss of Expression (LOE)."""
#     results = []
#     for _, row in group.iterrows():
#         # Exclude the current sample and calculate metrics from other samples
#         other_samples = group.drop(row.name)
#         min_cyclo_TPM = other_samples['Cyclo_TPM'].min()
#         median_cyclo_TPM = other_samples['Cyclo_TPM'].median()
        
#         # Calculate the test statistic
#         ratio = (min_cyclo_TPM + 1) / (row['Cyclo_TPM'] + 1)
#         test_statistic = np.log2(ratio) * np.log2(median_cyclo_TPM + 2)
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group


def Cyclo_Expression_Outlier_LOE(group):
    """Calculate test statistic for Cyclo Expression Outlier - Loss of Expression (LOE)."""
    group = group.copy()

    loo_min = _leave_one_out_min(group["Cyclo_TPM"])
    loo_median = _leave_one_out_median(group["Cyclo_TPM"])

    ratio = (loo_min + 1) / (group["Cyclo_TPM"] + 1)
    group["test_statistic"] = np.log2(ratio) * np.log2(loo_median + 2)

    return group

# def Cyclo_Expression_Outlier_GOE(group):
#     """Calculate test statistic for Cyclo Expression Outlier - Gain of Expression (GOE)."""
#     results = []
#     for _, row in group.iterrows():
#         # Exclude the current sample and calculate metrics from other samples
#         other_samples = group.drop(row.name)
#         max_cyclo_TPM = other_samples['Cyclo_TPM'].max()
#         median_cyclo_TPM = other_samples['Cyclo_TPM'].median()
        
#         # Calculate the test statistic
#         ratio = (row['Cyclo_TPM'] + 1) / (max_cyclo_TPM + 1)
#         test_statistic = np.log2(ratio) * np.log2(median_cyclo_TPM + 2)
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def Cyclo_Expression_Outlier_GOE(group):
    """Calculate test statistic for Cyclo Expression Outlier - Gain of Expression (GOE)."""
    group = group.copy()

    loo_max = _leave_one_out_max(group["Cyclo_TPM"])
    loo_median = _leave_one_out_median(group["Cyclo_TPM"])

    ratio = (group["Cyclo_TPM"] + 1) / (loo_max + 1)
    group["test_statistic"] = np.log2(ratio) * np.log2(loo_median + 2)

    return group


# def NMD_rare_steady_state_transcript(group):
#     """
#     Calculate test statistic for NMD Rare Steady State Transcript.
#     This function calculates a unique test statistic for each sample in the group.
#     """
#     results = []
#     for _, row in group.iterrows():
#         # Calculate the test statistic for this specific sample
#         cyclo_term = row['Total_bin_cyclo_count_Bin1_le'] / (row['cyclo_count'] + 1)
#         noncyclo_term = row['Total_bin_noncyclo_count_Bin1_le'] / (row['noncyclo_count'] + 1)
#         test_statistic = (
#             (cyclo_term - noncyclo_term) *
#             np.log2(row['Total_bin_cyclo_count_Bin1_le'] + 1) *
#             np.log2(row['Total_bin_noncyclo_count_Bin2_g'] + 1)
#         )
#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group


def NMD_rare_steady_state_transcript(group):
    """
    Calculate test statistic for NMD Rare Steady State Transcript.
    This function calculates a unique test statistic for each sample in the group.
    """
    group = group.copy()

    cyclo_term = group["Total_bin_cyclo_count_Bin1_le"] / (group["cyclo_count"] + 1)
    noncyclo_term = group["Total_bin_noncyclo_count_Bin1_le"] / (group["noncyclo_count"] + 1)

    group["test_statistic"] = (
        (cyclo_term - noncyclo_term) *
        np.log2(group["Total_bin_cyclo_count_Bin1_le"] + 1) *
        np.log2(group["Total_bin_noncyclo_count_Bin2_g"] + 1)
    )

    return group

def NMD_rare_steady_state_transcript_hap1(group):
    """
    Calculate test statistic for NMD Rare Steady State Transcript (Haplotype 1).
    Assumes the bin-aggregated columns (Total_bin_*) were computed using H1_* counts.
    """
    results = []
    for _, row in group.iterrows():
        cyclo_term = row["Total_bin_cyclo_count_Bin1_le"] / (row["H1_cyclo_count"] + 1)
        noncyclo_term = row["Total_bin_noncyclo_count_Bin1_le"] / (row["H1_noncyclo_count"] + 1)
        test_statistic = (
            (cyclo_term - noncyclo_term)
            * np.log2(row["Total_bin_cyclo_count_Bin1_le"] + 1)
            * np.log2(row["Total_bin_noncyclo_count_Bin2_g"] + 1)
        )
        results.append(test_statistic)

    group["test_statistic"] = results
    return group


def NMD_rare_steady_state_transcript_hap2(group):
    """
    Calculate test statistic for NMD Rare Steady State Transcript (Haplotype 2).
    Assumes the bin-aggregated columns (Total_bin_*) were computed using H2_* counts.
    """
    results = []
    for _, row in group.iterrows():
        cyclo_term = row["Total_bin_cyclo_count_Bin1_le"] / (row["H2_cyclo_count"] + 1)
        noncyclo_term = row["Total_bin_noncyclo_count_Bin1_le"] / (row["H2_noncyclo_count"] + 1)
        test_statistic = (
            (cyclo_term - noncyclo_term)
            * np.log2(row["Total_bin_cyclo_count_Bin1_le"] + 1)
            * np.log2(row["Total_bin_noncyclo_count_Bin2_g"] + 1)
        )
        results.append(test_statistic)

    group["test_statistic"] = results
    return group


# def Noncyclo_Allelic_Imbalance(group):
#     """
#     Calculate test statistic for Noncyclo allelic imbalance.
#     This function calculates a unique test statistic for each sample in the group.
#     """
#     results = []
#     for _, row in group.iterrows():
#         # Calculate the test statistic for this specific sample

#         reads_phased = (row['H1_noncyclo_count'] + row['H2_noncyclo_count'])
#         proportion_phased = (reads_phased + 1) / (reads_phased + row['H0_noncyclo_count'] + 1)

#         #If the proportion phased is too small, then the test stat can't be calculated.
#         if (proportion_phased < 0.1) | (reads_phased < 10):
#             test_statistic = 0
#         else:
#             test_statistic = abs(
#                 np.log2((row['H1_noncyclo_count']+1) / (row['H2_noncyclo_count']+1)) *
#                 np.log2(row['H1_noncyclo_count'] + row['H2_noncyclo_count'])
#             )

#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def Noncyclo_Allelic_Imbalance(group):
    """
    Calculate test statistic for Noncyclo allelic imbalance.
    This function calculates a unique test statistic for each sample in the group.
    """
    group = group.copy()

    reads_phased = group["H1_noncyclo_count"] + group["H2_noncyclo_count"]

    proportion_phased = (
        (reads_phased + 1) /
        (reads_phased + group["H0_noncyclo_count"] + 1)
    )

    stat = np.abs(
        np.log2((group["H1_noncyclo_count"] + 1) / (group["H2_noncyclo_count"] + 1)) *
        np.log2(reads_phased)
    )

    group["test_statistic"] = np.where(
        (proportion_phased < 0.1) | (reads_phased < 10),
        0,
        stat
    )

    return group

# def Cyclo_Allelic_Imbalance(group):
#     """
#     Calculate test statistic for Cyclo allelic imbalance.
#     This function calculates a unique test statistic for each sample in the group.
#     """
#     results = []
#     for _, row in group.iterrows():
#         # Calculate the test statistic for this specific sample

#         reads_phased = (row['H1_cyclo_count'] + row['H2_cyclo_count'])
#         proportion_phased = (reads_phased + 1) / (reads_phased + row['H0_cyclo_count'] + 1)

#         #If the proportion phased is too small, then the test stat can't be calculated.
#         if (proportion_phased < 0.1) | (reads_phased < 10):
#             test_statistic = 0
#         else:
#             test_statistic = abs(
#                 np.log2((row['H1_cyclo_count']+1) / (row['H2_cyclo_count']+1)) *
#                 np.log2(row['H1_cyclo_count'] + row['H2_cyclo_count'])
#             )

#         results.append(test_statistic)
    
#     # Assign the calculated test statistics back to the group
#     group['test_statistic'] = results
#     return group

def Cyclo_Allelic_Imbalance(group):
    """
    Calculate test statistic for Cyclo allelic imbalance.
    This function calculates a unique test statistic for each sample in the group.
    """
    group = group.copy()

    reads_phased = group["H1_cyclo_count"] + group["H2_cyclo_count"]

    proportion_phased = (
        (reads_phased + 1) /
        (reads_phased + group["H0_cyclo_count"] + 1)
    )

    stat = np.abs(
        np.log2((group["H1_cyclo_count"] + 1) / (group["H2_cyclo_count"] + 1)) *
        np.log2(reads_phased)
    )

    group["test_statistic"] = np.where(
        (proportion_phased < 0.1) | (reads_phased < 10),
        0,
        stat
    )

    return group

def Cyclo_Expression_Outlier_GOE_minor_isoforms(group):
    """Calculate test statistic for Cyclo Expression Outlier - Gain of Expression (GOE) only using minor isoforms for the gene."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        max_cyclo_TPM = other_samples['Minor_isoform_cyclo_reads'].max()
        median_cyclo_TPM = other_samples['Minor_isoform_cyclo_reads'].median()
        
        # Calculate the test statistic
        ratio = (row['Minor_isoform_cyclo_reads'] + 1) / (max_cyclo_TPM + 1)
        test_statistic = np.log2(ratio) * np.log2(median_cyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group

def Noncyclo_Expression_Outlier_GOE_minor_isoforms(group):
    """Calculate test statistic for Noncyclo Expression Outlier - Gain of Expression (GOE) only using minor isoforms for the gene."""
    results = []
    for _, row in group.iterrows():
        # Exclude the current sample and calculate metrics from other samples
        other_samples = group.drop(row.name)
        max_noncyclo_TPM = other_samples['Minor_isoform_noncyclo_reads'].max()
        median_noncyclo_TPM = other_samples['Minor_isoform_noncyclo_reads'].median()
        
        # Calculate the test statistic
        ratio = (row['Minor_isoform_noncyclo_reads'] + 1) / (max_noncyclo_TPM + 1)
        test_statistic = np.log2(ratio) * np.log2(median_noncyclo_TPM + 2)
        results.append(test_statistic)
    
    # Assign the calculated test statistics back to the group
    group['test_statistic'] = results
    return group


def process_hypothesis_test(filtered_data, group_col, test_statistic_func, gene_group_col=None, gene_level=True, bin_proportion=0.01, filter_before_ranking=True, filter_count_threshold=10):
    """
    Combine hypothesis testing, z-score calculation, ranking, and additional metrics into a single function.
    
    Parameters:
    - filtered_data (pd.DataFrame): The filtered data to process.
    - group_col (str): The column to group by (e.g., 'Isoform').
    - test_statistic_func (function): The hypothesis test function to apply.
    - gene_group_col (str, optional): The column to group by at the gene level. Defaults to group_col.
    - gene_level (bool): Whether to aggregate at the gene level before processing.
    - bin_proportion (float): Bin proportion threshold for low-abundance isoforms. Only used for NMD_rare_steady_state_transcript.
    
    Returns:
    - pd.DataFrame: Ranked data with additional calculated columns.
    """

    # Make a copy of the input data to avoid modifying the original
    filtered_data = filtered_data.copy()

    if gene_group_col is None:
        gene_group_col = group_col

    # Ensure the required columns are present
    required_columns = [
        "Isoform", "Sample", "cyclo_count", "noncyclo_count", "total_cyclo", "total_noncyclo",
        "Cyclo_TPM", "Noncyclo_TPM"
    ]
    missing_columns = [col for col in required_columns if col not in filtered_data.columns]
    if missing_columns:
        raise KeyError(f"The following required columns are missing from the input data: {missing_columns}")

    # Add gene-level metrics
    filtered_data["gene_cyclo_count"] = filtered_data.groupby([gene_group_col, "Sample"])["cyclo_count"].transform("sum")
    filtered_data["gene_noncyclo_count"] = filtered_data.groupby([gene_group_col, "Sample"])["noncyclo_count"].transform("sum")

    filtered_data["isoform_cyclo_proportion"] = filtered_data["cyclo_count"] / filtered_data["gene_cyclo_count"]
    filtered_data["isoform_noncyclo_proportion"] = filtered_data["noncyclo_count"] / filtered_data["gene_noncyclo_count"]

    # Gene-level aggregation if specified
    if gene_level:
        
        # Define all possible aggregation columns
        aggregation_dict = {
            "cyclo_count": ("cyclo_count", "sum"),
            "noncyclo_count": ("noncyclo_count", "sum"),
            "Cyclo_TPM": ("Cyclo_TPM", "sum"),
            "Noncyclo_TPM": ("Noncyclo_TPM", "sum"),
            "total_cyclo": ("total_cyclo", "first"),
            "total_noncyclo": ("total_noncyclo", "first"),
        }

        # Conditionally add H1/H2 columns if they exist
        optional_columns = [
            "H0_cyclo_count",
            "H1_cyclo_count",
            "H2_cyclo_count",
            "H0_noncyclo_count",
            "H1_noncyclo_count",
            "H2_noncyclo_count",
        ]
        for col in optional_columns:
            if col in filtered_data.columns:
                aggregation_dict[col] = (col, "sum")

        # Aggregate counts at the gene level using the dynamically constructed dictionary
        gene_level_data = (
            filtered_data.groupby([gene_group_col, "Sample"])
            .agg(**aggregation_dict)
            .reset_index()
        )

        # Recalculate Cyclo_TPM_rank and Noncyclo_TPM_rank
        # Calculate Cyclo_TPM_rank and Noncyclo_TPM_rank with average ranking for ties. Should go from 1 to number of patients. The lower the rank, the larger the TPM.
        gene_level_data["Cyclo_TPM_Rank"] = gene_level_data.groupby(gene_group_col)["Cyclo_TPM"].rank(ascending=False, method="average")
        gene_level_data["Noncyclo_TPM_Rank"] = gene_level_data.groupby(gene_group_col)["Noncyclo_TPM"].rank(ascending=False, method="average")
          
        # Recalculate mean and median TPM
        gene_level_data["Cyclo_TPM_Median"] = gene_level_data.groupby(gene_group_col)["Cyclo_TPM"].transform("median")
        gene_level_data["Noncyclo_TPM_Median"] = gene_level_data.groupby(gene_group_col)["Noncyclo_TPM"].transform("median")
        gene_level_data["Cyclo_TPM_Mean"] = gene_level_data.groupby(gene_group_col)["Cyclo_TPM"].transform("mean")
        gene_level_data["Noncyclo_TPM_Mean"] = gene_level_data.groupby(gene_group_col)["Noncyclo_TPM"].transform("mean")



        if (test_statistic_func == Cyclo_Expression_Outlier_GOE_minor_isoforms) | (test_statistic_func == Noncyclo_Expression_Outlier_GOE_minor_isoforms):
            # Create minor isoform counts for each gene
            
            # First, determine which isoforms are minor across all samples
            isoform_totals = (
                filtered_data.groupby(["Isoform", gene_group_col])["noncyclo_count"]
                .sum()
                .reset_index(name="isoform_noncyclo_total")
            )
            isoform_totals["gene_noncyclo_total"] = (
                isoform_totals.groupby(gene_group_col)["isoform_noncyclo_total"].transform("sum")
            )
            isoform_totals["isoform_fraction"] = (
                isoform_totals["isoform_noncyclo_total"] / isoform_totals["gene_noncyclo_total"]
            )
            
            isoform_totals["minor_isoform"] = (
                isoform_totals["isoform_fraction"] < 0.01
            ).fillna(False)

            # Merge this minor isoform status back to the sample-level data
            filtered_with_minor = filtered_data.merge(
                isoform_totals[["Isoform", gene_group_col, "minor_isoform"]],
                on=["Isoform", gene_group_col],
                how="left"
            )

            # After merge when some rows don’t get a minor_isoform value
            minor_mask = filtered_with_minor["minor_isoform"].fillna(False).astype(bool)


            # Aggregate using the mask. Total reads from minor isoforms for each (Sample, Gene)
            minor_reads = (
                filtered_with_minor.loc[minor_mask]
                .groupby([gene_group_col, "Sample"], as_index=False)
                .agg(
                    Minor_isoform_cyclo_reads=("cyclo_count", "sum"),
                    Minor_isoform_noncyclo_reads=("noncyclo_count", "sum"),
                )
            )

            # Merge into gene-level table
            gene_level_data = gene_level_data.merge(minor_reads, on=[gene_group_col, "Sample"], how="left")

            # Fill NaNs with 0 in case no minor isoforms exist for a gene in a sample
            gene_level_data[["Minor_isoform_cyclo_reads", "Minor_isoform_noncyclo_reads"]] = (
                gene_level_data[["Minor_isoform_cyclo_reads", "Minor_isoform_noncyclo_reads"]].fillna(0)
            )



        # if test_statistic_func == NMD_rare_steady_state_transcript:
        #     # Create bins and calculate aggregated values
        #     filtered_data["bin"] = filtered_data["isoform_noncyclo_proportion"].apply(
        #         lambda x: "Bin1_le" if x <= bin_proportion else "Bin2_g"
        #     )

        #     bin_aggregated = filtered_data.groupby(["Sample", gene_group_col, "bin"]).agg(
        #         Total_bin_cyclo_count=("cyclo_count", "sum"),
        #         Total_bin_noncyclo_count=("noncyclo_count", "sum")
        #     ).reset_index()

        #     # Pivot to wide format
        #     wide_result = bin_aggregated.pivot_table(
        #         index=["Sample", gene_group_col],
        #         columns="bin",
        #         values=["Total_bin_cyclo_count", "Total_bin_noncyclo_count"],
        #         fill_value=0
        #     )
        #     wide_result.columns = [
        #         f"{col[0]}_{col[1]}" for col in wide_result.columns.to_flat_index()
        #     ]
        #     wide_result.reset_index(inplace=True)

        #     # Merge with gene-level data
        #     gene_level_data = gene_level_data.merge(wide_result, on=["Sample", gene_group_col], how="left")

        #     # Calculate proportions and differences
        #     gene_level_data["proportion_in_Bin1_cyclo"] = gene_level_data["Total_bin_cyclo_count_Bin1_le"] / (
        #         gene_level_data["Total_bin_cyclo_count_Bin1_le"] + gene_level_data["Total_bin_cyclo_count_Bin2_g"]
        #     )
        #     gene_level_data["proportion_in_Bin1_noncyclo"] = gene_level_data["Total_bin_noncyclo_count_Bin1_le"] / (
        #         gene_level_data["Total_bin_noncyclo_count_Bin1_le"] + gene_level_data["Total_bin_noncyclo_count_Bin2_g"]
        #     )

        #     gene_level_data.fillna(0, inplace=True)
        #     gene_level_data["bin_proportion_difference"] = (
        #         gene_level_data["proportion_in_Bin1_cyclo"] - gene_level_data["proportion_in_Bin1_noncyclo"]
        #     ) / (
        #         gene_level_data["proportion_in_Bin1_cyclo"] + gene_level_data["proportion_in_Bin1_noncyclo"]
        #     )


        if test_statistic_func in [
            NMD_rare_steady_state_transcript,
            NMD_rare_steady_state_transcript_hap1,
            NMD_rare_steady_state_transcript_hap2,
        ]:
            # # Create bins and calculate aggregated values
            # filtered_data["bin"] = filtered_data["isoform_noncyclo_proportion"].apply(
            #     lambda x: "Bin1_le" if x <= bin_proportion else "Bin2_g"
            # )

            # Create bins and calculate aggregated values
            filtered_data["bin"] = np.where(
                filtered_data["isoform_noncyclo_proportion"] <= bin_proportion,
                "Bin1_le",
                "Bin2_g"
            )

            # Choose which columns to aggregate (total vs hap1 vs hap2)
            if test_statistic_func == NMD_rare_steady_state_transcript:
                cyclo_col = "cyclo_count"
                noncyclo_col = "noncyclo_count"
            elif test_statistic_func == NMD_rare_steady_state_transcript_hap1:
                cyclo_col = "H1_cyclo_count"
                noncyclo_col = "H1_noncyclo_count"
            else:  # NMD_rare_steady_state_transcript_hap2
                cyclo_col = "H2_cyclo_count"
                noncyclo_col = "H2_noncyclo_count"

            bin_aggregated = filtered_data.groupby(["Sample", gene_group_col, "bin"]).agg(
                Total_bin_cyclo_count=(cyclo_col, "sum"),
                Total_bin_noncyclo_count=(noncyclo_col, "sum")
            ).reset_index()

            # Pivot to wide format
            wide_result = bin_aggregated.pivot_table(
                index=["Sample", gene_group_col],
                columns="bin",
                values=["Total_bin_cyclo_count", "Total_bin_noncyclo_count"],
                fill_value=0
            )
            wide_result.columns = [
                f"{col[0]}_{col[1]}" for col in wide_result.columns.to_flat_index()
            ]
            wide_result.reset_index(inplace=True)

            # Merge with gene-level data
            gene_level_data = gene_level_data.merge(wide_result, on=["Sample", gene_group_col], how="left")

            # Calculate proportions and differences
            gene_level_data["proportion_in_Bin1_cyclo"] = gene_level_data["Total_bin_cyclo_count_Bin1_le"] / (
                gene_level_data["Total_bin_cyclo_count_Bin1_le"] + gene_level_data["Total_bin_cyclo_count_Bin2_g"]
            )
            gene_level_data["proportion_in_Bin1_noncyclo"] = gene_level_data["Total_bin_noncyclo_count_Bin1_le"] / (
                gene_level_data["Total_bin_noncyclo_count_Bin1_le"] + gene_level_data["Total_bin_noncyclo_count_Bin2_g"]
            )

            gene_level_data.fillna(0, inplace=True)
            gene_level_data["bin_proportion_difference"] = (
                gene_level_data["proportion_in_Bin1_cyclo"] - gene_level_data["proportion_in_Bin1_noncyclo"]
            ) / (
                gene_level_data["proportion_in_Bin1_cyclo"] + gene_level_data["proportion_in_Bin1_noncyclo"]
            )

        processed_data = gene_level_data

        group_col = gene_group_col

    else:
        processed_data = filtered_data


    # Add additional metrics based on the test statistic function
    if test_statistic_func == NMD_test_statistic:

        processed_data["CycloFraction"] = processed_data["cyclo_count"] / processed_data["total_cyclo"]
        processed_data["NoncycloFraction"] = processed_data["noncyclo_count"] / processed_data["total_noncyclo"]
        processed_data["NormalizedCycloFraction"] = processed_data["CycloFraction"] / (
            processed_data["CycloFraction"] + processed_data["NoncycloFraction"]
        )
        processed_data["NormalizedNoncycloFraction"] = processed_data["NoncycloFraction"] / (
            processed_data["CycloFraction"] + processed_data["NoncycloFraction"]
        )
        processed_data["NormalizedFractionDifference"] = (
            processed_data["NormalizedCycloFraction"] - processed_data["NormalizedNoncycloFraction"]
        )

    elif test_statistic_func == NMD_hap1_test_statistic:

        processed_data["CycloFraction_H1"] = processed_data["H1_cyclo_count"] / processed_data["total_cyclo"]
        processed_data["NoncycloFraction_H1"] = processed_data["H1_noncyclo_count"] / processed_data["total_noncyclo"]
        processed_data["NormalizedCycloFraction_H1"] = processed_data["CycloFraction_H1"] / (
            processed_data["CycloFraction_H1"] + processed_data["NoncycloFraction_H1"]
        )
        processed_data["NormalizedNoncycloFraction_H1"] = processed_data["NoncycloFraction_H1"] / (
            processed_data["CycloFraction_H1"] + processed_data["NoncycloFraction_H1"]
        )
        processed_data["NormalizedFractionDifference_H1"] = (
            processed_data["NormalizedCycloFraction_H1"] - processed_data["NormalizedNoncycloFraction_H1"]
        )

    elif test_statistic_func == NMD_hap2_test_statistic:

        processed_data["CycloFraction_H2"] = processed_data["H2_cyclo_count"] / processed_data["total_cyclo"]
        processed_data["NoncycloFraction_H2"] = processed_data["H2_noncyclo_count"] / processed_data["total_noncyclo"]
        processed_data["NormalizedCycloFraction_H2"] = processed_data["CycloFraction_H2"] / (
            processed_data["CycloFraction_H2"] + processed_data["NoncycloFraction_H2"]
        )
        processed_data["NormalizedNoncycloFraction_H2"] = processed_data["NoncycloFraction_H2"] / (
            processed_data["CycloFraction_H2"] + processed_data["NoncycloFraction_H2"]
        )
        processed_data["NormalizedFractionDifference_H2"] = (
            processed_data["NormalizedCycloFraction_H2"] - processed_data["NormalizedNoncycloFraction_H2"]
        )

    elif test_statistic_func in [Noncyclo_Expression_Outlier_LOE, Noncyclo_Expression_Outlier_GOE]:
        processed_data["Avg_Noncyclo_TPM"] = processed_data.groupby(group_col)["Noncyclo_TPM"].transform("mean")
        processed_data["SD_Noncyclo_TPM"] = processed_data.groupby(group_col)["Noncyclo_TPM"].transform("std")
        processed_data["Noncyclo_TPM_Z_Score"] = (
            processed_data["Noncyclo_TPM"] - processed_data["Avg_Noncyclo_TPM"]
        ) / processed_data["SD_Noncyclo_TPM"]

    elif test_statistic_func in [Cyclo_Expression_Outlier_LOE, Cyclo_Expression_Outlier_GOE]:
        processed_data["Avg_Cyclo_TPM"] = processed_data.groupby(group_col)["Cyclo_TPM"].transform("mean")
        processed_data["SD_Cyclo_TPM"] = processed_data.groupby(group_col)["Cyclo_TPM"].transform("std")
        processed_data["Cyclo_TPM_Z_Score"] = (
            processed_data["Cyclo_TPM"] - processed_data["Avg_Cyclo_TPM"]
        ) / processed_data["SD_Cyclo_TPM"]
    elif test_statistic_func == Noncyclo_Allelic_Imbalance:
        processed_data["NoncycloHaplotypeDifference"] = (
            (processed_data["H1_noncyclo_count"] - processed_data["H2_noncyclo_count"]) / 
            (processed_data["H1_noncyclo_count"] + processed_data["H2_noncyclo_count"])
        )
    elif test_statistic_func == Cyclo_Allelic_Imbalance:
        processed_data["CycloHaplotypeDifference"] = (
            (processed_data["H1_cyclo_count"] - processed_data["H2_cyclo_count"]) / 
            (processed_data["H1_cyclo_count"] + processed_data["H2_cyclo_count"])
        )
    elif test_statistic_func in [Cyclo_Expression_Outlier_GOE_minor_isoforms]:
        processed_data["Avg_Minor_isoform_cyclo_reads"] = processed_data.groupby(group_col)["Minor_isoform_cyclo_reads"].transform("mean")
        processed_data["SD_Minor_isoform_cyclo_reads"] = processed_data.groupby(group_col)["Minor_isoform_cyclo_reads"].transform("std")
        processed_data["Minor_isoform_cyclo_reads_Z_Score"] = (
            processed_data["Minor_isoform_cyclo_reads"] - processed_data["Avg_Minor_isoform_cyclo_reads"]
        ) / processed_data["SD_Minor_isoform_cyclo_reads"]
    elif test_statistic_func in [Noncyclo_Expression_Outlier_GOE_minor_isoforms]:
        processed_data["Avg_Minor_isoform_noncyclo_reads"] = processed_data.groupby(group_col)["Minor_isoform_noncyclo_reads"].transform("mean")
        processed_data["SD_Minor_isoform_noncyclo_reads"] = processed_data.groupby(group_col)["Minor_isoform_noncyclo_reads"].transform("std")
        processed_data["Minor_isoform_noncyclo_reads_Z_Score"] = (
            processed_data["Minor_isoform_noncyclo_reads"] - processed_data["Avg_Minor_isoform_noncyclo_reads"]
        ) / processed_data["SD_Minor_isoform_noncyclo_reads"]


    if gene_level == True:
        processed_data = filter_based_on_counts(processed_data, count_threshold=filter_count_threshold, group_col=gene_group_col)
    else:
        processed_data = filter_based_on_counts(processed_data, count_threshold=filter_count_threshold, group_col=group_col)

    # Apply hypothesis test
    tested_data = apply_hypothesis_test(processed_data, group_col=gene_group_col if gene_level else group_col, test_statistic_func=test_statistic_func)
    
    # Calculate z-scores
    z_scored_data = calculate_z_score_MAD(tested_data, group_col=gene_group_col if gene_level else group_col, stat_col="test_statistic")

    # Output the z_scored_data as an intermediate output. This includes the test stat and z-scores for all genes/isoforms, even unranked ones.
    test_name = test_statistic_func.__name__
    level = "gen" if gene_level else "iso"
    output_file = f"{test_name}_{level}_full_test_stat_and_z_scored_data.tsv.gz"
    z_scored_data.to_csv(output_file, index=False, compression="gzip", sep="\t")

    #Filter before ranking
    if filter_before_ranking == True:
        
        # Filter out minor isoforms (except when doing LOE analysis)
        #if not gene_level and test_statistic_func not in {Noncyclo_Expression_Outlier_LOE, Cyclo_Expression_Outlier_LOE}:
        #    z_scored_data = z_scored_data[(z_scored_data["isoform_cyclo_proportion"] > 0.1) | 
        #                                  (z_scored_data["isoform_noncyclo_proportion"] > 0.1)]


        if test_statistic_func in [
            NMD_rare_steady_state_transcript,
            NMD_rare_steady_state_transcript_hap1,
            NMD_rare_steady_state_transcript_hap2,
        ]:
            z_scored_data = z_scored_data[z_scored_data["bin_proportion_difference"] > 0]
            z_scored_data = z_scored_data[z_scored_data["Total_bin_cyclo_count_Bin1_le"] > 10]
        elif test_statistic_func == NMD_test_statistic:
            z_scored_data = z_scored_data[z_scored_data["NormalizedFractionDifference"] > 0]
        elif test_statistic_func == NMD_hap1_test_statistic:
            z_scored_data = z_scored_data[z_scored_data["NormalizedFractionDifference_H1"] > 0]
        elif test_statistic_func == NMD_hap2_test_statistic:
            z_scored_data = z_scored_data[z_scored_data["NormalizedFractionDifference_H2"] > 0]
        elif test_statistic_func == Noncyclo_Expression_Outlier_LOE:
            z_scored_data = z_scored_data[z_scored_data["Noncyclo_TPM_Z_Score"] < 0]
        elif test_statistic_func == Noncyclo_Expression_Outlier_GOE:
            z_scored_data = z_scored_data[z_scored_data["Noncyclo_TPM_Z_Score"] > 0]
        elif test_statistic_func == Cyclo_Expression_Outlier_LOE:
            z_scored_data = z_scored_data[z_scored_data["Cyclo_TPM_Z_Score"] < 0]
        elif test_statistic_func == Cyclo_Expression_Outlier_GOE:
            z_scored_data = z_scored_data[z_scored_data["Cyclo_TPM_Z_Score"] > 0]
        elif test_statistic_func == Cyclo_Expression_Outlier_GOE_minor_isoforms:
            z_scored_data = z_scored_data[z_scored_data["Minor_isoform_cyclo_reads_Z_Score"] > 0]
        elif test_statistic_func == Noncyclo_Expression_Outlier_GOE_minor_isoforms:
            z_scored_data = z_scored_data[z_scored_data["Minor_isoform_noncyclo_reads_Z_Score"] > 0]
            


    # Calculate ranks
    ranked_data = calculate_ranks_for_sample(z_scored_data, group_col=gene_group_col if gene_level else group_col)
    
    return ranked_data


# Wrapping because "process_test" might be a better name
def process_test(*args, **kwargs):
    return process_hypothesis_test(*args, **kwargs)