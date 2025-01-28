#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import pyreadr  # To save as RDS
from IsoRanker import (
    create_expression_matrix,
    create_long_format,
    filter_based_on_counts,
    process_hypothesis_test,
    NMD_test_statistic,
    Noncyclo_Expression_Outlier_LOE,
    Noncyclo_Expression_Outlier_GOE,
    Cyclo_Expression_Outlier_GOE,
    NMD_rare_steady_state_transcript,
)

def save_as_rds(df, file_path):
    """Save DataFrame as an RDS file."""
    try:
        pyreadr.write_rds(file_path, df, compress="gzip")
        print(f"RDS file saved to {file_path}", flush=True)
    except Exception as e:
        print(f"Error saving RDS file: {e}", flush=True)

def select_columns(df, test_stat_func, gene_level):
    """Select specific columns for saving as RDS based on the test_statistic_func and analysis level."""
    common_columns = [
        "associated_gene", "Sample", "Cyclo_TPM", "Noncyclo_TPM", "test_statistic",
        "z_score_of_test_stat", "rank_top_99_9_percentile", "rank_top_99_5_percentile",
        "rank_top_99_percentile", "rank_top_98_percentile", "rank_top_95_percentile", "Phenotypes"
    ]

    if test_stat_func == NMD_test_statistic:
        columns = common_columns + ["NormalizedFractionDifference"]
    elif test_stat_func in [Noncyclo_Expression_Outlier_LOE, Noncyclo_Expression_Outlier_GOE]:
        columns = common_columns + ["Cyclo_TPM_Rank", "Noncyclo_TPM_Rank", "Noncyclo_TPM_Z_Score"]
    elif test_stat_func in [Cyclo_Expression_Outlier_GOE, Cyclo_Expression_Outlier_LOE]:
        columns = common_columns + ["Cyclo_TPM_Rank", "Noncyclo_TPM_Rank", "Cyclo_TPM_Z_Score"]
    elif test_stat_func == NMD_rare_steady_state_transcript:
        columns = common_columns + ["proportion_in_Bin1_cyclo", "proportion_in_Bin1_noncyclo", "bin_proportion_difference"]
    else:
        columns = common_columns  # Default to common columns

    if not gene_level:
        # Include isoform-specific columns for isoform-level analysis
        columns += ["Isoform", "isoform_cyclo_proportion", "isoform_noncyclo_proportion"]

    # Select only the columns that exist in the DataFrame
    return df[[col for col in columns if col in df.columns]]

def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run analysis pipeline.")
    parser.add_argument("--read_stat_path", required=True, help="Path to the read_stat.txt file.")
    parser.add_argument("--sample_info_path", required=True, help="Path to the sample info CSV file.")
    parser.add_argument("--classification_path", required=True, help="Path to the pigeon classification file.")
    parser.add_argument("--genemap_path", required=True, help="Path to the genemap2.txt file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the output files.")

    args = parser.parse_args()

    # Assign command-line arguments to variables
    read_stat_path = args.read_stat_path
    sample_info_path = args.sample_info_path
    classification_path = args.classification_path
    genemap_path = args.genemap_path
    output_dir = args.output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load input files
    print("Reading input files", flush=True)
    sample_info = pd.read_csv(sample_info_path)
    classification_data = pd.read_csv(classification_path, sep="\t")
    genemap = pd.read_csv(genemap_path, sep='\t', skiprows=3) # Read the file, skipping the first 3 rows
    genemap = genemap[genemap['Approved Gene Symbol'].notnull()]

    # Create the expression matrix
    print("Creating expression matrix", flush=True)
    expression_matrix = create_expression_matrix(read_stat_path, output_file=os.path.join(output_dir, "expression_matrix.csv"))

    # Generate long-format DataFrame
    print("Creating long format expression matrix", flush=True)
    long_format_df = create_long_format(expression_matrix, sample_info)

    # Filter input
    long_format_df = filter_based_on_counts(long_format_df, count_threshold=10, group_col='Isoform')

    # Merge with classification data
    classification_subset = classification_data[['isoform', 'associated_gene']]
    long_format_annotated = long_format_df.merge(
        classification_subset,
        left_on="Isoform",
        right_on="isoform",
        how="left"
    ).drop(columns=["isoform"])

    long_format_annotated.to_csv(os.path.join(output_dir, "long_format_annotated.csv"), index=False)

    # Define hypothesis tests
    test_stat_funcs = [
        ("NMD", NMD_test_statistic),
        ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
        ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
        ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE),
        ("NMD_rare_steady_state_transcript", NMD_rare_steady_state_transcript)
    ]

    # Run gene-level analysis
    for test_name, test_func in test_stat_funcs:
        print(f"Processing test statistic: {test_name}", flush=True)
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated,
            group_col='Isoform',
            test_statistic_func=test_func,
            gene_group_col='associated_gene',
            gene_level=True,
            bin_proportion=0.01
        )

        filtered_ranked_data = ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],
            how='left',
            left_on='associated_gene',
            right_on='Approved Gene Symbol'
        ).drop(columns=['Approved Gene Symbol'])

        # Select specific columns and save as RDS
        selected_data = select_columns(filtered_ranked_data, test_func, gene_level=True)
        output_rds = os.path.join(output_dir, f"{test_name}_gene_data.rds")
        save_as_rds(selected_data, output_rds)


    # Define hypothesis tests
    test_stat_funcs = [
        ("NMD", NMD_test_statistic),
        ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
        ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
        ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE)
    ]

    # Run isoform-level analysis
    for test_name, test_func in test_stat_funcs:  # Exclude the last test for isoform-level analysis
        print(f"Processing test statistic: {test_name}", flush=True)
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated,
            group_col='Isoform',
            test_statistic_func=test_func,
            gene_group_col='associated_gene',
            gene_level=False,
            bin_proportion=0.01
        )

        filtered_ranked_data = ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],
            how='left',
            left_on='associated_gene',
            right_on='Approved Gene Symbol'
        ).drop(columns=['Approved Gene Symbol'])

        # Select specific columns and save as RDS
        selected_data = select_columns(filtered_ranked_data, test_func, gene_level=False)
        output_rds = os.path.join(output_dir, f"{test_name}_isoform_data.rds")
        save_as_rds(selected_data, output_rds)


if __name__ == "__main__":
    main()
