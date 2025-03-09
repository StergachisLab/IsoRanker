#!/usr/bin/env python3

import pandas as pd
import argparse
import os
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

    # Store full results to generate lookup table
    full_ranked_gene_data = []

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

        # Append tuple (test_name, ranked_data) to the list
        full_ranked_gene_data.append((test_name, ranked_data))

        #rank_columns = [col for col in ranked_data.columns if col.startswith('rank_top_')]
        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        filtered_ranked_data = filtered_ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],
            how='left',
            left_on='associated_gene',
            right_on='Approved Gene Symbol'
        ).drop(columns=['Approved Gene Symbol'])

        output_file = os.path.join(output_dir, f"{test_name}_gene_data.csv")
        filtered_ranked_data.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}", flush=True)

    # Run isoform-level analysis
    for test_name, test_func in test_stat_funcs[:-1]:  # Exclude the last test for isoform-level analysis
        print(f"Processing test statistic: {test_name}", flush=True)
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated,
            group_col='Isoform',
            test_statistic_func=test_func,
            gene_group_col='associated_gene',
            gene_level=False,
            bin_proportion=0.01
        )

        #rank_columns = [col for col in ranked_data.columns if col.startswith('rank_top_')]
        #filtered_ranked_data = ranked_data[ranked_data[rank_columns].le(1000).any(axis=1)]
        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        filtered_ranked_data = ranked_data

        filtered_ranked_data = filtered_ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],
            how='left',
            left_on='associated_gene',
            right_on='Approved Gene Symbol'
        ).drop(columns=['Approved Gene Symbol'])

        output_file = os.path.join(output_dir, f"{test_name}_isoform_data.csv")
        filtered_ranked_data.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}", flush=True)


    ################################
    # Create lookup tables
    ################################

    # Collapse long_format_annotated by summing TPM values per Sample and Gene
    sample_gene_rankings_lookup_table = long_format_annotated.groupby(["Sample", "associated_gene"], as_index=False).agg(
        {"Cyclo_TPM": "sum", "Noncyclo_TPM": "sum"}
    )

    # Now merge the ranked gene-level data
    for test_name, df in full_ranked_gene_data:
        # Rename rank column to be test-specific
        df_renamed = df.rename(columns={"rank_top_99_5_percentile": f"{test_name}_rank_top_99_5_percentile"})

        # Keep only relevant columns
        df_renamed = df_renamed[["Sample", "associated_gene", f"{test_name}_rank_top_99_5_percentile"]]

        # Merge into merged_df using outer join
        sample_gene_rankings_lookup_table = pd.merge(sample_gene_rankings_lookup_table, df_renamed, on=["Sample", "associated_gene"], how="outer")

    # Save as a compressed CSV (gzip format)
    sample_gene_rankings_lookup_table.to_csv("sample_gene_rankings_lookup_table.csv.gz", index=False, compression="gzip")

    # Group by gene (associated_gene) and compute median, Q1 (25th percentile), and Q3 (75th percentile)
    gene_coverage_lookup_table = sample_gene_rankings_lookup_table.groupby("associated_gene").agg(
        Cyclo_TPM_median=("Cyclo_TPM", "median"),
        Cyclo_TPM_Q1=("Cyclo_TPM", lambda x: x.quantile(0.25)),  # 25th percentile
        Cyclo_TPM_Q3=("Cyclo_TPM", lambda x: x.quantile(0.75)),  # 75th percentile
        Cyclo_TPM_min=("Cyclo_TPM", "min"),  # Minimum value
        Cyclo_TPM_max=("Cyclo_TPM", "max"),  # Maximum value

        Noncyclo_TPM_median=("Noncyclo_TPM", "median"),
        Noncyclo_TPM_Q1=("Noncyclo_TPM", lambda x: x.quantile(0.25)),  # 25th percentile
        Noncyclo_TPM_Q3=("Noncyclo_TPM", lambda x: x.quantile(0.75)),  # 75th percentile
        Noncyclo_TPM_min=("Noncyclo_TPM", "min"),  # Minimum value
        Noncyclo_TPM_max=("Noncyclo_TPM", "max")   # Maximum value
    ).reset_index()

    gene_coverage_lookup_table.to_csv("gene_coverage_lookup_table.csv.gz", index=False, compression="gzip")


if __name__ == "__main__":
    main()