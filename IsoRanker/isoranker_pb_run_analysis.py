#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import shutil
import argparse

from IsoRanker import (
    load_data,
    filter_based_on_counts,
    apply_hypothesis_test,
    calculate_z_score,
    NMD_test_statistic,
    Noncyclo_Expression_Outlier_LOE,
    Noncyclo_Expression_Outlier_GOE,
    Cyclo_Expression_Outlier_GOE,
    NMD_rare_steady_state_transcript,
    Noncyclo_Allelic_Imbalance,
    Cyclo_Allelic_Imbalance,
    calculate_ranks_for_sample,
    create_expression_matrix,
    create_long_format,
    process_hypothesis_test,
    update_files_with_haplotype_info,
    merge_csvs_by_keyword, 
    process_vep_vcf, 
    merge_haplotype_data, 
    process_phenotype_data,
    process_and_plot_pca,
    analyze_isoforms,
    process_pileup
)

def main():

    ################################################
    # Parse command-line arguments
    ################################################
    parser = argparse.ArgumentParser(description="Run analysis pipeline.")
    parser.add_argument("--read_stat_path", required=True, help="Path to the read_stat.txt file.")
    parser.add_argument("--sample_info_path", required=True, help="Path to the sample info CSV file.")
    parser.add_argument("--classification_path", required=True, help="Path to the pigeon classification file.")
    parser.add_argument("--genemap_path", required=True, help="Path to the genemap2.txt file.")
    parser.add_argument("--hpo_file_path", required=True, help="Path to the file that matches HPO terms to OMIM. Download from here: https://hpo.jax.org/data/annotations.")
    parser.add_argument("--probands_file_path", required=True, help="Path to the file that contains proband HPO terms.")
    parser.add_argument("--reference_fasta_path", required=True, help="Path to the hg38 reference fasta path.")
    #parser.add_argument("--output_dir", required=True, help="Directory to save the output files.")

    args = parser.parse_args()

    ################################################
    # Assign command-line arguments to variables
    ################################################
    read_stat_path = args.read_stat_path
    sample_info_path = args.sample_info_path
    classification_path = args.classification_path
    genemap_path = args.genemap_path
    hpo_file_path = args.hpo_file_path
    probands_file_path = args.probands_file_path
    reference_fasta_path = args.reference_fasta_path
    output_dir = "."
    #output_dir = args.output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    ################################################
    # Load input files
    ################################################
    print("Reading input files", flush=True)
    sample_info = pd.read_csv(sample_info_path)
    classification_data = pd.read_csv(classification_path, sep="\t")
    genemap = pd.read_csv(genemap_path, sep='\t', skiprows=3) # Read the file, skipping the first 3 rows
    genemap = genemap[genemap['Approved Gene Symbol'].notnull()]


    ################################################
    # Update input files with haplotype information
    ################################################
    print("Updating input files with haplotype information", flush=True)
    update_files_with_haplotype_info(sample_info_path, read_stat_path, output_dir)
    read_stat_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
    sample_info_path = os.path.join(output_dir, "updated_sample_info.csv.gz")
    sample_info = pd.read_csv(sample_info_path)

    ################################################
    # Create the expression matrix
    ################################################
    print("Creating expression matrix", flush=True)
    expression_matrix = create_expression_matrix(read_stat_path, output_file=os.path.join(output_dir, "expression_matrix.csv.gz"))

    ################################################
    # Generate long-format DataFrame
    ################################################
    print("Creating long format expression matrix", flush=True)
    long_format_df = create_long_format(expression_matrix, sample_info)

    # Merge with classification data
    classification_subset = classification_data[['isoform', 'associated_gene']]
    long_format_annotated = long_format_df.merge(
        classification_subset,
        left_on="Isoform",
        right_on="isoform",
        how="left"
    ).drop(columns=["isoform"])

    long_format_annotated.to_csv(os.path.join(output_dir, "long_format_annotated.csv.gz"), index=False, compression = "gzip")

    ################################################
    # Gene level hypothesis testing
    ################################################

    # Define hypothesis tests
    test_stat_funcs = [
        ("NMD", NMD_test_statistic),
        ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
        ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
        ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE),
        ("NMD_rare_steady_state_transcript", NMD_rare_steady_state_transcript),
        ("Nonyclo_Allelic_Imbalance", Noncyclo_Allelic_Imbalance),
        ("Cyclo_Allelic_Imbalance", Cyclo_Allelic_Imbalance)
    ]

    # Store full results to generate lookup table
    full_ranked_gene_data = []

    #Gene level
    for test_name, test_func in test_stat_funcs:
        print(f"Processing test statistic: {test_name}", flush=True)

        # Apply the process_hypothesis_test function
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated, 
            group_col='Isoform', 
            test_statistic_func=test_func, 
            gene_group_col='associated_gene', 
            gene_level=True, 
            bin_proportion=0.01, 
            filter_before_ranking=True, 
            filter_count_threshold=10)

        # Append tuple (test_name, ranked_data) to the list
        full_ranked_gene_data.append((test_name, ranked_data))

        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        # Add OMIM data to genes
        filtered_ranked_data = filtered_ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],  # Select relevant columns from genemap
            how='left',  # Perform a left join to keep all rows from filtered_ranked_data
            left_on='associated_gene',  # Column in filtered_ranked_data to join on
            right_on='Approved Gene Symbol'  # Column in genemap to join on
        )
        # Drop the 'Approved Gene Name' column if it is no longer needed
        filtered_ranked_data = filtered_ranked_data.drop(columns=['Approved Gene Symbol'])

        # Save the results to a CSV file
        output_file = os.path.join(output_dir, f"{test_name}_gene_top_ranked_data.csv.gz")
        filtered_ranked_data.to_csv(output_file, index=False, compression = "gzip")
        print(f"Results saved to {output_file}", flush=True)

    ################################################
    # Isoform level hypothesis testing
    ################################################

    test_stat_funcs = [
        ("NMD", NMD_test_statistic),
        ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
        ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
        ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE)
    ]

    #Isoform level
    for test_name, test_func in test_stat_funcs:
        print(f"Processing test statistic: {test_name}", flush=True)

        # Apply the process_hypothesis_test function
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated, 
            group_col='Isoform', 
            test_statistic_func=test_func, 
            gene_group_col='associated_gene', 
            gene_level=False, 
            bin_proportion=0.01, 
            filter_before_ranking=True, 
            filter_count_threshold=10)

        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        # Add OMIM data to genes
        filtered_ranked_data = filtered_ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],  # Select relevant columns from genemap
            how='left',  # Perform a left join to keep all rows from filtered_ranked_data
            left_on='associated_gene',  # Column in filtered_ranked_data to join on
            right_on='Approved Gene Symbol'  # Column in genemap to join on
        )
        # Drop the 'Approved Gene Name' column if it is no longer needed
        filtered_ranked_data = filtered_ranked_data.drop(columns=['Approved Gene Symbol'])

        # Save the results to a CSV file
        output_file = os.path.join(output_dir, f"{test_name}_isoform_top_ranked_data.csv.gz")
        filtered_ranked_data.to_csv(output_file, index=False, compression = "gzip")
        print(f"Results saved to {output_file}", flush=True)


    ################################################
    # Combine output files
    ################################################
    print("Combining output files", flush=True)

    # Isoform
    keyword = "isoform" 
    output_csv = os.path.join(output_dir, f"merged_ranked_{keyword}.csv.gz")
    merge_csvs_by_keyword(output_dir, keyword, output_csv)

    # Gene
    keyword = "gene"
    output_csv = os.path.join(output_dir, f"merged_ranked_{keyword}.csv.gz")
    merge_csvs_by_keyword(output_dir, keyword, output_csv)


    ################################################
    # Add patient phenotype information
    ################################################
    print("Adding patient phenotype information to combined output files", flush=True)

    all_comparisons, all_comparisons_long = process_phenotype_data(hpo_file_path, genemap_path, probands_file_path)

    # Master file gene
    master_file = pd.read_csv("merged_ranked_gene.csv.gz", compression = "gzip")

    # Merge `master_file` with `all_comparisons_long` based on `Sample` and gene name
    merged_data = master_file.merge(
        all_comparisons_long,
        left_on=['Sample', 'associated_gene'],  # Match Sample and associated gene
        right_on=['Sample', 'Approved Gene Symbol'],  # Match with Approved Gene Symbol
        how='left'  # Keep all rows from master_file
    )

    # Save merged output
    merged_data.to_csv("merged_ranked_gene_with_phenotype.csv.gz", index=False, compression = "gzip")


    # Master file isoform
    master_file = pd.read_csv("merged_ranked_isoform.csv.gz", compression = "gzip")

    # Merge `master_file` with `all_comparisons_long` based on `Sample` and gene name
    merged_data = master_file.merge(
        all_comparisons_long,
        left_on=['Sample', 'associated_gene'],  # Match Sample and associated gene
        right_on=['Sample', 'Approved Gene Symbol'],  # Match with Approved Gene Symbol
        how='left'  # Keep all rows from master_file
    )

    # Save merged output
    merged_data.to_csv("merged_ranked_isoform_with_phenotype.csv.gz", index=False, compression = "gzip")

    ################################################
    # Create lookup tables
    ################################################
    print("Creating lookup tables", flush=True)

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


    ################################################
    # QC
    ################################################

    ################################################
    # PCA
    ################################################

    print("Creating PCA plot", flush=True)

    pca_results = process_and_plot_pca(long_format_annotated, output_pdf="pca_plot.pdf", grouping_col = "associated_gene")

    ################################################
    # Gene diversity
    ################################################

    print("Analyzing gene diversity", flush=True)

    analyze_isoforms(long_format_annotated, "gene_diversity.csv.gz", "associated_gene")

    # Plotting:

    df = pd.read_csv("gene_diversity.csv.gz", compression = "gzip")

    # Drop 'Cyclo Total Reads' and 'Noncyclo Total Reads' columns
    df_filtered = df.drop(columns=["Cyclo Total Reads", "Noncyclo Total Reads"])

    # Define colors for Cyclo and Noncyclo categories
    cyclo_color = "red"
    noncyclo_color = "blue"

    # Determine a common y-axis limit for all plots
    y_max = df_filtered.drop(columns=["Sample"]).max().max()

    # Create subplots for each sample
    num_samples = len(df_filtered["Sample"])
    fig, axes = plt.subplots(nrows=num_samples, figsize=(12, num_samples * 3), sharex=True, sharey=True)

    # If only one sample, make axes iterable
    if num_samples == 1:
        axes = [axes]

    # Loop through each sample and create a separate bar plot
    for ax, sample in zip(axes, df_filtered["Sample"]):
        sample_data = df_filtered[df_filtered["Sample"] == sample].drop(columns=["Sample"]).T
        colors = [cyclo_color if "Cyclo" in col else noncyclo_color for col in sample_data.index]
        
        ax.bar(sample_data.index, sample_data.iloc[:, 0], color=colors)
        ax.set_title(f"Sample: {sample}")
        ax.set_ylabel("Count")
        ax.set_ylim(0, y_max)  # Set common y-axis limit
        ax.tick_params(axis="x", rotation=90)

    # Formatting
    plt.tight_layout()

    # Save to PDF
    plt.savefig("gene_diversity.pdf", format="pdf")


    ################################################
    # Isoform diversity
    ################################################

    print("Analzying isoform diversity", flush=True)

    analyze_isoforms(long_format_annotated, "isoform_diversity.csv.gz", "Isoform")

    df = pd.read_csv("isoform_diversity.csv.gz", compression = "gzip")

    # Drop 'Cyclo Total Reads' and 'Noncyclo Total Reads' columns
    df_filtered = df.drop(columns=["Cyclo Total Reads", "Noncyclo Total Reads"])

    # Define colors for Cyclo and Noncyclo categories
    cyclo_color = "red"
    noncyclo_color = "blue"

    # Determine a common y-axis limit for all plots
    y_max = df_filtered.drop(columns=["Sample"]).max().max()

    # Create subplots for each sample
    num_samples = len(df_filtered["Sample"])
    fig, axes = plt.subplots(nrows=num_samples, figsize=(12, num_samples * 3), sharex=True, sharey=True)

    # If only one sample, make axes iterable
    if num_samples == 1:
        axes = [axes]

    # Loop through each sample and create a separate bar plot
    for ax, sample in zip(axes, df_filtered["Sample"]):
        sample_data = df_filtered[df_filtered["Sample"] == sample].drop(columns=["Sample"]).T
        colors = [cyclo_color if "Cyclo" in col else noncyclo_color for col in sample_data.index]
        
        ax.bar(sample_data.index, sample_data.iloc[:, 0], color=colors)
        ax.set_title(f"Sample: {sample}")
        ax.set_ylabel("Count")
        ax.set_ylim(0, y_max)  # Set common y-axis limit
        ax.tick_params(axis="x", rotation=90)

    # Formatting
    plt.tight_layout()

    # Save to PDF
    plt.savefig("isoform_diversity.pdf", format="pdf")

    ###################################
    # SRSF6
    ###################################

    print("Analzying SRSF6 cassette exon inclusion", flush=True)

    process_pileup(df=sample_info, reference_fasta= reference_fasta_path, chromosome="chr20", position=43459200, output_file="SRSF6.csv.gz")

    # Plotting:

    SRSF6_df = pd.read_csv("SRSF6.csv.gz", compression = "gzip")

    # Sort DataFrame by source
    df_sorted = SRSF6_df.sort_values(by="Source")

    # Plot bar chart
    plt.figure(figsize=(15, 9))
    sns.barplot(data=df_sorted, x="Source", y="Exonic_Proportion", edgecolor="black")

    # Labeling and formatting
    plt.xlabel("Source")
    plt.ylabel("Exonic Proportion")
    plt.title("Bar Plot of Exonic Proportion by Source")
    plt.xticks(rotation=45, fontsize=7, ha="right")  # Adjust alignment to prevent cutoff

    # Save to PDF
    plt.savefig("SRSF6_exonic_proportion.pdf", format="pdf")


    ###################################
    # Organize output files
    ###################################

    # Define folder names
    QC_FOLDER = "qc"
    BROWSER_FOLDER = "browser"
    LOOKUP_TABLES_FOLDER = os.path.join(BROWSER_FOLDER, "lookup_tables")
    COMBINED_RESULTS_FOLDER = os.path.join(BROWSER_FOLDER, "combined_results")
    SEPARATED_RESULTS_FOLDER = os.path.join(BROWSER_FOLDER, "separated_results")
    INTERMEDIATE_FOLDER = "intermediate"

    # Ensure directories exist
    os.makedirs(QC_FOLDER, exist_ok=True)
    os.makedirs(BROWSER_FOLDER, exist_ok=True)
    os.makedirs(LOOKUP_TABLES_FOLDER, exist_ok=True)
    os.makedirs(COMBINED_RESULTS_FOLDER, exist_ok=True)
    os.makedirs(SEPARATED_RESULTS_FOLDER, exist_ok=True)
    os.makedirs(INTERMEDIATE_FOLDER, exist_ok=True)

    # Define file categories
    qc_files = {
        "pca_plot.pdf",
        "gene_diversity.csv.gz",
        "isoform_diversity.csv.gz",
        "SRSF6.csv.gz",
        "SRSF6_exonic_proportion.pdf",
        "isoform_diversity.pdf",
        "gene_diversity.pdf"
    }

    lookup_table_files = {
        "sample_gene_rankings_lookup_table.csv.gz",
        "gene_coverage_lookup_table.csv.gz"
    }

    combined_results_files = {
        "merged_ranked_gene_with_phenotype.csv.gz",
        "merged_ranked_isoform_with_phenotype.csv.gz"
    }

    separated_results_files = {
        "Cyclo_Allelic_Imbalance_gene_top_ranked_data.csv.gz",
        "Cyclo_GOE_gene_top_ranked_data.csv.gz",
        "Cyclo_GOE_isoform_top_ranked_data.csv.gz",
        "NMD_gene_top_ranked_data.csv.gz",
        "NMD_isoform_top_ranked_data.csv.gz",
        "NMD_rare_steady_state_transcript_gene_top_ranked_data.csv.gz",
        "Noncyclo_GOE_gene_top_ranked_data.csv.gz",
        "Noncyclo_GOE_isoform_top_ranked_data.csv.gz",
        "Noncyclo_LOE_gene_top_ranked_data.csv.gz",
        "Noncyclo_LOE_isoform_top_ranked_data.csv.gz",
        "Nonyclo_Allelic_Imbalance_gene_top_ranked_data.csv.gz"
    }

    browser_files = lookup_table_files | combined_results_files | separated_results_files

    # Get all `.gz` files in the current directory
    all_gz_files = {f for f in os.listdir() if f.endswith(".gz")}

    # Find files that should go into intermediate (everything not in QC or Browser)
    intermediate_files = all_gz_files - qc_files - browser_files

    # Function to move files
    def move_files(file_list, destination_folder):
        for file in file_list:
            if os.path.exists(file):  # Ensure the file exists before moving
                shutil.move(file, os.path.join(destination_folder, file))
            else:
                print(f"Warning: {file} not found, skipping.")

    # Move files to their respective folders
    move_files(qc_files, QC_FOLDER)
    move_files(lookup_table_files, LOOKUP_TABLES_FOLDER)
    move_files(combined_results_files, COMBINED_RESULTS_FOLDER)
    move_files(separated_results_files, SEPARATED_RESULTS_FOLDER)
    move_files(intermediate_files, INTERMEDIATE_FOLDER)

    print("File organization complete!")


if __name__ == "__main__":
    main()