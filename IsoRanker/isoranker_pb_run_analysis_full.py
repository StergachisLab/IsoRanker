#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import shutil
import argparse

from IsoRanker import (
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
    merge_tsvs_by_keyword, 
    process_vep_vcf, 
    merge_haplotype_data, 
    process_phenotype_data,
    process_and_plot_pca,
    analyze_isoforms,
    process_pileup,
    split_fusion_genes,
    write_sample_gene_lists,
    passes_max_af_filter,
    filter_multiple_vcfs
)

def main():

    ################################################
    # Parse command-line arguments
    ################################################
    parser = argparse.ArgumentParser(description="Run analysis pipeline.")
    parser.add_argument("--read_stat_path", required=True, help="Path to the read_stat.txt file.")
    parser.add_argument("--sample_info_path", required=True, help="Path to the sample info tsv file.")
    parser.add_argument("--classification_path", required=True, help="Path to the pigeon classification file.")
    parser.add_argument("--genemap_path", required=True, help="Path to the genemap2.txt file.")
    parser.add_argument("--hpo_file_path", required=True, help="Path to the file that matches HPO terms to OMIM. Download from here: https://hpo.jax.org/data/annotations.")
    parser.add_argument("--probands_file_path", required=True, help="Path to the file that contains proband HPO terms.")
    parser.add_argument("--reference_fasta_path", required=True, help="Path to the hg38 reference fasta path.")
    parser.add_argument("--final_output_dir", required=True, help="Directory to save the output files.")
    parser.add_argument("--gtf_path_input", required=True, help="Path to gencode reference file.")


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
    final_output_dir = args.final_output_dir
    gtf_path_input = args.gtf_path_input

    output_dir = "."

    #output_dir = final_output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    ################################################
    # Load input files
    ################################################
    print("Reading input files", flush=True)
    sample_info = pd.read_csv(sample_info_path, sep="\t")
    classification_data = pd.read_csv(classification_path, sep="\t")
    genemap = pd.read_csv(genemap_path, sep='\t', skiprows=3) # Read the file, skipping the first 3 rows
    genemap = genemap[genemap['Approved Gene Symbol'].notnull()]


    ################################################
    # Update input files with haplotype information
    ################################################
    print("Updating input files with haplotype information", flush=True)
    update_files_with_haplotype_info(sample_info_path, read_stat_path, output_dir)
    read_stat_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
    sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")
    sample_info = pd.read_csv(sample_info_path, compression = "gzip", sep="\t")

    ################################################
    # Create the isoform-level expression matrix
    ################################################
    print("Creating expression matrix", flush=True)
    expression_matrix = create_expression_matrix(read_stat_path, output_file=os.path.join(output_dir, "expression_matrix.tsv.gz"))

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

    long_format_annotated.to_csv(os.path.join(output_dir, "long_format_annotated.tsv.gz"), index=False, compression = "gzip", sep="\t")

    ################################################
    # Create the gene-level expression matrix
    ################################################
    # This matrix it not used in any of the analysis, but might be a helpful intermediate file to output. 

    print("Creating gene expression matrix", flush=True)

    # Load the expression matrix (rows: isoforms, columns: samples)
    #expr_df = pd.read_csv("/mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker_testing/Command_line/4.9.25/Output/intermediate/expression_matrix.tsv.gz", sep="\t", index_col=0, compression="gzip")
    
    # Merge mapping with expression data
    gene_expr_df = classification_subset.merge(expression_matrix, left_on="isoform", right_index=True)

    # Group by gene and sum expression values across isoforms
    gene_expr_df = gene_expr_df.drop(columns=["isoform"]).groupby("associated_gene").sum()

    # Save to output file
    gene_expr_df.to_csv("gene_expression_matrix.tsv.gz", sep="\t", compression="gzip")



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
        ("Noncyclo_Allelic_Imbalance", Noncyclo_Allelic_Imbalance),
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


        # Add OMIM data to genes
        ranked_data = ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],  # Select relevant columns from genemap
            how='left',  # Perform a left join to keep all rows from filtered_ranked_data
            left_on='associated_gene',  # Column in filtered_ranked_data to join on
            right_on='Approved Gene Symbol'  # Column in genemap to join on
        )
        # Drop the 'Approved Gene Name' column if it is no longer needed
        ranked_data = ranked_data.drop(columns=['Approved Gene Symbol'])

        ranked_data = split_fusion_genes(ranked_data)

        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        # Save the results to a tsv file
        output_file = os.path.join(output_dir, f"{test_name}_gene_top_ranked_data.tsv.gz")
        filtered_ranked_data.to_csv(output_file, index=False, compression = "gzip", sep="\t")
        print(f"Results saved to {output_file}", flush=True)

        # Save the full results to a tsv file. 
        output_file = os.path.join(output_dir, f"{test_name}_gene_full_ranked_data.tsv.gz")
        ranked_data.to_csv(output_file, index=False, compression = "gzip", sep="\t")
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
            filter_count_threshold=5)

        # Add OMIM data to genes
        ranked_data = ranked_data.merge(
            genemap[['Approved Gene Symbol', 'Phenotypes']],  # Select relevant columns from genemap
            how='left',  # Perform a left join to keep all rows from filtered_ranked_data
            left_on='associated_gene',  # Column in filtered_ranked_data to join on
            right_on='Approved Gene Symbol'  # Column in genemap to join on
        )
        # Drop the 'Approved Gene Name' column if it is no longer needed
        ranked_data = ranked_data.drop(columns=['Approved Gene Symbol'])

        ranked_data = split_fusion_genes(ranked_data)

        filtered_ranked_data = ranked_data[ranked_data["rank_top_99_5_percentile"] <= 25]

        # Save the results to a tsv file
        output_file = os.path.join(output_dir, f"{test_name}_isoform_top_ranked_data.tsv.gz")
        filtered_ranked_data.to_csv(output_file, index=False, compression = "gzip", sep="\t")
        print(f"Results saved to {output_file}", flush=True)

        # Save the full results to a tsv file. 
        output_file = os.path.join(output_dir, f"{test_name}_isoform_full_ranked_data.tsv.gz")
        ranked_data.to_csv(output_file, index=False, compression = "gzip", sep="\t")
        print(f"Results saved to {output_file}", flush=True)



    ################################################
    # Combine output files
    ################################################

    # Combine the ranked files

    print("Combining output files", flush=True)

    # Isoform
    keyword = "_isoform_top_" 
    output_tsv = os.path.join(output_dir, f"merged_ranked_isoform.tsv.gz")
    merge_tsvs_by_keyword(output_dir, keyword, output_tsv)

    # Gene
    keyword = "_gene_top_"
    output_tsv = os.path.join(output_dir, f"merged_ranked_gene.tsv.gz")
    merge_tsvs_by_keyword(output_dir, keyword, output_tsv)


    ################################################
    # Add patient phenotype information
    ################################################
    print("Adding patient phenotype information to combined output files", flush=True)

    all_comparisons, all_comparisons_long = process_phenotype_data(hpo_file_path, genemap_path, probands_file_path)

    # Master file gene
    master_file = pd.read_csv("merged_ranked_gene.tsv.gz", compression = "gzip", sep="\t")

    # Merge `master_file` with `all_comparisons_long` based on `Sample` and gene name
    merged_data = master_file.merge(
        all_comparisons_long,
        left_on=['Sample', 'associated_gene'],  # Match Sample and associated gene
        right_on=['Sample', 'Approved Gene Symbol'],  # Match with Approved Gene Symbol
        how='left'  # Keep all rows from master_file
    )

    # Save merged output
    merged_data.to_csv("merged_ranked_gene_with_phenotype.tsv.gz", index=False, compression = "gzip", sep="\t")


    # Master file isoform
    master_file = pd.read_csv("merged_ranked_isoform.tsv.gz", compression = "gzip", sep="\t")

    # Merge `master_file` with `all_comparisons_long` based on `Sample` and gene name
    merged_data = master_file.merge(
        all_comparisons_long,
        left_on=['Sample', 'associated_gene'],  # Match Sample and associated gene
        right_on=['Sample', 'Approved Gene Symbol'],  # Match with Approved Gene Symbol
        how='left'  # Keep all rows from master_file
    )

    # Save merged output
    merged_data.to_csv("merged_ranked_isoform_with_phenotype.tsv.gz", index=False, compression = "gzip", sep="\t")

    ################################################
    # Add genetic variant information
    ################################################

    ### Make gene list
    gene_file_input = "merged_ranked_gene_with_phenotype.tsv.gz"
    isoform_file_input = "merged_ranked_isoform_with_phenotype.tsv.gz"

    write_sample_gene_lists(
        gene_file=gene_file_input,
        isoform_file=isoform_file_input,
        output_dir="gene_lists_by_sample"
    )

    ### Subset vcfs based on gene lists and frequency
    filter_multiple_vcfs(pair_list_path=sample_info_path, 
                         gene_list_dir="gene_lists_by_sample",
                         flank=50, 
                         max_af_cutoff=0.01, 
                         gtf_path=gtf_path_input,
                         output_dir="subsetted_vcfs")


    ### Extract annotations and genotype into table format so that the information can be added to the IsoRanker tables
    input_vcf_dir = "subsetted_vcfs"
    output_base_dir = "variant_annotations_tables"

    os.makedirs(output_base_dir, exist_ok=True)

    for individual in sample_info['individual'].unique():
        vcf_path = os.path.join(input_vcf_dir, f"{individual}.isoranker_subsetted.vcf.gz")
        os.makedirs(output_base_dir, exist_ok=True)

        if os.path.exists(vcf_path):
            print(f"Processing {individual}...")
            process_vep_vcf(vcf_path, output_base_dir, individual)
        else:
            print(f"VCF not found for {individual}: {vcf_path}")


    ### Add annotations to isoranker tables
    gene_level_merged_df = pd.read_csv(gene_file_input, sep="\t")

    isoform_level_merged_df = pd.read_csv(isoform_file_input, sep="\t")

    # Placeholder for merged haplotype data
    haplotype_entries = []

    # For each sample in gene_level_merged_df
    for sample in gene_level_merged_df['Sample'].unique():
        tsv_path = f"variant_annotations_tables/{sample}_gene_haplotype_split.tsv"
        
        if os.path.exists(tsv_path):
            # Load haplotype data
            haplo_df = pd.read_csv(tsv_path, sep="\t")
            haplo_df["Sample"] = sample  # tag with sample for merge
            
            haplotype_entries.append(haplo_df)
        else:
            print(f"Haplotype TSV not found for {sample}")

    # Combine all haplotype info
    all_haplo = pd.concat(haplotype_entries, ignore_index=True)

    # Rename for clarity to match gene_level_merged_df
    all_haplo.rename(columns={"Gene": "associated_gene"}, inplace=True)

    # Merge into your existing gene-level dataframe
    gene_level_merged_df_variant_annotated = gene_level_merged_df.merge(
        all_haplo, how="left", on=["Sample", "associated_gene"]
    )
    gene_level_merged_df_variant_annotated.to_csv('merged_ranked_gene_with_phenotype_with_variant.tsv.gz', index=False, sep="\t", compression="gzip")

    # Merge into your existing gene-level dataframe
    isoform_level_merged_df_variant_annotated = isoform_level_merged_df.merge(
        all_haplo, how="left", on=["Sample", "associated_gene"]
    )
    isoform_level_merged_df_variant_annotated.to_csv('merged_ranked_isoform_with_phenotype_with_variant.tsv.gz', index=False, sep="\t", compression="gzip")

    classification_data = classification_data[['isoform', 'structural_category', 'subcategory']]

    isoform_level_merged_df_variant_annotated_pigeon = isoform_level_merged_df_variant_annotated.merge(
        classification_data,
        left_on="Isoform",   # Match isoform IDs in long_format_df
        right_on="isoform",  # Match isoform IDs in classification_subset
        how="left"           # Keep all rows from long_format_df, even if there's no match in classification_subset
    ).drop(columns=["isoform"])  # Drop redundant 'isoform' column from classification_subset

    isoform_level_merged_df_variant_annotated_pigeon.to_csv('merged_ranked_isoform_with_phenotype_with_variant_with_pigeon.tsv.gz', index=False, sep="\t", compression="gzip")

    # Variants are annotated in this format: # (CHROM_POS_REF/ALT,MAX_AF,ClinVar_CLNSIG,ClinVar_CLNDN,Consequence,SpliceAI_Score,SpliceAI_Flag,Symbol,CADD_PHRED,CADD_Flag)

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

    # Save as a compressed tsv (gzip format)
    sample_gene_rankings_lookup_table.to_csv("sample_gene_rankings_lookup_table.tsv.gz", index=False, compression="gzip", sep="\t")

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

    gene_coverage_lookup_table.to_csv("gene_coverage_lookup_table.tsv.gz", index=False, compression="gzip", sep="\t")


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

    analyze_isoforms(long_format_annotated, "gene_diversity.tsv.gz", "associated_gene")

    # Plotting:

    df = pd.read_csv("gene_diversity.tsv.gz", compression = "gzip", sep="\t")

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

    analyze_isoforms(long_format_annotated, "isoform_diversity.tsv.gz", "Isoform")

    df = pd.read_csv("isoform_diversity.tsv.gz", compression = "gzip", sep="\t")

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

    process_pileup(df=sample_info, reference_fasta= reference_fasta_path, chromosome="chr20", position=43459200, output_file="SRSF6.tsv.gz")

    # Plotting:

    SRSF6_df = pd.read_csv("SRSF6.tsv.gz", compression = "gzip", sep="\t")

    # Sort DataFrame by source
    df_sorted = SRSF6_df.sort_values(by="Source")

    # Define colors based on condition
    colors = ["red" if "_cyclo" in source else "blue" for source in df_sorted["Source"]]

    # Plot bar chart
    plt.figure(figsize=(15, 9))
    sns.barplot(data=df_sorted, x="Source", y="Exonic_Proportion", edgecolor="black", palette=colors)

    # Labeling and formatting
    plt.xlabel("Source")
    plt.ylabel("Exonic Proportion")
    plt.title("Bar Plot of Exonic Proportion by Source")
    plt.xticks(rotation=45, fontsize=7, ha="right")  # Adjust alignment to prevent cutoff

    plt.tight_layout()

    # Save to PDF
    plt.savefig("SRSF6_exonic_proportion.pdf", format="pdf")

    ###################################
    # SRSF6 deltas
    ###################################

    # Compute deltas between every pair of rows (0-1, 2-3, etc.)
    deltas = []
    for i in range(0, len(df_sorted), 2):
        source = df_sorted.iloc[i]["Source"] + "_" + df_sorted.iloc[i + 1]["Source"]
        delta = (df_sorted.iloc[i]["Exonic_Proportion"] - df_sorted.iloc[i + 1]["Exonic_Proportion"])/df_sorted.iloc[i]["Exonic_Proportion"] 
        deltas.append({"Source": source, "Exonic_Delta": delta})

    # Create new DataFrame
    delta_df = pd.DataFrame(deltas)

    # Optional: Set a style
    sns.set(style="whitegrid")

    # Create the boxplot
    plt.figure(figsize=(8, 5))
    ax = sns.boxplot(data=delta_df, y="Exonic_Delta", color="lightgray", width=0.3)

    # Overlay the swarmplot
    sns.swarmplot(data=delta_df, y="Exonic_Delta", color="black", size=8)

    # Set y-axis limits
    ax.set_ylim(-1, 1.1)

    # Optional: Add x-axis label
    ax.set_xlabel("")
    ax.set_ylabel("Normalized Exonic Delta")
    ax.set_title("Normalized Exonic Delta per Patient")
    plt.tight_layout()

    # Save to PDF
    plt.savefig("SRSF6_exonic_proportion_normalized_delta.pdf", format="pdf")


    ###################################
    # Organize output files
    ###################################

    # Define folder names
    #OUTPUT_FOLDER = "Output"
    OUTPUT_FOLDER = final_output_dir
    QC_FOLDER = os.path.join(OUTPUT_FOLDER, "qc")
    BROWSER_FOLDER = os.path.join(OUTPUT_FOLDER, "browser")
    LOOKUP_TABLES_FOLDER = os.path.join(BROWSER_FOLDER, "lookup_tables")
    COMBINED_RESULTS_FOLDER = os.path.join(BROWSER_FOLDER, "combined_results")
    SEPARATED_RESULTS_FOLDER = os.path.join(BROWSER_FOLDER, "separated_results")
    INTERMEDIATE_FOLDER = os.path.join(OUTPUT_FOLDER, "intermediate")

    # Ensure directories exist
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    os.makedirs(QC_FOLDER, exist_ok=True)
    os.makedirs(BROWSER_FOLDER, exist_ok=True)
    os.makedirs(LOOKUP_TABLES_FOLDER, exist_ok=True)
    os.makedirs(COMBINED_RESULTS_FOLDER, exist_ok=True)
    os.makedirs(SEPARATED_RESULTS_FOLDER, exist_ok=True)
    os.makedirs(INTERMEDIATE_FOLDER, exist_ok=True)

    # Define file categories
    qc_files = {
        "pca_plot.pdf",
        "gene_diversity.tsv.gz",
        "isoform_diversity.tsv.gz",
        "SRSF6.tsv.gz",
        "SRSF6_exonic_proportion.pdf",
        "isoform_diversity.pdf",
        "gene_diversity.pdf",
        "SRSF6_exonic_proportion_normalized_delta.pdf"
    }

    lookup_table_files = {
        "sample_gene_rankings_lookup_table.tsv.gz",
        "gene_coverage_lookup_table.tsv.gz"
    }

    combined_results_files = {
        "merged_ranked_gene_with_phenotype.tsv.gz",
        "merged_ranked_isoform_with_phenotype.tsv.gz",
        "merged_ranked_gene_with_phenotype_with_variant.tsv.gz",
        "merged_ranked_isoform_with_phenotype_with_variant.tsv.gz",
        "merged_ranked_isoform_with_phenotype_with_variant_with_pigeon.tsv.gz"
    }

    separated_results_files = {
        "Cyclo_Allelic_Imbalance_gene_top_ranked_data.tsv.gz",
        "Cyclo_GOE_gene_top_ranked_data.tsv.gz",
        "Cyclo_GOE_isoform_top_ranked_data.tsv.gz",
        "NMD_gene_top_ranked_data.tsv.gz",
        "NMD_isoform_top_ranked_data.tsv.gz",
        "NMD_rare_steady_state_transcript_gene_top_ranked_data.tsv.gz",
        "Noncyclo_GOE_gene_top_ranked_data.tsv.gz",
        "Noncyclo_GOE_isoform_top_ranked_data.tsv.gz",
        "Noncyclo_LOE_gene_top_ranked_data.tsv.gz",
        "Noncyclo_LOE_isoform_top_ranked_data.tsv.gz",
        "Noncyclo_Allelic_Imbalance_gene_top_ranked_data.tsv.gz"
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
                print(f"Warning: {file} not found, skipping.", flush=True)

    # Move files to their respective folders
    move_files(qc_files, QC_FOLDER)
    move_files(lookup_table_files, LOOKUP_TABLES_FOLDER)
    move_files(combined_results_files, COMBINED_RESULTS_FOLDER)
    move_files(separated_results_files, SEPARATED_RESULTS_FOLDER)
    move_files(intermediate_files, INTERMEDIATE_FOLDER)


    if os.path.isdir("gene_lists_by_sample"):
        shutil.move("gene_lists_by_sample", INTERMEDIATE_FOLDER)

    if os.path.isdir("subsetted_vcfs"):
        shutil.move("subsetted_vcfs", INTERMEDIATE_FOLDER)

    if os.path.isdir("variant_annotations_tables"):
        shutil.move("variant_annotations_tables", INTERMEDIATE_FOLDER)


    print("File organization complete!", flush=True)


if __name__ == "__main__":
    main()