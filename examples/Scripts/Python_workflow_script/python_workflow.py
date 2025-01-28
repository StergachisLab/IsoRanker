# Import required libraries
import pandas as pd
import matplotlib.pyplot as plt
import os

# Import functions from the package
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
    calculate_ranks_for_sample,
    create_expression_matrix,
    create_long_format,
    process_hypothesis_test
)

#########################################################
# Define input files
#########################################################

# Provides every single read, which sample each read came from, and the isoform asssociated with each read
read_stat_path = "IsoRanker/examples/Input/read_stat.txt"

# Assigns samples to patients, cyclo/non-cyclo, and haplotypes
sample_info_path = "IsoRanker/examples/Input/Sample_info.csv"
sample_info = pd.read_csv(sample_info_path)

# Provides pigeon annotation for isoforms
classification_path = "IsoRanker/examples/Input/pigeon_classification.txt"
classification_data = pd.read_csv(classification_path, sep="\t")

# Provide omim information for genes
genemap_path = "IsoRanker/examples/Input/genemap2.txt"
# Read the file, skipping the first 3 rows
genemap = pd.read_csv(genemap_path, sep='\t', skiprows=3)
genemap = genemap[genemap['Approved Gene Symbol'].notnull()]

#########################################################
# Use input files to create expression matrix
#########################################################

# Create the expression matrix and save it to a file
expression_matrix = create_expression_matrix(read_stat_path, output_file="expression_matrix.csv")

# Generate the long-format dataFrame and adding cyclo/noncyclo as well as haplotype labels and calculate TPM
long_format_df = create_long_format(expression_matrix, sample_info)

# Filter input
long_format_df = filter_based_on_counts(long_format_df, count_threshold=10, group_col='Isoform')

# Select only the 'isoform' and 'associated_gene' columns from classification_data
classification_subset = classification_data[['isoform', 'associated_gene']]
# Merge the classification subset with the long_format_df
long_format_annotated = long_format_df.merge(
    classification_subset,
    left_on="Isoform",   # Match isoform IDs in long_format_df
    right_on="isoform",  # Match isoform IDs in classification_subset
    how="left"           # Keep all rows from long_format_df, even if there's no match in classification_subset
).drop(columns=["isoform"])  # Drop redundant 'isoform' column from classification_subset

long_format_annotated.to_csv("long_format_annotated.csv", index=False)

#########################################################
# Use long format of expression matrix to calculate test statistic, z-scores, and ranks
#########################################################

# Calculate for all hypothesis tests
test_stat_funcs = [
    ("NMD", NMD_test_statistic),
    ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
    ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
    ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE),
    ("NMD_rare_steady_state_transcript", NMD_rare_steady_state_transcript)
]

# Gene level
for test_name, test_func in test_stat_funcs:
        print(f"Processing test statistic: {test_name}")

        # Apply the process_hypothesis_test function
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated, 
            group_col='Isoform', 
            test_statistic_func=test_func, 
            gene_group_col='associated_gene', 
            gene_level=True, 
            bin_proportion=0.01)

        # Filter for Ranks ≤ 1000
        rank_columns = [col for col in ranked_data.columns if col.startswith('rank_top_')]
        filtered_ranked_data = ranked_data[ranked_data[rank_columns].le(1000).any(axis=1)]

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
        output_dir = ""
        output_file = os.path.join(output_dir, f"{test_name}_gene_top_ranked_data.csv")
        filtered_ranked_data.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")


test_stat_funcs = [
    ("NMD", NMD_test_statistic),
    ("Noncyclo_LOE", Noncyclo_Expression_Outlier_LOE),
    ("Noncyclo_GOE", Noncyclo_Expression_Outlier_GOE),
    ("Cyclo_GOE", Cyclo_Expression_Outlier_GOE)
]

# Isoform level
for test_name, test_func in test_stat_funcs:
        print(f"Processing test statistic: {test_name}")

        # Apply the process_hypothesis_test function
        ranked_data = process_hypothesis_test(
            filtered_data=long_format_annotated, 
            group_col='Isoform', 
            test_statistic_func=test_func, 
            gene_group_col='associated_gene', 
            gene_level=False, 
            bin_proportion=0.01)

        # Filter for Ranks ≤ 1000
        rank_columns = [col for col in ranked_data.columns if col.startswith('rank_top_')]
        filtered_ranked_data = ranked_data[ranked_data[rank_columns].le(1000).any(axis=1)]

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
        output_dir = ""
        output_file = os.path.join(output_dir, f"{test_name}_isoform_top_ranked_data.csv")
        filtered_ranked_data.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")