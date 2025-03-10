import os
import sys
import gzip
import pandas as pd
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def process_and_plot_pca(df, output_pdf="pca_plot.pdf", grouping_col = "associated_gene"):

    df_filtered = filter_based_on_counts(df, count_threshold=10, group_col=grouping_col)

    # Group by Sample and associated_gene and sum Cyclo_TPM and Noncyclo_TPM
    grouped_df = df_filtered.groupby(["Sample", grouping_col])[["Cyclo_TPM", "Noncyclo_TPM"]].sum().reset_index()
    
    # Create two separate DataFrames for Cyclo and Noncyclo
    cyclo_df = grouped_df[["Sample", grouping_col, "Cyclo_TPM"]].copy()
    noncyclo_df = grouped_df[["Sample", grouping_col, "Noncyclo_TPM"]].copy()
    
    # Rename columns
    cyclo_df.columns = ["Sample", grouping_col, "TPM"]
    noncyclo_df.columns = ["Sample", grouping_col, "TPM"]
    
    # Append _Cyclo and _Noncyclo to Sample names
    cyclo_df["Sample"] = cyclo_df["Sample"] + "_Cyclo"
    noncyclo_df["Sample"] = noncyclo_df["Sample"] + "_Noncyclo"
    
    # Concatenate the two DataFrames
    final_df = pd.concat([cyclo_df, noncyclo_df], ignore_index=True)
    
    # Pivot table to have Samples as rows and associated_gene as columns
    pivot_df = final_df.pivot(index="Sample", columns=grouping_col, values="TPM").fillna(0)

    # Standardize the data to match R's prcomp(center=TRUE, scale=TRUE)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pivot_df)
    
    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)
    
    # Create DataFrame with PC1 and PC2
    pca_df = pd.DataFrame(principal_components, columns=["PC1", "PC2"], index=pivot_df.index).reset_index()
    
    # Add Cyclo/Noncyclo category
    pca_df["Condition"] = pca_df["Sample"].apply(lambda x: "Cyclo" if "_Cyclo" in x else "Noncyclo")
    
    # Define color mapping
    color_mapping = {"Cyclo": "red", "Noncyclo": "blue"}
    
    # Plot PCA results
    plt.figure(figsize=(8, 6))
    ax = sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="Condition", style="Condition", s=100,
        palette=color_mapping, markers={"Cyclo": "o", "Noncyclo": "o"}
    )
    
    # Add labels to each point
    for i, row in pca_df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["Sample"], fontsize=5, ha='right', va='bottom')
    
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of Samples")
    plt.legend()
    
    # Save the plot to a PDF
    plt.savefig(output_pdf)
    plt.close()
    
    return pca_df


def analyze_isoforms(df, output_file, grouping_column):
    """
    Analyzes isoform counts for each sample, computing the number of unique genes
    that pass various read count thresholds.
    
    Parameters:
    - df (pd.DataFrame): Input DataFrame with columns ['Sample', grouping_column, 'cyclo_count', 'noncyclo_count']
    - output_file (str): Path to save the results as a compressed CSV.
    - grouping_column (str): Column name used to group the counts (e.g., 'associated_gene' or 'Isoform').
    
    Returns:
    - None: Saves a compressed CSV file.
    """

    # Step 1: Collapse data by Sample and grouping_column by summing counts
    grouped_df = df.groupby(["Sample", grouping_column], as_index=False)[["cyclo_count", "noncyclo_count"]].sum()

    # Initialize results list
    results = []

    # Step 2: Iterate through each unique Sample
    for sample in grouped_df['Sample'].unique():
        sample_data = grouped_df[grouped_df['Sample'] == sample]
        
        # Cyclo analysis after collapsing
        cyclo_data = sample_data[sample_data['cyclo_count'] > 0]
        cyclo_total_reads = cyclo_data['cyclo_count'].sum()
        cyclo_unique = cyclo_data[grouping_column].nunique()
        cyclo_gt1 = cyclo_data[cyclo_data['cyclo_count'] > 1][grouping_column].nunique()
        cyclo_gt10 = cyclo_data[cyclo_data['cyclo_count'] > 10][grouping_column].nunique()
        cyclo_gt100 = cyclo_data[cyclo_data['cyclo_count'] > 100][grouping_column].nunique()
        cyclo_gt1000 = cyclo_data[cyclo_data['cyclo_count'] > 1000][grouping_column].nunique()
        
        # Noncyclo analysis after collapsing
        noncyclo_data = sample_data[sample_data['noncyclo_count'] > 0]
        noncyclo_total_reads = noncyclo_data['noncyclo_count'].sum()
        noncyclo_unique = noncyclo_data[grouping_column].nunique()
        noncyclo_gt1 = noncyclo_data[noncyclo_data['noncyclo_count'] > 1][grouping_column].nunique()
        noncyclo_gt10 = noncyclo_data[noncyclo_data['noncyclo_count'] > 10][grouping_column].nunique()
        noncyclo_gt100 = noncyclo_data[noncyclo_data['noncyclo_count'] > 100][grouping_column].nunique()
        noncyclo_gt1000 = noncyclo_data[noncyclo_data['noncyclo_count'] > 1000][grouping_column].nunique()
        
        # Append results
        results.append([
            sample, cyclo_unique, cyclo_total_reads, cyclo_gt1,
            cyclo_gt10, cyclo_gt100, cyclo_gt1000,
            noncyclo_unique, noncyclo_total_reads, noncyclo_gt1,
            noncyclo_gt10, noncyclo_gt100, noncyclo_gt1000
        ])
    
    # Step 3: Create a DataFrame for output
    output_df = pd.DataFrame(results, columns=[
        'Sample',
        f'Cyclo Unique {grouping_column}', 'Cyclo Total Reads', f'Cyclo Unique {grouping_column} >1',
        f'Cyclo Unique {grouping_column} >10', f'Cyclo Unique {grouping_column} >100', f'Cyclo Unique {grouping_column} >1000',
        f'Noncyclo Unique {grouping_column}', 'Noncyclo Total Reads', f'Noncyclo Unique {grouping_column} >1',
        f'Noncyclo Unique {grouping_column} >10', f'Noncyclo Unique {grouping_column} >100', f'Noncyclo Unique {grouping_column} >1000'
    ])
    
    # Step 4: Save to compressed CSV
    output_df.to_csv(output_file, index=False, compression="gzip")


def process_pileup(df, reference_fasta, chromosome, position, output_file):
    unique_bams = df.drop_duplicates(subset=['BAM_FILE'])[['SAMPLE', 'ID', 'CYCLO_NONCYCLO', 'BAM_FILE']]
    
    results = [["Source", "Chromosome", "Position", "Reference_Base", "Original_Read_Depth", "Adjusted_Read_Depth", "Exonic_Proportion", "Read_Bases", "Base_Qualities"]]
    
    for _, row in unique_bams.iterrows():
        source = f"{row['SAMPLE']}_{row['ID']}_{row['CYCLO_NONCYCLO']}"
        bam_file = row['BAM_FILE']
        print(f"Processing BAM file: {bam_file} from source: {source}")
        
        if pd.notna(bam_file) and os.path.exists(bam_file):
            try:
                samfile = pysam.AlignmentFile(bam_file, "rb")
                pileup_data = []
                
                for pileupcolumn in samfile.pileup(chromosome, position - 1, position, truncate=True, fasta=pysam.FastaFile(reference_fasta)):
                    ref_base = pileupcolumn.reference_pos
                    read_depth = pileupcolumn.n
                    read_bases = ''.join([pileupread.alignment.query_sequence[pileupread.query_position] if pileupread.query_position is not None else '' for pileupread in pileupcolumn.pileups])
                    base_qualities = ''.join([chr(pileupread.alignment.query_qualities[pileupread.query_position] + 33) if pileupread.query_position is not None else '' for pileupread in pileupcolumn.pileups])
                    
                    skip_count = read_bases.count('<') + read_bases.count('>')
                    adjusted_read_depth = read_depth - skip_count
                    exonic_proportion = adjusted_read_depth / read_depth if read_depth > 0 else 0
                    
                    pileup_data.append([source, chromosome, position, ref_base, read_depth, adjusted_read_depth, round(exonic_proportion, 2), read_bases, base_qualities])
                
                samfile.close()
                results.extend(pileup_data if pileup_data else [[source, chromosome, position, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"]])
            except Exception as e:
                results.append([source, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", str(e)])
        else:
            print(f"{bam_file} not found")
            results.append([source, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"])
    
    df_output = pd.DataFrame(results[1:], columns=results[0])
    df_output.to_csv(output_file, sep='\t', index=False)
    print(f"Processing complete. Results saved to {output_file}")
