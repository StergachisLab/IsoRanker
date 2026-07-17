import os
import sys
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from IsoRanker import (
    filter_based_on_counts
)


def process_and_plot_pca(df, output_pdf="pca_plot.pdf", grouping_col="associated_gene", count_filter=10):
    """
    Performs principal component analysis (PCA) on transcript expression data and generates a PCA plot.

    Parameters:
    - df (pd.DataFrame): Input DataFrame containing transcript expression data with columns:
        - 'Sample' (str): Identifier for each sample.
        - grouping_col (str): Column used for grouping (e.g., 'associated_gene' or 'Isoform').
        - 'Cyclo_TPM' (float): TPM values for cycloheximide-treated samples.
        - 'Noncyclo_TPM' (float): TPM values for untreated samples.
    - output_pdf (str, optional): Path to save the PCA plot as a PDF. Default is "pca_plot.pdf".
    - grouping_col (str, optional): Column name used to group TPM values. Default is "associated_gene".

    Returns:
    - pd.DataFrame: A DataFrame containing PCA results with the following columns:
        - 'Sample' (str): Sample identifiers.
        - 'PC1' (float): First principal component values.
        - 'PC2' (float): Second principal component values.
        - 'Condition' (str): 'Cyclo' or 'Noncyclo', based on the sample name.

    Notes:
    - The function filters data to include only genes with at least 10 counts using `filter_based_on_counts()`.
    - TPM values are grouped by 'Sample' and `grouping_col`, then standardized using `StandardScaler()`.
    - PCA is performed with two components.
    - The PCA plot is colored by treatment condition (Cyclo in red, Noncyclo in blue).
    - The plot includes variance explained on axis labels and is saved as a PDF.
    """

    df_filtered = filter_based_on_counts(df, count_threshold=count_filter, group_col=grouping_col)
    

    # Group by Sample and grouping_col and sum TPM values
    grouped_df = df_filtered.groupby(["Sample", grouping_col])[["Cyclo_TPM", "Noncyclo_TPM"]].sum().reset_index()

    # Create DataFrames for Cyclo and Noncyclo
    cyclo_df = grouped_df[["Sample", grouping_col, "Cyclo_TPM"]].copy()
    noncyclo_df = grouped_df[["Sample", grouping_col, "Noncyclo_TPM"]].copy()

    # Rename columns
    cyclo_df.columns = ["Sample", grouping_col, "TPM"]
    noncyclo_df.columns = ["Sample", grouping_col, "TPM"]

    # Label samples
    cyclo_df["Sample"] = cyclo_df["Sample"] + "_Cyclo"
    noncyclo_df["Sample"] = noncyclo_df["Sample"] + "_Noncyclo"

    # Concatenate
    final_df = pd.concat([cyclo_df, noncyclo_df], ignore_index=True)

    # Pivot: rows = Sample, columns = genes/isoforms, values = TPM
    pivot_df = final_df.pivot(index="Sample", columns=grouping_col, values="TPM").fillna(0)

    # Remove samples with all-zero TPMs
    pivot_df = pivot_df.loc[pivot_df.any(axis=1)]

    # Standardize
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pivot_df)

    # PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)

    # Create PCA result DataFrame
    pca_df = pd.DataFrame(principal_components, columns=["PC1", "PC2"], index=pivot_df.index).reset_index()
    pca_df["Condition"] = pca_df["Sample"].apply(lambda x: "Cyclo" if "_Cyclo" in x else "Noncyclo")

    # Plotting
    color_mapping = {"Cyclo": "red", "Noncyclo": "blue"}
    plt.figure(figsize=(8, 6))
    ax = sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="Condition", style="Condition", s=100,
        palette=color_mapping, markers={"Cyclo": "o", "Noncyclo": "o"}
    )

    # Add sample labels
    for _, row in pca_df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["Sample"], fontsize=5, ha='right', va='bottom')

    # Axis labels with explained variance
    explained_var = pca.explained_variance_ratio_ * 100
    plt.xlabel(f"PC1 ({explained_var[0]:.2f}%)")
    plt.ylabel(f"PC2 ({explained_var[1]:.2f}%)")
    plt.title("PCA of Samples")
    plt.legend()

    # Save
    plt.savefig(output_pdf)
    plt.close()

    return pca_df



def analyze_isoforms(df, output_file, grouping_column):
    """
    Analyzes isoform counts for each sample, computing the number of unique genes
    that pass various read count thresholds.
    
    Parameters:
    - df (pd.DataFrame): Input DataFrame with columns ['Sample', grouping_column, 'cyclo_count', 'noncyclo_count']
    - output_file (str): Path to save the results as a compressed TSV.
    - grouping_column (str): Column name used to group the counts (e.g., 'associated_gene' or 'Isoform').
    
    Returns:
    - None: Saves a compressed TSV file.
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
    
    # Step 4: Save to compressed TSV
    output_df.to_csv(output_file, index=False, compression="gzip", sep="\t")


def process_pileup(df, reference_fasta, chromosome, position, output_file):
    """
    Processes pileup data from BAM files for a specific genomic position, extracting read depth, 
    nucleotide composition, and base qualities.

    Parameters:
    - df (pd.DataFrame): Input DataFrame containing sample and BAM file information. 
                         Expected columns: ['sample', 'patient', 'cyclo', 'bam_file']
    - reference_fasta (str): Path to the reference FASTA file, used to retrieve the reference base.
    - chromosome (str): Chromosome name (e.g., 'chr1', '2', 'X') for the pileup analysis.
    - position (int): genomic position at which pileup data is collected.
    - output_file (str): Path to save the results as a compressed tsv (gzip format).

    Returns:
    - None: Saves a compressed tsv file with columns:
        - 'Source': Unique identifier combining 'sample', 'patient', and 'cyclo'.
        - 'Chromosome': Chromosome name.
        - 'Position': Genomic position.
        - 'Reference_Base': Reference nucleotide at this position.
        - 'Original_Read_Depth': Total number of reads covering the position.
        - 'Exon_Read_Depth': Number of reads containing an exonic nucleotide (A, C, T, or G).
        - 'Exonic_Proportion': Proportion of exonic reads (Exon_Read_Depth / Original_Read_Depth).
        - 'Read_Bases': String of bases observed at this position.
        - 'Base_Qualities': ASCII-encoded Phred quality scores.

    Notes:
    - Uses `pysam.AlignmentFile.pileup()` to extract pileup data from BAM files.
    - Reference bases are obtained from the provided FASTA file.
    - Reads that contain deletions or are soft-clipped at this position are ignored.
    - Assumes that **A, C, T, G** bases represent exonic reads.
    - If a BAM file does not exist, an "N/A" row is added to the results.
    - Errors encountered during processing are logged in the output file.
    """

    #unique_bams = df.drop_duplicates(subset=['bam_file'])[['sample', 'individual', 'condition', 'bam_file']]
    unique_bams = df.drop_duplicates(subset=['individual', 'condition'])[['sample', 'individual', 'condition', 'bam_file']]

    
    results = [["Source", "Chromosome", "Position", "Reference_Base", "Original_Read_Depth", "Exon_Read_Depth", "Exonic_Proportion", "Read_Bases", "Base_Qualities"]]
    
    for _, row in unique_bams.iterrows():
        source = f"{row['individual']}_{row['condition']}"
        bam_file = row['bam_file']
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
                    
                    exon_read_count = read_bases.count('A') + read_bases.count('C') + read_bases.count('T') + read_bases.count('G')
                    exonic_proportion = exon_read_count / read_depth if read_depth > 0 else 0
                    
                    pileup_data.append([source, chromosome, position, ref_base, read_depth, exon_read_count, round(exonic_proportion, 2), read_bases, base_qualities])
                
                samfile.close()
                results.extend(pileup_data if pileup_data else [[source, chromosome, position, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"]])
            except Exception as e:
                results.append([source, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", str(e)])
        else:
            print(f"{bam_file} not found")
            results.append([source, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"])
    
        
    df_output = pd.DataFrame(results[1:], columns=results[0])
    df_output.to_csv(output_file, index=False, compression = "gzip", sep="\t")
    print(f"Processing complete. Results saved to {output_file}")


def plot_haplotype_skew_noncyclo(
    data,
    output_table="haplotype_skew_summary_noncyclo_individual.tsv",
    output_figure="haplotype_skew_plot_noncyclo_jitter_individual.pdf",
    genes_of_interest=None,
    min_total=10,
    phased_prop_min=0.1,
    random_seed=42,
):
    """
    Calculate and plot noncyclo haplotype skew for selected genes.

    Haplotype skew is defined as::

        abs(H1 - H2) / (H1 + H2)

    A value is calculated only when both of the following are true:

    - H1 + H2 >= ``min_total``
    - (H1 + H2) / (H0 + H1 + H2) >= ``phased_prop_min``

    Point area is globally scaled by the phased-read proportion across all
    plotted sample-gene observations.

    Parameters
    ----------
    data : pandas.DataFrame or str or os.PathLike
        Input DataFrame or path to a tab-delimited table. Gzip-compressed input
        is detected by pandas from the filename extension.
    output_table : str or os.PathLike, optional
        Path for the tidy output table.
    output_figure : str or os.PathLike, optional
        Path for the PDF figure.
    genes_of_interest : sequence of str, optional
        Genes to plot. Defaults to the eight imprinted genes used in the
        original QC script.
    min_total : int, optional
        Minimum H1 + H2 count required to calculate skew.
    phased_prop_min : float, optional
        Minimum fraction of reads assigned to H1 or H2.
    random_seed : int or None, optional
        Seed for horizontal jitter. Use None for non-reproducible jitter.

    Returns
    -------
    pandas.DataFrame
        Tidy table containing counts, haplotype skew, and phased proportion.

    Raises
    ------
    ValueError
        If required columns are missing, thresholds are invalid, duplicate
        sample-gene rows are present, or no requested genes are found.
    """
    if genes_of_interest is None:
        genes_of_interest = [
            "ZNF597", "ZNF331", "MEG3", "NAP1L5",
            "PEG10", "PLAGL1", "ZDBF2", "DLK1",
        ]
    else:
        genes_of_interest = list(genes_of_interest)


    if min_total < 0:
        raise ValueError("min_total must be non-negative.")
    if not 0 <= phased_prop_min <= 1:
        raise ValueError("phased_prop_min must be between 0 and 1.")

    if isinstance(data, (str, os.PathLike)):
        df = pd.read_csv(data, sep="\t", compression="infer")
    elif isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        raise TypeError("data must be a pandas DataFrame or a file path.")

    sample_col = "Sample"
    gene_col = "associated_gene"
    h0_col = "H0_noncyclo_count"
    h1_col = "H1_noncyclo_count"
    h2_col = "H2_noncyclo_count"
    required_columns = {sample_col, gene_col, h0_col, h1_col, h2_col}
    missing_columns = sorted(required_columns.difference(df.columns))
    if missing_columns:
        raise ValueError(
            "Input is missing required column(s): " + ", ".join(missing_columns)
        )

    # Preserve the original ordering of Individual 1, Individual 2, etc.
    sample_num = pd.to_numeric(
        df[sample_col].astype(str).str.extract(r"Individual (\d+)", expand=False),
        errors="coerce",
    )
    sample_order = (
        df.assign(_sample_num=sample_num)
        .sort_values(["_sample_num", sample_col], na_position="last")[sample_col]
        .drop_duplicates()
        .tolist()
    )

    df = df[df[gene_col].isin(genes_of_interest)].copy()
    if df.empty:
        raise ValueError("None of the requested genes were found in the input data.")

    duplicate_mask = df.duplicated([sample_col, gene_col], keep=False)
    if duplicate_mask.any():
        duplicate_pairs = (
            df.loc[duplicate_mask, [sample_col, gene_col]]
            .drop_duplicates()
            .head(10)
            .apply(lambda row: f"{row[sample_col]} / {row[gene_col]}", axis=1)
            .tolist()
        )
        raise ValueError(
            "Duplicate sample-gene rows would make the plot ambiguous. "
            "Collapse them before calling this function. Examples: "
            + "; ".join(duplicate_pairs)
        )

    for column in (h0_col, h1_col, h2_col):
        df[column] = pd.to_numeric(df[column], errors="coerce")

    denominator = df[h1_col] + df[h2_col]
    total_reads = denominator + df[h0_col]
    numerator = (df[h1_col] - df[h2_col]).abs()
    phased_proportion = denominator.div(total_reads).where(total_reads > 0)

    valid_skew = (
        denominator.ge(min_total)
        & phased_proportion.ge(phased_prop_min)
        & denominator.gt(0)
    )
    df["haplotype_skew"] = numerator.div(denominator).where(valid_skew)
    df["phased_proportion"] = phased_proportion

    tidy = df.rename(
        columns={sample_col: "Sample", gene_col: "Gene"}
    )[
        [
            "Sample",
            "Gene",
            h0_col,
            h1_col,
            h2_col,
            "haplotype_skew",
            "phased_proportion",
        ]
    ].copy()

    sample_rank = {sample: rank for rank, sample in enumerate(sample_order)}
    gene_rank = {gene: rank for rank, gene in enumerate(genes_of_interest)}
    tidy["_sample_rank"] = tidy["Sample"].map(sample_rank)
    tidy["_gene_rank"] = tidy["Gene"].map(gene_rank)
    tidy = (
        tidy.sort_values(["_sample_rank", "_gene_rank", "Sample", "Gene"])
        .drop(columns=["_sample_rank", "_gene_rank"])
        .reset_index(drop=True)
    )

    output_table = Path(output_table)
    output_figure = Path(output_figure)
    output_table.parent.mkdir(parents=True, exist_ok=True)
    output_figure.parent.mkdir(parents=True, exist_ok=True)
    tidy.to_csv(output_table, sep="\t", index=False)
    print(f"Wrote table: {output_table}")

    pivot_skew = tidy.pivot(index="Sample", columns="Gene", values="haplotype_skew")
    pivot_phase = tidy.pivot(index="Sample", columns="Gene", values="phased_proportion")

    plotted_sample_order = [sample for sample in sample_order if sample in pivot_skew.index]
    pivot_skew = pivot_skew.reindex(plotted_sample_order)
    pivot_phase = pivot_phase.reindex(plotted_sample_order)

    valid_phase = (
        tidy.loc[tidy["haplotype_skew"].notna(), "phased_proportion"]
        .dropna()
        .astype(float)
    )
    if valid_phase.empty:
        global_min, global_max = 0.0, 0.0
    else:
        global_min, global_max = valid_phase.min(), valid_phase.max()

    rng = np.random.default_rng(random_seed)
    figure, axis = plt.subplots(figsize=(12, 6), dpi=150)
    size_min, size_max = 36, 196
    x_positions = np.arange(len(pivot_skew.index))
    plotted_any = False

    for gene in genes_of_interest:
        if gene not in pivot_skew.columns or gene not in pivot_phase.columns:
            continue

        y_values = pivot_skew[gene]
        phase_values = pivot_phase[gene]
        mask = y_values.notna() & phase_values.notna()
        if not mask.any():
            continue

        plotted_any = True
        phase_subset = phase_values.loc[mask].astype(float).to_numpy()
        if global_min == global_max:
            sizes = np.full(mask.sum(), (size_min + size_max) / 2.0)
        else:
            sizes = np.interp(
                phase_subset,
                (global_min, global_max),
                (size_min, size_max),
            )

        jitter = rng.uniform(-0.2, 0.2, size=mask.sum())
        x_jittered = x_positions[mask.to_numpy()] + jitter
        axis.scatter(
            x_jittered,
            y_values.loc[mask].to_numpy(),
            s=sizes,
            label=gene,
            alpha=0.8,
            edgecolors="none",
            linewidth=0.4,
        )

    axis.set_ylim(0, 1)
    axis.set_ylabel(
        f"Haplotype skew  |H1 − H2| / (H1 + H2)\n"
        f"(only if H1 + H2 ≥ {min_total} and "
        f"(H1 + H2)/(H1 + H2 + H0) ≥ {phased_prop_min}; "
        f"point size ∝ phased proportion, globally scaled)"
    )
    axis.set_xlabel("Sample")
    axis.set_title(
        "Haplotype skew per gene (Noncyclo counts), "
        "bubble size = phased proportion (global scaling)"
    )
    axis.set_xticks(x_positions)
    axis.set_xticklabels(pivot_skew.index, rotation=60, ha="right")

    if plotted_any:
        axis.legend(
            title="Gene",
            ncols=3,
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            borderaxespad=0.0,
        )

    figure.tight_layout()
    figure.savefig(output_figure, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote figure: {output_figure}")

    return tidy