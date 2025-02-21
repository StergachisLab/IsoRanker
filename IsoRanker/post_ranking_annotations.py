import os
import pysam
import pandas as pd
import numpy as np
import re
from pyhpo import stats, Ontology, HPOSet


def merge_csvs_by_keyword(directory, keyword, output_csv):
    """
    Merges CSV files in the specified directory that contain a given keyword in their filename.
    The function preserves column order from the first file and includes all columns across files.
    Adds a column `Source_File` to track the original file for each row.

    Parameters:
    - directory (str): Path to the directory containing CSV files.
    - keyword (str): Keyword to match files (e.g., "gene" or "isoform").
    - output_csv (str): Output file path for the merged CSV.

    Returns:
    - None: Saves the merged CSV to the specified output path.
    """

    csv_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".csv") and keyword.lower() in f.lower()]

    if not csv_files:
        print(f"No matching files found for keyword '{keyword}'. Skipping...")
        return

    first_df = pd.read_csv(csv_files[0])
    all_columns = list(first_df.columns)
    all_columns.append("Source_File")

    for file in csv_files[1:]:
        df = pd.read_csv(file)
        new_columns = [col for col in df.columns if col not in all_columns]
        all_columns.extend(new_columns)

    merged_df = pd.concat(
        [pd.read_csv(f).reindex(columns=all_columns).assign(Source_File=os.path.basename(f)) for f in csv_files],
        ignore_index=True
    )

    merged_df.to_csv(output_csv, index=False)
    print(f"Merged {len(csv_files)} '{keyword}' files into: {output_csv}")


def process_vep_vcf(vcf_file, output_dir, output_filename):
    """
    Processes a VEP-annotated VCF file to filter variants, extract genotypes,
    and generate structured summary tables including the highest SpliceAI score.

    Parameters:
    - vcf_file (str): Path to the input VCF file.
    - output_dir (str): Directory where output CSV files will be saved.
    - output_filename (str): Base name for the output files (without extension).

    Returns:
    - df (pd.DataFrame): DataFrame containing all filtered variants.
    - df_final (pd.DataFrame): DataFrame summarizing variants grouped by gene and haplotype.
    """

    os.makedirs(output_dir, exist_ok=True)
    vcf = pysam.VariantFile(vcf_file)
    csq_fields = vcf.header.info['CSQ'].description.split(": ")[1].split("|")

    variant_data = []

    for record in vcf:
        if record.qual is None or record.qual < 30:
            continue

        variant_info = {key: record.info.get(key, None) for key in vcf.header.info.keys() if key != "CSQ"}

        if "CSQ" in record.info:
            for csq_entry in record.info["CSQ"]:
                csq_values = csq_entry.split("|")
                csq_dict = dict(zip(csq_fields, csq_values))

                max_af = float(csq_dict.get("MAX_AF", "1"))
                if max_af >= 0.01:
                    continue

                spliceai_scores = [float(csq_dict.get(col, "0")) for col in 
                                   ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]]
                highest_spliceai_score = max(spliceai_scores)

                spliceai_high_score_flag = (
                    "SpliceAI_High_Score" if highest_spliceai_score > 0.5 else
                    "SpliceAI_Moderate_Score" if highest_spliceai_score > 0.2 else ""
                )

                variant_entry = {
                    "CHROM": record.chrom,
                    "POS": record.pos,
                    "REF": record.ref,
                    "ALT": record.alts[0],
                    "QUAL": record.qual,
                    "MAX_AF": max_af,
                    "Highest_SpliceAI_Score": highest_spliceai_score,
                    "SpliceAI_High_Score_Flag": spliceai_high_score_flag,
                    **variant_info,
                    **csq_dict
                }

                for sample in vcf.header.samples:
                    genotype_tuple = record.samples[sample]["GT"]
                    genotype_str = "|".join(map(str, genotype_tuple)) if record.samples[sample].phased else "/".join(map(str, genotype_tuple))
                    variant_entry[f"GT_{sample}"] = genotype_str

                variant_data.append(variant_entry)

    df = pd.DataFrame(variant_data)
    df.to_csv(os.path.join(output_dir, f"{output_filename}_vcf_vep_filtered.csv"), index=False)

    return df


def merge_haplotype_data(sample_gene_file, haplotype_patient_file, output_file):
    """
    Merges haplotype data into a master file containing Sample and associated_gene.

    Parameters:
    - sample_gene_file (str): Path to the CSV file containing 'Sample' and 'associated_gene' columns.
    - haplotype_patient_file (str): Path to the CSV file containing 'Gene', 'Hap1', 'Hap2', 'Hap0', and 'Patient'.
    - output_file (str): Path to save the merged CSV.

    Returns:
    - None: Saves the merged CSV to the specified output file.
    """

    sample_gene_df = pd.read_csv(sample_gene_file)
    haplotype_df = pd.read_csv(haplotype_patient_file)

    merged_df = sample_gene_df.merge(
        haplotype_df,
        left_on=["Sample", "associated_gene"],
        right_on=["Patient", "Gene"],
        how="left"
    ).drop(columns=["Patient", "Gene"], errors="ignore")

    merged_df.to_csv(output_file, index=False)
    print(f"Merged data saved to: {output_file}")


def process_phenotype_data(hpo_file, genemap_file, probands_file, output_prefix="all_proband_by_omim_comparison_ic"):
    """
    Computes HPO similarity scores between probands and OMIM diseases.

    Parameters:
    - hpo_file (str): Path to the phenotype.hpoa file.
    - genemap_file (str): Path to the genemap2 file.
    - probands_file (str): Path to the probands file containing 'Sample' and 'Phenotype' columns.
    - output_prefix (str): Prefix for output files. Default is 'all_proband_by_omim_comparison_ic'.

    Returns:
    - all_comparisons (pd.DataFrame)
    - all_comparisons_long (pd.DataFrame)
    """

    Ontology()

    probands = pd.read_csv(probands_file)[['Sample', 'Phenotype']].drop_duplicates()
    probands['hpo_terms'] = probands['Phenotype'].apply(lambda x: x.split(","))
    probands['hpo_set'] = probands['hpo_terms'].apply(HPOSet.from_queries)

    genemap = pd.read_csv(genemap_file, sep='\t', skiprows=3)
    genemap = genemap[genemap["PhenotypesExtracted_code"] == "3"]
    genemap['HPO_terms'] = genemap['PhenotypesExtracted_OMIMnum'].map(hpo_file)
    genemap['HPO_set'] = genemap['HPO_terms'].apply(HPOSet.from_queries)

    all_comparisons = pd.DataFrame(index=genemap['PhenotypesExtracted_OMIMnum'], columns=probands['Sample'])

    for o in genemap.itertuples():
        for p in probands.itertuples():
            similarity = p.hpo_set.similarity(o.HPO_set, method='ic') if pd.notna(o.HPO_set) else np.nan
            all_comparisons.at[o.PhenotypesExtracted_OMIMnum, p.Sample] = similarity

    all_comparisons_long = all_comparisons.reset_index().melt(id_vars=['PhenotypesExtracted_OMIMnum'], var_name='Sample', value_name='similarity')
    
    all_comparisons.to_csv(f"{output_prefix}.tsv", sep='\t')
    all_comparisons_long.to_csv(f"{output_prefix}_longFormat.tsv", sep='\t', index=False)

    return all_comparisons, all_comparisons_long
