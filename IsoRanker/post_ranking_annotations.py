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

    csv_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".csv.gz") and keyword.lower() in f.lower()]

    if not csv_files:
        print(f"No matching files found for keyword '{keyword}'. Skipping...")
        return

    first_df = pd.read_csv(csv_files[0], compression="gzip")
    all_columns = list(first_df.columns)
    all_columns.append("Source_File")

    for file in csv_files[1:]:
        df = pd.read_csv(file, compression="gzip")
        new_columns = [col for col in df.columns if col not in all_columns]
        all_columns.extend(new_columns)

    merged_df = pd.concat(
        [pd.read_csv(f, compression="gzip").reindex(columns=all_columns).assign(Source_File=os.path.basename(f)) for f in csv_files],
        ignore_index=True
    )

    merged_df.to_csv(output_csv, index=False, compression="gzip")
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

    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Open VCF
    vcf = pysam.VariantFile(vcf_file)

    # Extract CSQ header fields
    csq_fields = vcf.header.info['CSQ'].description.split(": ")[1].split("|")

    variant_data = []

    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alts = record.alts
        qual = record.qual  # Variant quality score

        if qual is None or qual < 30:  # Filter out low-quality variants
            continue

        # Extract all INFO fields except CSQ
        variant_info = {key: record.info.get(key, None) for key in vcf.header.info.keys() if key != "CSQ"}

        # Extract CSQ field and create separate rows
        if "CSQ" in record.info:
            for csq_entry in record.info["CSQ"]:
                csq_values = csq_entry.split("|")
                csq_dict = dict(zip(csq_fields, csq_values))  # Convert CSQ into dictionary

                # Extract MAX_AF safely (handling scientific notation)
                max_af_str = csq_dict.get("MAX_AF", "1")  # Default to "1" if missing
                try:
                    max_af = float(max_af_str)  # Converts both decimal & scientific notation safely
                except ValueError:
                    max_af = 1.0  # If conversion fails, default to 1.0

                # Apply MAX_AF filter (removes common variants)
                if max_af >= 0.01:
                    continue  # Skip high-frequency variants

                # Extract SpliceAI fields and get the highest SpliceAI score
                spliceai_columns = [
                    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"
                ]

                spliceai_scores = []
                for col in spliceai_columns:
                    value = csq_dict.get(col, "0")  # Get value, default to "0" if missing
                    try:
                        spliceai_scores.append(float(value))  # Convert safely
                    except ValueError:
                        spliceai_scores.append(0)  # Default to 0 if conversion fails
                
                highest_spliceai_score = max(spliceai_scores)  # Get the highest score

                # Determine SpliceAI Score Category
                if highest_spliceai_score > 0.5:
                    spliceai_high_score_flag = "SpliceAI_High_Score"
                elif highest_spliceai_score > 0.2:
                    spliceai_high_score_flag = "SpliceAI_Moderate_Score"
                else:
                    spliceai_high_score_flag = ""  # Empty if = 0.2

                # Extract SpliceAI_pred_SYMBOL
                spliceai_pred_symbol = csq_dict.get("SpliceAI_pred_SYMBOL", "")

                # Initialize row with variant and CSQ data
                variant_entry = {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alts[0],
                    "QUAL": qual,
                    "MAX_AF": max_af,
                    "Highest_SpliceAI_Score": highest_spliceai_score,
                    "SpliceAI_High_Score_Flag": spliceai_high_score_flag,
                    "SpliceAI_pred_SYMBOL": spliceai_pred_symbol,
                    **variant_info,  # Add all INFO fields
                    **csq_dict,  # Add all CSQ fields
                }

                # Extract genotypes for each sample (handling phased/unphased)
                for sample in vcf.header.samples:
                    genotype_tuple = record.samples[sample]["GT"]  # Get GT as a tuple (e.g., (0,1))
                    is_phased = record.samples[sample].phased  # Check if phased

                    if genotype_tuple is not None:
                        separator = "|" if is_phased else "/"  # Use "|" for phased, "/" for unphased
                        genotype_str = separator.join(map(str, genotype_tuple))  # Convert to string
                    else:
                        genotype_str = "NA"  # Handle missing genotype

                    variant_entry[f"GT_{sample}"] = genotype_str  # Store GT per sample

                # Append the variant entry to the list
                variant_data.append(variant_entry)

    # Convert to DataFrame
    df = pd.DataFrame(variant_data)

    # Save df.csv
    df_csv_path = os.path.join(output_dir, f"{output_filename}_vcf_vep_filtered.csv")
    df.to_csv(df_csv_path, index=False)
    print(f"Saved: {df_csv_path}")

    # Process df to generate df_final
    required_columns = {"SYMBOL", "HGVSg", "MAX_AF", "CLIN_SIG", "Consequence"}
    genotype_columns = [col for col in df.columns if col.startswith("GT_")]  # Identify genotype columns

    # Initialize new structure
    haplotype_data = {}

    # Process each variant in df
    for _, row in df.iterrows():
        gene = row["SYMBOL"]
        hgvs = row["HGVSg"]
        max_af = row["MAX_AF"]
        clin_sig = row["CLIN_SIG"]
        consequence = row["Consequence"]
        highest_spliceai_score = row["Highest_SpliceAI_Score"]
        spliceai_high_score_flag = row["SpliceAI_High_Score_Flag"]
        spliceai_pred_symbol = row["SpliceAI_pred_SYMBOL"]

        # Format the variant string with SpliceAI information
        variant_info = f"({hgvs},{max_af},{clin_sig},{consequence},{highest_spliceai_score},{spliceai_high_score_flag},{spliceai_pred_symbol})"

        if gene not in haplotype_data:
            haplotype_data[gene] = {"Hap1": [], "Hap2": [], "Hap0": []}

        for sample in genotype_columns:
            genotype = row[sample]

            if genotype == "1|0":  # Heterozygous (haplotype 1)
                haplotype_data[gene]["Hap1"].append(variant_info)

            elif genotype == "0|1":  # Heterozygous (haplotype 2)
                haplotype_data[gene]["Hap2"].append(variant_info)

            elif genotype == "0/1":  # Unphased heterozygous
                haplotype_data[gene]["Hap0"].append(variant_info)

            elif genotype == "1/1":  # Homozygous ? Include in both Hap1 and Hap2
                haplotype_data[gene]["Hap1"].append(variant_info)
                haplotype_data[gene]["Hap2"].append(variant_info)

    # Convert to DataFrame
    final_data = [
        {"Gene": gene, 
         "Hap1": "; ".join(set(data["Hap1"])), 
         "Hap2": "; ".join(set(data["Hap2"])), 
         "Hap0": "; ".join(set(data["Hap0"]))}
        for gene, data in haplotype_data.items()
    ]

    df_final = pd.DataFrame(final_data)

    # Save df_final.csv
    df_final_csv_path = os.path.join(output_dir, f"{output_filename}_gene_haplotype_split.csv")
    df_final.to_csv(df_final_csv_path, index=False)
    print(f"Saved: {df_final_csv_path}")

    return df, df_final


# I think this is deprecated
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
    Processes phenotype, gene, and disease data to compute HPO similarity scores between probands and OMIM diseases.

    Parameters:
    hpo_file (str): Path to the phenotype.hpoa file.
    genemap_file (str): Path to the genemap2 file.
    probands_file (str): Path to the probands file containing 'Sample' and 'Phenotype' columns.
    output_prefix (str): Prefix for output files. Default is 'all_proband_by_omim_comparison_ic'.

    Outputs:
    - all_proband_by_omim_comparison_ic.tsv
    - all_proband_by_omim_comparison_ic_longFormat.tsv
    """

    # Load HPO annotations
    hpo = pd.read_csv(hpo_file, sep='\t', comment='#')
    hpo = hpo[hpo['aspect'].isin(['P', 'H'])][['database_id', 'hpo_id']].drop_duplicates()
    hpo_dict = hpo.groupby('database_id')['hpo_id'].apply(lambda x: x.unique().tolist()).to_dict()

    # Load and process probands data
    probands = pd.read_csv(probands_file)[['Sample', 'Phenotype']].drop_duplicates()
    probands = probands[probands['Phenotype'].notnull()]

    Ontology()  # Load HPO Ontology

    def string_to_hpo(input_list, silent=False):
        """Convert phenotype terms to HPO terms, handling unknown terms gracefully."""
        output_list = []
        for term in input_list:
            try:
                hpo_term = str(Ontology.get_hpo_object(term.strip())).split(" | ")[0]
                output_list.append(hpo_term)
            except RuntimeError:
                if not silent:
                    print(f"Warning: Proband HPO term not found for '{term.strip()}'", flush=True)
        return output_list if output_list else np.nan  # Return NaN if no terms were found

    probands['hpo_terms'] = probands['Phenotype'].apply(lambda x: string_to_hpo(x.strip().split(",")))
    probands['hpo_set'] = probands['hpo_terms'].apply(HPOSet.from_queries)

    # Load and process genemap file
    genemap = pd.read_csv(genemap_file, sep='\t', skiprows=3)  # Read the file, skipping the first 3 rows

    # Extract phenotype details from the 'Phenotypes' column
    def extract_phenotype_details(phenotype):
        """Extract phenotype details using regex."""
        phenotype_pattern = r'^(?:\{|\[)?(.*?)(?:\}|\])?, (\d{6}) \((\d+)\)(?:, )?(.*)?$'
        match = re.match(phenotype_pattern, phenotype)
        
        if match:
            phenotype_name, omim_num, code, inheritance = match.groups()
        else:
            phenotype_name, omim_num, code, inheritance = np.nan, np.nan, np.nan, np.nan

        if inheritance == "":
            inheritance = np.nan
        return phenotype_name, omim_num, code, inheritance

    new_rows = []
    for _, row in genemap.iterrows():
        phenotypes = row['Phenotypes']

        if pd.isna(phenotypes):
            new_row = row.to_dict()
            new_row.update({
                'PhenotypesExtracted_PhenotypeName': np.nan,
                'PhenotypesExtracted_OMIMnum': np.nan,
                'PhenotypesExtracted_code': np.nan,
                'PhenotypesExtracted_inheritance': np.nan
            })
            new_rows.append(new_row)
        else:
            for phenotype in phenotypes.split(';'):
                phenotype = phenotype.strip()
                phenotype_name, omim_num, code, inheritance = extract_phenotype_details(phenotype)
                new_row = row.to_dict()
                new_row.update({
                    'PhenotypesExtracted_PhenotypeName': phenotype_name,
                    'PhenotypesExtracted_OMIMnum': omim_num,
                    'PhenotypesExtracted_code': code,
                    'PhenotypesExtracted_inheritance': inheritance
                })
                new_rows.append(new_row)

    # Create a new DataFrame with extracted phenotype details
    genemap = pd.DataFrame(new_rows)

    # Keep only genes with known molecular basis (code 3)
    genemap = genemap[genemap["PhenotypesExtracted_code"] == "3"]
    genemap = genemap[genemap['Approved Gene Symbol'].notnull()]
    genemap = genemap[['Approved Gene Symbol', 'PhenotypesExtracted_OMIMnum', 'PhenotypesExtracted_code']].drop_duplicates()

    genemap['PhenotypesExtracted_OMIMnum'] = genemap['PhenotypesExtracted_OMIMnum'].apply(
        lambda x: f"OMIM:{int(x)}" if pd.notna(x) else x
    )
    genemap['HPO_terms'] = genemap['PhenotypesExtracted_OMIMnum'].map(hpo_dict).fillna(np.nan)

    def filter_valid_hpo_terms(hpo_terms):
        """Filters out invalid HPO terms before creating an HPOSet, logging invalid terms."""
        if not isinstance(hpo_terms, list) or len(hpo_terms) == 0:
            return np.nan  # Return NaN if empty or invalid
        
        valid_terms = []
        invalid_terms = []
    
        for term in hpo_terms:
            try:
                if Ontology.get_hpo_object(term):  # Ensure term exists
                    valid_terms.append(term)
            except RuntimeError:
                invalid_terms.append(term)  # Log invalid terms
    
        # Print invalid terms if any
        if invalid_terms:
            print(f"Warning: Invalid genemap HPO terms found and ignored: {invalid_terms}", flush=True)
    
        return HPOSet.from_queries(valid_terms) if valid_terms else np.nan  # Return NaN if no valid terms

    # Apply the function safely
    genemap['HPO_set'] = genemap['HPO_terms'].apply(filter_valid_hpo_terms)

    # Convert HPOSet to a string for deduplication
    genemap['HPO_set_str'] = genemap['HPO_set'].apply(lambda x: str(x) if isinstance(x, HPOSet) else np.nan)

    # Drop duplicates using both OMIM number and HPO_set_str
    genemap_by_omim = genemap[['PhenotypesExtracted_OMIMnum', 'HPO_set', 'HPO_set_str']].drop_duplicates(subset=['PhenotypesExtracted_OMIMnum', 'HPO_set_str'])


    # Compute similarity scores
    all_comparisons = pd.DataFrame(index=genemap_by_omim['PhenotypesExtracted_OMIMnum'], columns=probands['Sample'])

    for o in genemap_by_omim.itertuples():
        for p in probands.itertuples():
            p_set = p.hpo_set
            o_set = o.HPO_set
            similarity = p_set.similarity(o_set, method='ic') if pd.notna(o_set) else np.nan
            all_comparisons.at[o.PhenotypesExtracted_OMIMnum, p.Sample] = similarity

    all_comparisons_long = all_comparisons.reset_index().melt(
        id_vars=['PhenotypesExtracted_OMIMnum'], var_name='Sample', value_name='similarity'
    )

    # Convert OMIM-associated HPO terms to phenotype strings
    def hpo_codes_to_strings(hpo_terms):
        """
        Convert a list of HPO codes to their human-readable phenotype descriptions.
        
        Parameters:
        hpo_terms (list): List of HPO term codes (e.g., ['HP:0001250', 'HP:0004322']).
        
        Returns:
        list: List of human-readable HPO phenotype descriptions.
        """
        if not isinstance(hpo_terms, list):
            return np.nan  # Return NaN if not a list
        
        phenotype_strings = []
        for hpo_code in hpo_terms:
            try:
                hpo_object = Ontology.get_hpo_object(hpo_code.strip())
                phenotype_strings.append(hpo_object.name)  # Extract human-readable name
            except RuntimeError:
                continue  # Skip invalid terms
        
        return phenotype_strings if phenotype_strings else np.nan  # Return NaN if empty

    # **Add Gene Name, Sample HPO Terms, and OMIM HPO Terms to `all_comparisons_long`**
    genemap['HPO_terms_OMIM'] = genemap['HPO_terms'].apply(hpo_codes_to_strings)
    probands['HPO_terms_Sample'] = probands['hpo_terms'].apply(hpo_codes_to_strings)
    all_comparisons_long = all_comparisons_long.merge(
        genemap[['PhenotypesExtracted_OMIMnum', 'Approved Gene Symbol', 'HPO_terms_OMIM']],
        on='PhenotypesExtracted_OMIMnum',
        how='left'
    ).merge(
        probands[['Sample', 'HPO_terms_Sample']],
        on='Sample',
        how='left'
    )

    # Save results
    all_comparisons.to_csv(f"{output_prefix}.csv.gz", compression = "gzip")
    all_comparisons_long.to_csv(f"{output_prefix}_longFormat.csv.gz", index=False, compression = "gzip")

    return all_comparisons, all_comparisons_long