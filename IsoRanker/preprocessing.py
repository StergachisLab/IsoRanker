import os
import pandas as pd
import gzip


# def update_files_with_haplotype_info(sample_info_with_haplotype_location, read_stats_path, output_dir):
#     """
#     Updates read statistics and sample metadata with haplotype information.

#     This function performs the following steps:
#     - Loads `sample_info_with_haplotype_location`, which contains paths to haplotype assignment files.
#     - Loads a read-level summary table (`read_stats.txt`) with IDs for each read.
#     - For each sample, appends haplotype assignments (H0/H1/H2) to the read ID in `read_stats`.
#     - Generates a new version of `sample_info` where each sample is duplicated per haplotype.
#     - Writes two outputs:
#         - `updated_read_stats.txt.gz`: updated read stats with haplotype-annotated IDs.
#         - `updated_sample_info.tsv.gz`: expanded sample metadata for each haplotype.

#     Parameters:
#     - sample_info_with_haplotype_location (str): Path to a TSV file containing columns `sample` and `haplotype`, where `haplotype` points to a file with read-level haplotype assignments.
#     - read_stats_path (str): Path to the original read stats file (`read_stats.txt`), with a column `id`.
#     - output_dir (str): Directory where updated output files will be written.

#     Returns:
#     - None: Writes updated read stats and sample info files to the specified output directory.
#     """
    
#     import os
#     import pandas as pd

#     os.makedirs(output_dir, exist_ok=True)

#     updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
#     updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")

#     sample_info = pd.read_csv(sample_info_with_haplotype_location, sep="\t")

#     def write_updated_sample_info():
#         updated_sample_info = []

#         for _, row in sample_info.iterrows():
#             row_dict = row.to_dict()

#             sample = row_dict["sample"]
#             haplotype = row_dict.get("haplotype", None)

#             if pd.notna(haplotype) and str(haplotype).strip():
#                 for hap in ["H0", "H1", "H2"]:
#                     new_row = row_dict.copy()
#                     new_row["sample"] = f"{sample}{hap}"
#                     new_row["haplotype"] = hap
#                     updated_sample_info.append(new_row)
#             else:
#                 new_row = row_dict.copy()
#                 new_row["sample"] = f"{sample}H0"
#                 new_row["haplotype"] = "H0"
#                 updated_sample_info.append(new_row)

#         updated_sample_info_df = pd.DataFrame(updated_sample_info)

#         updated_sample_info_df.to_csv(
#             updated_sample_info_path,
#             sep="\t",
#             index=False,
#             compression="gzip"
#         )

#         print(
#             f"Updated sample info saved to {updated_sample_info_path}",
#             flush=True
#         )


#     if os.path.exists(updated_read_stats_path):
#         print(
#             f"{updated_read_stats_path} already exists. "
#             "Skipping read_stats generation and creating updated_sample_info only.",
#             flush=True
#         )
#         write_updated_sample_info()
#         return

#     read_stats = pd.read_csv(
#         read_stats_path,
#         sep="\t",
#         dtype={0: str},
#         compression="infer"
#     )

#     read_stats.rename(columns={read_stats.columns[0]: "id"}, inplace=True)

#     all_haplotypes = []

#     for _, row in sample_info.iterrows():
#         sample = row["sample"]
#         haplotype_file = row["haplotype"]

#         print(f"Processing haplotype assignment file for {sample}", flush=True)

#         if pd.notna(haplotype_file) and str(haplotype_file).strip():
#             try:
#                 with open(haplotype_file, "r") as f:
#                     lines = f.readlines()
#                     if len(lines) <= 1:
#                         print(
#                             f"Warning: {haplotype_file} contains only headers or is empty. Skipping.",
#                             flush=True
#                         )
#                         continue

#                 for chunk in pd.read_csv(
#                     haplotype_file,
#                     sep="\t",
#                     chunksize=500000,
#                     dtype=str
#                 ):
#                     if chunk.empty:
#                         continue

#                     column_map = {}

#                     if "#readname" in chunk.columns:
#                         column_map["#readname"] = "id"
#                     elif "ReadName" in chunk.columns:
#                         column_map["ReadName"] = "id"
#                     else:
#                         raise ValueError(
#                             f"{haplotype_file} must contain either '#readname' or 'ReadName'"
#                         )

#                     if "haplotype" in chunk.columns:
#                         column_map["haplotype"] = "haplotype"
#                     elif "Haplotype" in chunk.columns:
#                         column_map["Haplotype"] = "haplotype"
#                     else:
#                         raise ValueError(
#                             f"{haplotype_file} must contain either 'haplotype' or 'Haplotype'"
#                         )

#                     chunk = chunk.rename(columns=column_map)[["id", "haplotype"]]

#                     chunk["haplotype"] = (
#                         chunk["haplotype"]
#                         .fillna("H0")
#                         .astype(str)
#                         .str.strip()
#                         .replace({
#                             "0": "H0",
#                             "1": "H1",
#                             "2": "H2",
#                             "H0": "H0",
#                             "H1": "H1",
#                             "H2": "H2",
#                             "none": "H0",
#                             "None": "H0",
#                             "nan": "H0",
#                             "": "H0",
#                         })
#                     )

#                     all_haplotypes.append(chunk)

#             except FileNotFoundError:
#                 print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
#             except pd.errors.EmptyDataError:
#                 print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)

#     if all_haplotypes:
#         haplotype_df = pd.concat(all_haplotypes, ignore_index=True)
#     else:
#         haplotype_df = pd.DataFrame(columns=["id", "haplotype"])

#     read_stats = read_stats.merge(haplotype_df, on="id", how="left")

#     read_stats["haplotype"] = (
#         read_stats["haplotype"]
#         .fillna("H0")
#         .astype(str)
#         .str.strip()
#         .replace({
#             "0": "H0",
#             "1": "H1",
#             "2": "H2",
#             "H0": "H0",
#             "H1": "H1",
#             "H2": "H2",
#             "none": "H0",
#             "None": "H0",
#             "nan": "H0",
#             "": "H0",
#         })
#     )

#     split_id = read_stats["id"].str.split("_m", n=1)

#     read_stats["id"] = (
#         split_id.str[0]
#         + read_stats["haplotype"]
#         + "_m"
#         + split_id.str[1]
#     )

#     read_stats.drop(columns=["haplotype"], inplace=True)

#     read_stats.to_csv(
#         updated_read_stats_path,
#         sep="\t",
#         index=False,
#         compression="gzip"
#     )

#     print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

#     write_updated_sample_info()




# New version that should take less RAM. Hasn't been tested yet though.
def update_files_with_haplotype_info(sample_info_path, read_stats_path, output_dir):
    import os
    import gzip
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)

    updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
    updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")

    sample_info = pd.read_csv(sample_info_path, sep="\t", dtype=str)

    def normalize_haplotype(s):
        return (
            s.fillna("H0")
            .astype(str)
            .str.strip()
            .replace({
                "0": "H0",
                "1": "H1",
                "2": "H2",
                "0.0": "H0",
                "1.0": "H1",
                "2.0": "H2",
                "H0": "H0",
                "H1": "H1",
                "H2": "H2",
                "h0": "H0",
                "h1": "H1",
                "h2": "H2",
                "none": "H0",
                "None": "H0",
                "nan": "H0",
                "NaN": "H0",
                "": "H0",
            })
        )

    def write_updated_sample_info():
        updated_rows = []

        for _, row in sample_info.iterrows():
            row_dict = row.to_dict()
            sample = row_dict["sample"]
            haplotype_file = row_dict.get("haplotype", None)

            if pd.notna(haplotype_file) and str(haplotype_file).strip():
                for hap in ["H0", "H1", "H2"]:
                    new_row = row_dict.copy()
                    new_row["sample"] = f"{sample}{hap}"
                    new_row["haplotype"] = hap
                    updated_rows.append(new_row)
            else:
                new_row = row_dict.copy()
                new_row["sample"] = f"{sample}H0"
                new_row["haplotype"] = "H0"
                updated_rows.append(new_row)

        pd.DataFrame(updated_rows).to_csv(
            updated_sample_info_path,
            sep="\t",
            index=False,
            compression="gzip",
        )

        print(f"Updated sample info saved to {updated_sample_info_path}", flush=True)

    if os.path.exists(updated_read_stats_path):
        print(
            f"{updated_read_stats_path} already exists. "
            "Skipping read_stats generation and creating updated_sample_info only.",
            flush=True,
        )
        write_updated_sample_info()
        return

    # ------------------------------------------------------------
    # Build read_id -> haplotype map in chunks
    # Flexible input columns:
    #   read name: #readname, ReadName, readname, id
    #   haplotype: haplotype, Haplotype, HP
    # ------------------------------------------------------------
    hap_map = {}

    for _, row in sample_info.iterrows():
        sample = row["sample"]
        haplotype_file = row.get("haplotype", None)

        print(f"Processing haplotype assignment file for {sample}", flush=True)

        if pd.isna(haplotype_file) or not str(haplotype_file).strip():
            print(f"No haplotype file for {sample}; reads will default to H0.", flush=True)
            continue

        haplotype_file = str(haplotype_file).strip()

        try:
            header = pd.read_csv(
                haplotype_file,
                sep="\t",
                nrows=0,
                dtype=str,
            ).columns.tolist()

            if "#readname" in header:
                read_col = "#readname"
            elif "ReadName" in header:
                read_col = "ReadName"
            elif "readname" in header:
                read_col = "readname"
            elif "id" in header:
                read_col = "id"
            else:
                raise ValueError(
                    f"{haplotype_file} must contain one of: "
                    "'#readname', 'ReadName', 'readname', or 'id'. "
                    f"Found columns: {header}"
                )

            if "haplotype" in header:
                hap_col = "haplotype"
            elif "Haplotype" in header:
                hap_col = "Haplotype"
            elif "HP" in header:
                hap_col = "HP"
            else:
                raise ValueError(
                    f"{haplotype_file} must contain one of: "
                    "'haplotype', 'Haplotype', or 'HP'. "
                    f"Found columns: {header}"
                )

            for chunk in pd.read_csv(
                haplotype_file,
                sep="\t",
                usecols=[read_col, hap_col],
                dtype=str,
                chunksize=2_000_000,
            ):
                if chunk.empty:
                    continue

                chunk = chunk.rename(columns={read_col: "id", hap_col: "haplotype"})
                chunk["id"] = chunk["id"].astype(str).str.strip()
                chunk["haplotype"] = normalize_haplotype(chunk["haplotype"])

                # Remove rows with missing/empty read IDs.
                chunk = chunk[chunk["id"].notna() & (chunk["id"] != "")]

                hap_map.update(zip(chunk["id"], chunk["haplotype"]))

        except FileNotFoundError:
            print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
        except pd.errors.EmptyDataError:
            print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)
        except ValueError as e:
            print(f"Warning: {e}. Skipping.", flush=True)

    print(f"Total haplotyped reads loaded: {len(hap_map):,}", flush=True)

    # ------------------------------------------------------------
    # Stream read_stats instead of loading full file into RAM
    # ------------------------------------------------------------
    first_chunk = True
    total_reads = 0
    total_h1_h2 = 0

    with gzip.open(updated_read_stats_path, "wt") as out:
        for read_chunk in pd.read_csv(
            read_stats_path,
            sep="\t",
            dtype=str,
            compression="infer",
            chunksize=2_000_000,
        ):
            if read_chunk.empty:
                continue

            # Original code treats the first column as the read ID.
            read_chunk.rename(columns={read_chunk.columns[0]: "id"}, inplace=True)

            read_id = read_chunk["id"].astype(str)

            hap = read_id.map(hap_map).fillna("H0")
            hap = normalize_haplotype(hap)

            split_id = read_id.str.split("_m", n=1)

            read_chunk["id"] = (
                split_id.str[0]
                + hap
                + "_m"
                + split_id.str[1]
            )

            total_reads += len(read_chunk)
            total_h1_h2 += hap.isin(["H1", "H2"]).sum()

            read_chunk.to_csv(
                out,
                sep="\t",
                index=False,
                header=first_chunk,
            )

            first_chunk = False

    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)
    print(f"Total reads processed: {total_reads:,}", flush=True)
    print(f"Reads assigned H1/H2: {total_h1_h2:,}", flush=True)

    write_updated_sample_info()


def filter_based_on_counts(df, count_threshold=10, group_col='Isoform'):
    """
    Filter isoforms based on count thresholds.
    
    Parameters:
    - df (pd.DataFrame): The input long-format DataFrame.
    - count_threshold (int): The threshold for filtering counts.
    - group_col (str): The column to group by (e.g., 'Isoform', 'Isoform_PBid').
    
    Returns:
    - pd.DataFrame: Filtered DataFrame with only the groups meeting the count threshold.
    """
    # Determine isoforms/groups to keep based on the threshold
    isoforms_to_keep = df.groupby(group_col).apply(
        lambda group: any(group['cyclo_count'] >= count_threshold) or any(group['noncyclo_count'] >= count_threshold),
        include_groups=False
    )
    isoforms_to_keep = isoforms_to_keep[isoforms_to_keep].index.tolist()

    # Filter the DataFrame to include only the isoforms/groups to keep
    return df[df[group_col].isin(isoforms_to_keep)]
