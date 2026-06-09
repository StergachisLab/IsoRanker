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


#Read only needed read_stats columns, not the full file.
#Read only needed haplotype columns, not full haplotype files.
#Avoid readlines(), which loads each whole haplotype file just to check if empty.
#Delete large intermediates after merge/write.
def update_files_with_haplotype_info(
    sample_info_with_haplotype_location,
    read_stats_path,
    output_dir,
    isoform_col=None,
):
    import os
    import gc
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)

    updated_read_stats_path = os.path.join(output_dir, "updated_read_stats.txt.gz")
    updated_sample_info_path = os.path.join(output_dir, "updated_sample_info.tsv.gz")

    sample_info = pd.read_csv(sample_info_with_haplotype_location, sep="\t", dtype=str)

    def normalize_haplotype(s):
        return (
            s.fillna("H0")
            .astype(str)
            .str.strip()
            .replace({
                "0": "H0", "1": "H1", "2": "H2",
                "0.0": "H0", "1.0": "H1", "2.0": "H2",
                "H0": "H0", "H1": "H1", "H2": "H2",
                "h0": "H0", "h1": "H1", "h2": "H2",
                "none": "H0", "None": "H0",
                "nan": "H0", "NaN": "H0",
                "": "H0",
            })
        )

    def find_isoform_col(columns):
        if isoform_col is not None:
            if isoform_col not in columns:
                raise ValueError(
                    f"isoform_col='{isoform_col}' not found in read_stats. "
                    f"Columns are: {list(columns)}"
                )
            return isoform_col

        for c in [
            "isoform", "Isoform",
            "pbid", "PBID",
            "associated_transcript",
            "transcript_id",
            "isoform_id",
        ]:
            if c in columns:
                return c

        # Fallback: old behavior-ish, assume second column is isoform assignment.
        return columns[1]

    def find_read_col(columns):
        for c in ["#readname", "ReadName", "readname", "id", "ID"]:
            if c in columns:
                return c
        raise ValueError(f"No read-name column found. Columns: {list(columns)}")

    def find_hap_col(columns):
        for c in ["haplotype", "Haplotype", "HP"]:
            if c in columns:
                return c
        raise ValueError(f"No haplotype column found. Columns: {list(columns)}")

    def write_updated_sample_info():
        updated_sample_info = []

        for _, row in sample_info.iterrows():
            row_dict = row.to_dict()
            sample = row_dict["sample"]
            haplotype = row_dict.get("haplotype", None)

            if pd.notna(haplotype) and str(haplotype).strip():
                for hap in ["H0", "H1", "H2"]:
                    new_row = row_dict.copy()
                    new_row["sample"] = f"{sample}{hap}"
                    new_row["haplotype"] = hap
                    updated_sample_info.append(new_row)
            else:
                new_row = row_dict.copy()
                new_row["sample"] = f"{sample}H0"
                new_row["haplotype"] = "H0"
                updated_sample_info.append(new_row)

        pd.DataFrame(updated_sample_info).to_csv(
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
    # Read only id + isoform from read_stats.
    # This saves a lot of RAM compared with loading all columns.
    # ------------------------------------------------------------
    read_stats_header = pd.read_csv(
        read_stats_path,
        sep="\t",
        nrows=0,
        compression="infer",
    ).columns.tolist()

    read_id_col = read_stats_header[0]
    iso_col = find_isoform_col(read_stats_header)

    print(f"Using read_stats read column: {read_id_col}", flush=True)
    print(f"Using read_stats isoform column: {iso_col}", flush=True)

    read_stats = pd.read_csv(
        read_stats_path,
        sep="\t",
        usecols=[read_id_col, iso_col],
        dtype=str,
        compression="infer",
    )

    read_stats = read_stats.rename(columns={read_id_col: "id", iso_col: "isoform"})
    read_stats["id"] = read_stats["id"].astype(str).str.strip()
    read_stats = read_stats[read_stats["id"].notna() & (read_stats["id"] != "")]

    print(f"Loaded compact read_stats rows: {len(read_stats):,}", flush=True)

    # ------------------------------------------------------------
    # Read only id + haplotype from each haplotype file.
    # ------------------------------------------------------------
    all_haplotypes = []

    for _, row in sample_info.iterrows():
        sample = row["sample"]
        haplotype_file = row.get("haplotype", None)

        print(f"Processing haplotype assignment file for {sample}", flush=True)

        if pd.isna(haplotype_file) or not str(haplotype_file).strip():
            continue

        haplotype_file = str(haplotype_file).strip()

        try:
            # Do not use readlines(); just read the header.
            hap_header = pd.read_csv(
                haplotype_file,
                sep="\t",
                nrows=0,
                dtype=str,
            ).columns.tolist()

            hap_read_col = find_read_col(hap_header)
            hap_col = find_hap_col(hap_header)

            for chunk in pd.read_csv(
                haplotype_file,
                sep="\t",
                usecols=[hap_read_col, hap_col],
                chunksize=1_000_000,
                dtype=str,
            ):
                if chunk.empty:
                    continue

                chunk = chunk.rename(columns={hap_read_col: "id", hap_col: "haplotype"})
                chunk["id"] = chunk["id"].astype(str).str.strip()
                chunk = chunk[chunk["id"].notna() & (chunk["id"] != "")]
                chunk["haplotype"] = normalize_haplotype(chunk["haplotype"])

                all_haplotypes.append(chunk[["id", "haplotype"]])

        except FileNotFoundError:
            print(f"Warning: File {haplotype_file} not found. Skipping.", flush=True)
        except pd.errors.EmptyDataError:
            print(f"Warning: {haplotype_file} is empty. Skipping.", flush=True)
        except ValueError as e:
            print(f"Warning: {e}. Skipping {haplotype_file}.", flush=True)

    if all_haplotypes:
        haplotype_df = pd.concat(all_haplotypes, ignore_index=True)
        del all_haplotypes
    else:
        haplotype_df = pd.DataFrame(columns=["id", "haplotype"])

    print(f"Loaded haplotype rows: {len(haplotype_df):,}", flush=True)

    # Optional but often helpful if duplicate read IDs exist.
    # Keeps old dict-like behavior less likely to duplicate rows.
    #haplotype_df = haplotype_df.drop_duplicates(subset=["id"], keep="last")

    #print(f"Haplotype rows after duplicate removal: {len(haplotype_df):,}", flush=True)

    # ------------------------------------------------------------
    # Merge compact read_stats with haplotypes.
    # Use left join to preserve old behavior: unassigned reads become H0.
    # ------------------------------------------------------------
    read_stats = read_stats.merge(haplotype_df, on="id", how="left")

    del haplotype_df
    gc.collect()

    read_stats["haplotype"] = normalize_haplotype(read_stats["haplotype"])

    split_id = read_stats["id"].str.split("_m", n=1, expand=True)

    read_stats["id"] = (
        split_id[0].to_numpy()
        + read_stats["haplotype"].to_numpy()
        + "_m"
        + split_id[1].to_numpy()
    )

    read_stats = read_stats[["id", "isoform"]]

    read_stats.to_csv(
        updated_read_stats_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(f"Updated read_stats saved to {updated_read_stats_path}", flush=True)

    del read_stats
    gc.collect()

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
