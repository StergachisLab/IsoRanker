import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def replace_original_with_corrected_tpm_and_scale_counts(
    long_df,
    cyclo_tpm_col="Cyclo_TPM",
    noncyclo_tpm_col="Noncyclo_TPM",
    cyclo_tpm_corrected_col="Cyclo_TPM_N50_corrected",
    noncyclo_tpm_corrected_col="Noncyclo_TPM_N50_corrected",
):
    """
    Replace original TPM and count columns with corrected values.

    Original TPM/count columns are preserved using an '_original' suffix.
    Count columns are scaled by the same correction factor used for TPM.
    """

    long_df = long_df.copy()

    cyclo_count_cols = [
        "H0_cyclo_count",
        "H1_cyclo_count",
        "H2_cyclo_count",
        "cyclo_count",
    ]

    noncyclo_count_cols = [
        "H0_noncyclo_count",
        "H1_noncyclo_count",
        "H2_noncyclo_count",
        "noncyclo_count",
    ]

    # Preserve original TPMs
    long_df[cyclo_tpm_col + "_original"] = long_df[cyclo_tpm_col]
    long_df[noncyclo_tpm_col + "_original"] = long_df[noncyclo_tpm_col]

    # Convert TPM columns to numeric
    long_df[cyclo_tpm_col + "_original"] = pd.to_numeric(
        long_df[cyclo_tpm_col + "_original"], errors="coerce"
    ).fillna(0)

    long_df[noncyclo_tpm_col + "_original"] = pd.to_numeric(
        long_df[noncyclo_tpm_col + "_original"], errors="coerce"
    ).fillna(0)

    long_df[cyclo_tpm_corrected_col] = pd.to_numeric(
        long_df[cyclo_tpm_corrected_col], errors="coerce"
    )

    long_df[noncyclo_tpm_corrected_col] = pd.to_numeric(
        long_df[noncyclo_tpm_corrected_col], errors="coerce"
    )

    # Scaling factors
    long_df["Cyclo_scaling_factor"] = np.where(
        long_df[cyclo_tpm_col + "_original"] > 0,
        long_df[cyclo_tpm_corrected_col] / long_df[cyclo_tpm_col + "_original"],
        1.0,
    )

    long_df["Noncyclo_scaling_factor"] = np.where(
        long_df[noncyclo_tpm_col + "_original"] > 0,
        long_df[noncyclo_tpm_corrected_col] / long_df[noncyclo_tpm_col + "_original"],
        1.0,
    )

    long_df["Cyclo_scaling_factor"] = (
        long_df["Cyclo_scaling_factor"]
        .replace([np.inf, -np.inf], 1.0)
        .fillna(1.0)
    )

    long_df["Noncyclo_scaling_factor"] = (
        long_df["Noncyclo_scaling_factor"]
        .replace([np.inf, -np.inf], 1.0)
        .fillna(1.0)
    )

    # Replace TPM columns with corrected TPMs
    long_df[cyclo_tpm_col] = long_df[cyclo_tpm_corrected_col]
    long_df[noncyclo_tpm_col] = long_df[noncyclo_tpm_corrected_col]

    # Scale cyclo count columns
    for col in cyclo_count_cols:
        if col in long_df.columns:
            long_df[col + "_original"] = long_df[col]
            long_df[col] = (
                pd.to_numeric(long_df[col], errors="coerce").fillna(0)
                * long_df["Cyclo_scaling_factor"]
            )

    # Scale noncyclo count columns
    for col in noncyclo_count_cols:
        if col in long_df.columns:
            long_df[col + "_original"] = long_df[col]
            long_df[col] = (
                pd.to_numeric(long_df[col], errors="coerce").fillna(0)
                * long_df["Noncyclo_scaling_factor"]
            )

    return long_df


def correct_tpm_by_n50_with_pca(
    long_df,
    n50_df,
    grouping_col="Isoform",
    sample_col="Sample",
    cyclo_tpm_col="Cyclo_TPM",
    noncyclo_tpm_col="Noncyclo_TPM",
    n50_sample_col="Individual",
    n50_condition_col="Sample type",
    n50_value_col="Read length N50",
    tpm_filter=10,
    output_dir=None,
    output_suffix="_N50_corrected",
    cmap="viridis",
    make_plots=True,
    return_matrices=False,
):
    """
    Correct TPM values for read length N50, replace original TPM/count columns
    with corrected values, preserve originals with '_original', and run PCA
    before and after correction.

    This function is designed for IsoRanker-style long-format dataframes where
    each row corresponds to one sample-feature pair, and Cyclo/Noncyclo values
    are stored in separate columns.

    The correction model is:

        log2(TPM + 1) ~ Condition + Read_length_N50

    The function removes the linear N50-associated component while preserving
    the Cyclo vs Noncyclo condition effect.

    Parameters
    ----------
    long_df : pandas.DataFrame
        Long-format expression dataframe.

        Required columns by default:
            - "Sample"
            - grouping_col, usually "Isoform" or "associated_gene"
            - "Cyclo_TPM"
            - "Noncyclo_TPM"

        If count columns are present, they will also be scaled and replaced:
            - "H0_cyclo_count"
            - "H1_cyclo_count"
            - "H2_cyclo_count"
            - "cyclo_count"
            - "H0_noncyclo_count"
            - "H1_noncyclo_count"
            - "H2_noncyclo_count"
            - "noncyclo_count"

        Original TPM and count values are preserved with an "_original" suffix.

    n50_df : pandas.DataFrame
        Dataframe containing N50 annotations.

        Required columns by default:
            - "Individual"
            - "Sample type"
            - "Read length N50"

        There should be one row per sample-condition pair.

    grouping_col : str, default "Isoform"
        Feature column used to build the sample-by-feature matrix.

        Use:
            - "Isoform" for isoform-level correction/PCA
            - "associated_gene" for gene-level correction/PCA

    sample_col : str, default "Sample"
        Sample identifier column in `long_df`.

        These values must match `n50_df[n50_sample_col]`.

    cyclo_tpm_col : str, default "Cyclo_TPM"
        TPM column for cycloheximide-treated samples.

        This column is replaced with N50-corrected TPM-like values in the returned
        dataframe. Original values are saved as "Cyclo_TPM_original".

    noncyclo_tpm_col : str, default "Noncyclo_TPM"
        TPM column for untreated/noncyclo samples.

        This column is replaced with N50-corrected TPM-like values in the returned
        dataframe. Original values are saved as "Noncyclo_TPM_original".

    n50_sample_col : str, default "Individual"
        Sample identifier column in `n50_df`.

        These values must match `long_df[sample_col]`.

    n50_condition_col : str, default "Sample type"
        Condition column in `n50_df`.

        Accepted values are case-insensitive:
            - "cyclo"
            - "noncyclo"

    n50_value_col : str, default "Read length N50"
        Numeric column containing read length N50.

        This is treated as a continuous technical covariate and z-scored before
        regression.

    tpm_filter : float or int, default 10
        Feature filter used only for PCA.

        A feature is included in PCA if at least one sample-condition has:

            TPM > tpm_filter

        This does not affect which features are corrected. N50 correction is run
        on all nonzero features.

    output_dir : str or None, default None
        Directory where outputs are saved.

        If provided, the function writes:
            - corrected IsoRanker-ready dataframe
            - metadata with N50
            - PCA table before correction
            - PCA table after correction
            - PCA PDF before correction
            - PCA PDF after correction

        If None, no files are written.

    output_suffix : str, default "_N50_corrected"
        Suffix used for intermediate corrected TPM columns.

        With default TPM column names, intermediate columns are:
            - "Cyclo_TPM_N50_corrected"
            - "Noncyclo_TPM_N50_corrected"

        The final returned dataframe also replaces "Cyclo_TPM" and
        "Noncyclo_TPM" with these corrected values.

    cmap : str, default "viridis"
        Matplotlib colormap used to color PCA points by read length N50.

    make_plots : bool, default True
        Whether to display PCA plots with `plt.show()`.

        Set to False for non-interactive cluster jobs. PDFs are still saved if
        `output_dir` is provided.

    return_matrices : bool, default False
        Whether to return intermediate matrices.

        If True, results["matrices"] includes:
            - raw_tpm_matrix
            - log2_tpm_matrix
            - corrected_log2_tpm_matrix
            - corrected_tpm_like_matrix
            - log_matrix_for_pca
            - corrected_log_matrix_for_pca

    Returns
    -------
    results : dict
        Dictionary containing corrected dataframe, PCA outputs, metadata, and
        optionally matrices.

        Main keys:
            corrected_long_df : pandas.DataFrame
                IsoRanker-ready dataframe. Original TPM/count columns are replaced
                with corrected/scaled values. Original values are preserved with
                "_original" suffix.

            metadata : pandas.DataFrame
                Sample-condition metadata with N50 annotations.

            pca_before : pandas.DataFrame
                PCA coordinates before N50 correction.

            pca_after : pandas.DataFrame
                PCA coordinates after N50 correction.

            pca_model_before : sklearn.decomposition.PCA
                Fitted PCA object before correction.

            pca_model_after : sklearn.decomposition.PCA
                Fitted PCA object after correction.

            explained_before : numpy.ndarray
                Percent variance explained by PC1 and PC2 before correction.

            explained_after : numpy.ndarray
                Percent variance explained by PC1 and PC2 after correction.

            features_used_for_pca : list
                Features passing `tpm_filter` and used for PCA.

    Notes
    -----
    The returned `corrected_long_df` can be passed directly to IsoRanker because
    the standard TPM and count columns are replaced with corrected values.
    """

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    required_cols = {
        sample_col,
        grouping_col,
        cyclo_tpm_col,
        noncyclo_tpm_col,
    }

    missing = required_cols - set(long_df.columns)
    if missing:
        raise ValueError(f"long_df is missing required columns: {missing}")

    n50_required = {
        n50_sample_col,
        n50_condition_col,
        n50_value_col,
    }

    missing_n50 = n50_required - set(n50_df.columns)
    if missing_n50:
        raise ValueError(f"n50_df is missing required columns: {missing_n50}")

    # =========================
    # Clean N50 table
    # =========================

    n50 = n50_df.copy()
    n50.columns = n50.columns.str.strip()

    n50 = n50.rename(columns={
        n50_sample_col: "Base_Sample",
        n50_condition_col: "Condition",
        n50_value_col: "Read_length_N50",
    })

    n50["Condition"] = (
        n50["Condition"]
        .astype(str)
        .str.strip()
        .str.lower()
        .map({"cyclo": "Cyclo", "noncyclo": "Noncyclo"})
    )

    n50["Read_length_N50"] = pd.to_numeric(
        n50["Read_length_N50"],
        errors="coerce",
    )

    if n50["Condition"].isna().any():
        bad = n50[n50["Condition"].isna()]
        raise ValueError(
            "Some N50 condition values could not be mapped to Cyclo/Noncyclo:\n"
            f"{bad}"
        )

    # =========================
    # Build TPM matrix
    # =========================

    tmp = long_df.copy()

    tmp[cyclo_tpm_col] = pd.to_numeric(
        tmp[cyclo_tpm_col],
        errors="coerce",
    ).fillna(0)

    tmp[noncyclo_tpm_col] = pd.to_numeric(
        tmp[noncyclo_tpm_col],
        errors="coerce",
    ).fillna(0)

    grouped = (
        tmp.groupby([sample_col, grouping_col])[
            [cyclo_tpm_col, noncyclo_tpm_col]
        ]
        .sum()
        .reset_index()
    )

    cyclo = grouped[[sample_col, grouping_col, cyclo_tpm_col]].copy()
    cyclo.columns = ["Base_Sample", grouping_col, "TPM"]
    cyclo["Condition"] = "Cyclo"

    noncyclo = grouped[[sample_col, grouping_col, noncyclo_tpm_col]].copy()
    noncyclo.columns = ["Base_Sample", grouping_col, "TPM"]
    noncyclo["Condition"] = "Noncyclo"

    long_expr = pd.concat([cyclo, noncyclo], ignore_index=True)

    long_expr["Sample_Condition"] = (
        long_expr["Base_Sample"].astype(str)
        + "_"
        + long_expr["Condition"]
    )

    raw_matrix = long_expr.pivot_table(
        index="Sample_Condition",
        columns=grouping_col,
        values="TPM",
        aggfunc="sum",
        fill_value=0,
    )

    raw_matrix = raw_matrix.loc[raw_matrix.any(axis=1)]
    raw_matrix = raw_matrix.loc[:, raw_matrix.any(axis=0)]

    metadata = (
        long_expr[["Sample_Condition", "Base_Sample", "Condition"]]
        .drop_duplicates()
        .set_index("Sample_Condition")
        .loc[raw_matrix.index]
        .reset_index()
    )

    metadata = metadata.merge(
        n50[["Base_Sample", "Condition", "Read_length_N50"]],
        on=["Base_Sample", "Condition"],
        how="left",
    )

    if metadata["Read_length_N50"].isna().any():
        missing = metadata.loc[
            metadata["Read_length_N50"].isna(),
            ["Sample_Condition", "Base_Sample", "Condition"],
        ]
        raise ValueError(
            "Some sample-condition rows are missing N50 annotations:\n"
            f"{missing}"
        )

    # =========================
    # N50 correction
    # =========================

    log_matrix = np.log2(raw_matrix + 1)

    meta = metadata.set_index("Sample_Condition").loc[log_matrix.index].copy()

    n50_sd = meta["Read_length_N50"].std()

    if n50_sd == 0 or pd.isna(n50_sd):
        raise ValueError("Read_length_N50 has zero or undefined variance.")

    meta["N50_z"] = (
        meta["Read_length_N50"] - meta["Read_length_N50"].mean()
    ) / n50_sd

    design_condition = pd.get_dummies(
        meta["Condition"],
        prefix="Condition",
        drop_first=True,
        dtype=float,
    )

    intercept = pd.DataFrame(
        {"Intercept": np.ones(log_matrix.shape[0])},
        index=log_matrix.index,
    )

    design_n50 = meta[["N50_z"]].astype(float)

    X_full = pd.concat(
        [intercept, design_condition, design_n50],
        axis=1,
    )

    X_keep = pd.concat(
        [intercept, design_condition],
        axis=1,
    )

    Y = log_matrix.values

    beta_full = np.linalg.pinv(X_full.values) @ Y
    beta_keep = np.linalg.pinv(X_keep.values) @ Y

    fitted_full = X_full.values @ beta_full
    fitted_keep = X_keep.values @ beta_keep

    residuals = Y - fitted_full
    corrected_log = fitted_keep + residuals

    corrected_log_matrix = pd.DataFrame(
        corrected_log,
        index=log_matrix.index,
        columns=log_matrix.columns,
    )

    corrected_tpm_like_matrix = 2 ** corrected_log_matrix

    # =========================
    # PCA before/after correction
    # =========================

    keep_features = raw_matrix.max(axis=0)
    keep_features = keep_features[keep_features > tpm_filter].index

    log_matrix_for_pca = log_matrix.loc[:, keep_features]
    corrected_log_matrix_for_pca = corrected_log_matrix.loc[:, keep_features]

    def run_pca(matrix):
        scaled = StandardScaler().fit_transform(matrix)

        pca = PCA(n_components=2)
        pcs = pca.fit_transform(scaled)

        explained = pca.explained_variance_ratio_ * 100

        pca_df = pd.DataFrame(
            pcs,
            columns=["PC1", "PC2"],
            index=matrix.index,
        ).reset_index(names="Sample_Condition")

        pca_df = pca_df.merge(
            metadata,
            on="Sample_Condition",
            how="left",
        )

        return pca_df, pca, explained

    pca_before, pca_model_before, explained_before = run_pca(log_matrix_for_pca)

    pca_after, pca_model_after, explained_after = run_pca(
        corrected_log_matrix_for_pca
    )

    def plot_pca_colored_by_n50(pca_df, explained, title, output_pdf=None):
        fig, ax = plt.subplots(figsize=(10, 8))

        markers = {
            "Cyclo": "o",
            "Noncyclo": "X",
        }

        vmin = pca_df["Read_length_N50"].min()
        vmax = pca_df["Read_length_N50"].max()

        scatter_for_colorbar = None

        for condition, marker in markers.items():
            subset = pca_df[pca_df["Condition"] == condition].copy()

            sc = ax.scatter(
                subset["PC1"],
                subset["PC2"],
                c=subset["Read_length_N50"],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                marker=marker,
                s=120,
                edgecolor="black",
                linewidth=0.5,
                label=condition,
            )

            scatter_for_colorbar = sc

        for _, row in pca_df.iterrows():
            ax.text(
                row["PC1"],
                row["PC2"],
                row["Base_Sample"],
                fontsize=6,
                ha="right",
                va="bottom",
            )

        cbar = fig.colorbar(scatter_for_colorbar, ax=ax)
        cbar.set_label("Read length N50")

        ax.legend(title="Condition")
        ax.set_xlabel(f"PC1 ({explained[0]:.2f}%)")
        ax.set_ylabel(f"PC2 ({explained[1]:.2f}%)")
        ax.set_title(title)

        plt.tight_layout()

        if output_pdf is not None:
            plt.savefig(output_pdf, bbox_inches="tight")

        if make_plots:
            plt.show()
        else:
            plt.close()

    if output_dir is not None:
        before_pdf = os.path.join(
            output_dir,
            f"pca_before_N50_correction_{grouping_col}_any_TPM_gt_{tpm_filter}.pdf",
        )

        after_pdf = os.path.join(
            output_dir,
            f"pca_after_N50_correction_{grouping_col}_any_TPM_gt_{tpm_filter}.pdf",
        )
    else:
        before_pdf = None
        after_pdf = None

    plot_pca_colored_by_n50(
        pca_before,
        explained_before,
        title=f"Before N50 correction: {grouping_col}, any TPM > {tpm_filter}",
        output_pdf=before_pdf,
    )

    plot_pca_colored_by_n50(
        pca_after,
        explained_after,
        title=f"After N50 correction: {grouping_col}, any TPM > {tpm_filter}",
        output_pdf=after_pdf,
    )

    # =========================
    # Convert corrected matrix back to long format
    # =========================

    corrected_long = (
        corrected_tpm_like_matrix
        .reset_index(names="Sample_Condition")
        .melt(
            id_vars="Sample_Condition",
            var_name=grouping_col,
            value_name="Corrected_TPM",
        )
    )

    corrected_long = corrected_long.merge(
        metadata[["Sample_Condition", "Base_Sample", "Condition"]],
        on="Sample_Condition",
        how="left",
    )

    corrected_wide = (
        corrected_long
        .pivot_table(
            index=["Base_Sample", grouping_col],
            columns="Condition",
            values="Corrected_TPM",
            aggfunc="first",
        )
        .reset_index()
    )

    corrected_wide.columns.name = None

    corrected_wide = corrected_wide.rename(columns={
        "Base_Sample": sample_col,
        "Cyclo": cyclo_tpm_col + output_suffix,
        "Noncyclo": noncyclo_tpm_col + output_suffix,
    })

    corrected_long_df = long_df.merge(
        corrected_wide,
        on=[sample_col, grouping_col],
        how="left",
    )

    # Replace original TPM/count columns with corrected/scaled values
    corrected_long_df = replace_original_with_corrected_tpm_and_scale_counts(
        corrected_long_df,
        cyclo_tpm_col=cyclo_tpm_col,
        noncyclo_tpm_col=noncyclo_tpm_col,
        cyclo_tpm_corrected_col=cyclo_tpm_col + output_suffix,
        noncyclo_tpm_corrected_col=noncyclo_tpm_col + output_suffix,
    )

    # =========================
    # Save outputs
    # =========================

    if output_dir is not None:
        corrected_long_df.to_csv(
            os.path.join(
                output_dir,
                f"long_format_for_isoranker_N50_corrected_{grouping_col}.tsv.gz",
            ),
            sep="\t",
            index=False,
            compression="gzip",
        )

        metadata.to_csv(
            os.path.join(
                output_dir,
                f"sample_metadata_with_N50_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

        pca_before.to_csv(
            os.path.join(
                output_dir,
                f"pca_before_N50_correction_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

        pca_after.to_csv(
            os.path.join(
                output_dir,
                f"pca_after_N50_correction_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

    results = {
        "corrected_long_df": corrected_long_df,
        "metadata": metadata,
        "pca_before": pca_before,
        "pca_after": pca_after,
        "pca_model_before": pca_model_before,
        "pca_model_after": pca_model_after,
        "explained_before": explained_before,
        "explained_after": explained_after,
        "features_used_for_pca": list(keep_features),
    }

    if return_matrices:
        results["matrices"] = {
            "raw_tpm_matrix": raw_matrix,
            "log2_tpm_matrix": log_matrix,
            "corrected_log2_tpm_matrix": corrected_log_matrix,
            "corrected_tpm_like_matrix": corrected_tpm_like_matrix,
            "log_matrix_for_pca": log_matrix_for_pca,
            "corrected_log_matrix_for_pca": corrected_log_matrix_for_pca,
        }

    return results