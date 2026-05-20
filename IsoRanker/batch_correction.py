import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ============================================================
# Parse metadata with automatic covariate detection
# ============================================================

def clean_covariate_name(name):
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_]", "", name)
    return name


def parse_covariate_metadata(
    metadata_df,
    sample_col="Individual",
    condition_col="Sample type",
):
    """
    Detects columns named like:

        CONTINUOUS: Read length N50
        CATEGORICAL: CRDR

    and returns:
        metadata_df
        continuous_correction_cols
        categorical_correction_cols
    """

    metadata_df = metadata_df.copy()
    metadata_df.columns = metadata_df.columns.str.strip()

    continuous_cols = []
    categorical_cols = {}
    rename_dict = {}

    for col in metadata_df.columns:

        if col.startswith("CONTINUOUS:"):
            clean_name = clean_covariate_name(
                col.replace("CONTINUOUS:", "").strip()
            )
            rename_dict[col] = clean_name
            continuous_cols.append(clean_name)

        elif col.startswith("CATEGORICAL:"):
            clean_name = clean_covariate_name(
                col.replace("CATEGORICAL:", "").strip()
            )
            rename_dict[col] = clean_name
            categorical_cols[col] = clean_name

    metadata_df = metadata_df.rename(columns=rename_dict)

    metadata_df = metadata_df.rename(columns={
        sample_col: "Base_Sample",
        condition_col: "Condition",
    })

    metadata_df["Condition"] = (
        metadata_df["Condition"]
        .astype(str)
        .str.strip()
        .str.lower()
        .map({"cyclo": "Cyclo", "noncyclo": "Noncyclo"})
    )

    if metadata_df["Condition"].isna().any():
        bad = metadata_df[metadata_df["Condition"].isna()]
        raise ValueError(
            "Some condition values could not be mapped to Cyclo/Noncyclo:\n"
            f"{bad}"
        )

    for col in continuous_cols:
        metadata_df[col] = pd.to_numeric(metadata_df[col], errors="coerce")

    categorical_cols = list(categorical_cols.values())

    return metadata_df, continuous_cols, categorical_cols


# ============================================================
# General correction helper
# ============================================================

def remove_covariate_effects_preserve_variables(
    matrix,
    metadata,
    sample_key_col="Sample_Condition",
    preserve_cols=("Condition",),
    continuous_correction_cols=None,
    categorical_correction_cols=None,
):
    continuous_correction_cols = list(continuous_correction_cols or [])
    categorical_correction_cols = list(categorical_correction_cols or [])
    preserve_cols = list(preserve_cols or [])

    meta = metadata.set_index(sample_key_col).loc[matrix.index].copy()

    intercept = pd.DataFrame(
        {"Intercept": np.ones(matrix.shape[0])},
        index=matrix.index,
    )

    keep_parts = [intercept]
    remove_parts = []

    # Variables to preserve
    for col in preserve_cols:
        if col not in meta.columns:
            raise ValueError(f"Preserve column not found in metadata: {col}")

        if pd.api.types.is_numeric_dtype(meta[col]):
            x = pd.to_numeric(meta[col], errors="coerce")

            if x.isna().any():
                raise ValueError(f"Preserve column has missing values: {col}")

            sd = x.std()

            if sd == 0 or pd.isna(sd):
                raise ValueError(f"Preserve column has zero variance: {col}")

            keep_parts.append(
                pd.DataFrame(
                    {f"{col}_z": (x - x.mean()) / sd},
                    index=matrix.index,
                )
            )

        else:
            keep_parts.append(
                pd.get_dummies(
                    meta[col].astype(str),
                    prefix=col,
                    drop_first=True,
                    dtype=float,
                )
            )

    # Continuous covariates to remove
    for col in continuous_correction_cols:
        if col not in meta.columns:
            raise ValueError(f"Continuous correction column not found: {col}")

        x = pd.to_numeric(meta[col], errors="coerce")

        if x.isna().any():
            raise ValueError(f"Continuous correction column has missing values: {col}")

        sd = x.std()

        if sd == 0 or pd.isna(sd):
            raise ValueError(f"Continuous correction column has zero variance: {col}")

        remove_parts.append(
            pd.DataFrame(
                {f"{col}_z": (x - x.mean()) / sd},
                index=matrix.index,
            )
        )

    # Categorical covariates to remove
    for col in categorical_correction_cols:
        if col not in meta.columns:
            raise ValueError(f"Categorical correction column not found: {col}")

        if meta[col].isna().any():
            raise ValueError(f"Categorical correction column has missing values: {col}")

        remove_parts.append(
            pd.get_dummies(
                meta[col].astype(str),
                prefix=col,
                drop_first=True,
                dtype=float,
            )
        )

    X_keep = pd.concat(keep_parts, axis=1)

    if remove_parts:
        X_remove = pd.concat(remove_parts, axis=1)
        X_full = pd.concat([X_keep, X_remove], axis=1)
    else:
        X_full = X_keep.copy()

    Y = matrix.values

    beta_full = np.linalg.pinv(X_full.values) @ Y
    beta_keep = np.linalg.pinv(X_keep.values) @ Y

    fitted_full = X_full.values @ beta_full
    fitted_keep = X_keep.values @ beta_keep

    corrected = fitted_keep + (Y - fitted_full)

    return pd.DataFrame(
        corrected,
        index=matrix.index,
        columns=matrix.columns,
    )


# ============================================================
# Replace TPM and count columns
# ============================================================

def replace_original_with_corrected_tpm_and_scale_counts(
    long_df,
    cyclo_tpm_col="Cyclo_TPM",
    noncyclo_tpm_col="Noncyclo_TPM",
    cyclo_tpm_corrected_col="Cyclo_TPM_corrected",
    noncyclo_tpm_corrected_col="Noncyclo_TPM_corrected",
):
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

    long_df[cyclo_tpm_col + "_original"] = long_df[cyclo_tpm_col]
    long_df[noncyclo_tpm_col + "_original"] = long_df[noncyclo_tpm_col]

    long_df[cyclo_tpm_col + "_original"] = pd.to_numeric(
        long_df[cyclo_tpm_col + "_original"],
        errors="coerce",
    ).fillna(0)

    long_df[noncyclo_tpm_col + "_original"] = pd.to_numeric(
        long_df[noncyclo_tpm_col + "_original"],
        errors="coerce",
    ).fillna(0)

    long_df[cyclo_tpm_corrected_col] = pd.to_numeric(
        long_df[cyclo_tpm_corrected_col],
        errors="coerce",
    )

    long_df[noncyclo_tpm_corrected_col] = pd.to_numeric(
        long_df[noncyclo_tpm_corrected_col],
        errors="coerce",
    )

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

    long_df[cyclo_tpm_col] = long_df[cyclo_tpm_corrected_col]
    long_df[noncyclo_tpm_col] = long_df[noncyclo_tpm_corrected_col]

    for col in cyclo_count_cols:
        if col in long_df.columns:
            long_df[col + "_original"] = long_df[col]
            long_df[col] = (
                pd.to_numeric(long_df[col], errors="coerce").fillna(0)
                * long_df["Cyclo_scaling_factor"]
            )

    for col in noncyclo_count_cols:
        if col in long_df.columns:
            long_df[col + "_original"] = long_df[col]
            long_df[col] = (
                pd.to_numeric(long_df[col], errors="coerce").fillna(0)
                * long_df["Noncyclo_scaling_factor"]
            )

    return long_df


# ============================================================
# Main function
# ============================================================

def correct_tpm_by_metadata_covariates_with_pca(
    long_df,
    metadata_df,
    grouping_col="Isoform",
    sample_col="Sample",
    cyclo_tpm_col="Cyclo_TPM",
    noncyclo_tpm_col="Noncyclo_TPM",
    metadata_sample_col="Individual",
    metadata_condition_col="Sample type",
    preserve_cols=("Condition",),
    tpm_filter=10,
    output_dir=None,
    output_suffix="_corrected",
    cmap="viridis",
    make_plots=True,
    return_matrices=False,
):
    """
    Correct TPM values using automatically detected metadata covariates and run PCA
    before and after correction.

    This function is designed for IsoRanker-style long-format expression data.
    It builds a sample-by-feature TPM matrix from Cyclo_TPM and Noncyclo_TPM,
    corrects log2(TPM + 1) values for technical covariates, converts corrected
    values back to TPM-like values, and returns an IsoRanker-ready dataframe.

    Metadata covariates are detected automatically using column prefixes:

        CONTINUOUS: Read length N50
        CATEGORICAL: CRDR

    The correction model is:

        log2(TPM + 1) ~ preserved variables + correction covariates

    The returned dataframe:
        - replaces Cyclo_TPM and Noncyclo_TPM with corrected TPM-like values
        - scales count columns using the same TPM correction factors
        - preserves original TPM/count columns with "_original" suffix

    Parameters
    ----------
    long_df : pandas.DataFrame
        Long-format expression dataframe.

        Required columns:
            - sample_col
            - grouping_col
            - cyclo_tpm_col
            - noncyclo_tpm_col

        Typical IsoRanker input columns:
            - "Sample"
            - "Isoform"
            - "Cyclo_TPM"
            - "Noncyclo_TPM"

        If present, these count columns are scaled:
            - H0_cyclo_count
            - H1_cyclo_count
            - H2_cyclo_count
            - cyclo_count
            - H0_noncyclo_count
            - H1_noncyclo_count
            - H2_noncyclo_count
            - noncyclo_count

    metadata_df : pandas.DataFrame
        Sample-condition metadata dataframe.

        Required columns by default:
            - "Individual"
            - "Sample type"

        Correction covariates should be encoded in the column names:

            CONTINUOUS: Read length N50
            CATEGORICAL: CRDR

        Continuous covariates are z-scored before regression.
        Categorical covariates are one-hot encoded with drop_first=True.

    grouping_col : str, default "Isoform"
        Feature column used to build the sample-by-feature matrix.

        Use:
            - "Isoform" for isoform-level correction
            - "associated_gene" for gene-level correction

    sample_col : str, default "Sample"
        Sample identifier column in `long_df`.

        Must match `metadata_df[metadata_sample_col]`.

    cyclo_tpm_col : str, default "Cyclo_TPM"
        TPM column for cycloheximide-treated samples.

        In the returned dataframe:
            - this column is replaced with corrected values
            - original values are saved as "Cyclo_TPM_original"

    noncyclo_tpm_col : str, default "Noncyclo_TPM"
        TPM column for untreated samples.

        In the returned dataframe:
            - this column is replaced with corrected values
            - original values are saved as "Noncyclo_TPM_original"

    metadata_sample_col : str, default "Individual"
        Sample identifier column in `metadata_df`.

    metadata_condition_col : str, default "Sample type"
        Condition column in `metadata_df`.

        Accepted values are case-insensitive:
            - cyclo
            - noncyclo

    preserve_cols : tuple or list of str, default ("Condition",)
        Variables to preserve during correction.

        These variables are included in both the full model and the keep model.
        Usually this should include "Condition" so Cyclo vs Noncyclo differences
        are not regressed out.

    tpm_filter : int or float, default 10
        Feature filter used only for PCA.

        A feature is included in PCA if at least one sample-condition has:

            TPM > tpm_filter

        This does not affect which features are corrected.

    output_dir : str or None, default None
        Directory where outputs are saved.

        If provided, saves:
            - IsoRanker-ready corrected dataframe
            - metadata with detected covariates
            - PCA table before correction
            - PCA table after correction
            - PCA PDF before correction
            - PCA PDF after correction

        If None, no files are written.

    output_suffix : str, default "_corrected"
        Suffix for intermediate corrected TPM columns.

        With default TPM names:
            - Cyclo_TPM_corrected
            - Noncyclo_TPM_corrected

    cmap : str, default "viridis"
        Matplotlib colormap for PCA coloring.

        The first detected continuous covariate is used for PCA coloring.

    make_plots : bool, default True
        Whether to display PCA plots with plt.show().

        Set to False for cluster or batch jobs.

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
        Dictionary containing corrected dataframe, PCA outputs, detected covariates,
        and optional matrices.

        Main keys:
            corrected_long_df
            metadata
            pca_before
            pca_after
            pca_model_before
            pca_model_after
            explained_before
            explained_after
            features_used_for_pca
            continuous_correction_cols
            categorical_correction_cols

    Notes
    -----
    The returned `corrected_long_df` can be used directly as IsoRanker input.
    """

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    required_long_cols = {
        sample_col,
        grouping_col,
        cyclo_tpm_col,
        noncyclo_tpm_col,
    }

    missing_long = required_long_cols - set(long_df.columns)

    if missing_long:
        raise ValueError(f"long_df is missing required columns: {missing_long}")

    covariates, continuous_correction_cols, categorical_correction_cols = (
        parse_covariate_metadata(
            metadata_df=metadata_df,
            sample_col=metadata_sample_col,
            condition_col=metadata_condition_col,
        )
    )

    if not continuous_correction_cols and not categorical_correction_cols:
        raise ValueError(
            "No correction covariates detected. "
            "Use columns like 'CONTINUOUS: Read length N50' or 'CATEGORICAL: CRDR'."
        )

    print("Detected continuous correction covariates:")
    print(continuous_correction_cols)

    print("Detected categorical correction covariates:")
    print(categorical_correction_cols)

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

    analysis_metadata = (
        long_expr[["Sample_Condition", "Base_Sample", "Condition"]]
        .drop_duplicates()
        .set_index("Sample_Condition")
        .loc[raw_matrix.index]
        .reset_index()
    )

    analysis_metadata = analysis_metadata.merge(
        covariates,
        on=["Base_Sample", "Condition"],
        how="left",
    )

    needed_cols = (
        list(preserve_cols)
        + list(continuous_correction_cols)
        + list(categorical_correction_cols)
    )

    for col in needed_cols:
        if col not in analysis_metadata.columns:
            raise ValueError(f"Required metadata column not found: {col}")

        if analysis_metadata[col].isna().any():
            missing_rows = analysis_metadata.loc[
                analysis_metadata[col].isna(),
                ["Sample_Condition", "Base_Sample", "Condition"],
            ]

            raise ValueError(
                f"Metadata column has missing values: {col}\n{missing_rows}"
            )

    log_matrix = np.log2(raw_matrix + 1)

    corrected_log_matrix = remove_covariate_effects_preserve_variables(
        matrix=log_matrix,
        metadata=analysis_metadata,
        sample_key_col="Sample_Condition",
        preserve_cols=preserve_cols,
        continuous_correction_cols=continuous_correction_cols,
        categorical_correction_cols=categorical_correction_cols,
    )

    corrected_tpm_like_matrix = 2 ** corrected_log_matrix

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
            analysis_metadata,
            on="Sample_Condition",
            how="left",
        )

        return pca_df, pca, explained

    pca_before, pca_model_before, explained_before = run_pca(log_matrix_for_pca)
    pca_after, pca_model_after, explained_after = run_pca(corrected_log_matrix_for_pca)

    color_col = continuous_correction_cols[0] if continuous_correction_cols else None

    def plot_pca(pca_df, explained, title, output_pdf=None):
        fig, ax = plt.subplots(figsize=(10, 8))

        markers = {
            "Cyclo": "o",
            "Noncyclo": "X",
        }

        scatter_for_colorbar = None

        for condition, marker in markers.items():
            subset = pca_df[pca_df["Condition"] == condition].copy()

            if color_col is not None:
                sc = ax.scatter(
                    subset["PC1"],
                    subset["PC2"],
                    c=subset[color_col],
                    cmap=cmap,
                    vmin=pca_df[color_col].min(),
                    vmax=pca_df[color_col].max(),
                    marker=marker,
                    s=120,
                    edgecolor="black",
                    linewidth=0.5,
                    label=condition,
                )
                scatter_for_colorbar = sc
            else:
                ax.scatter(
                    subset["PC1"],
                    subset["PC2"],
                    marker=marker,
                    s=120,
                    edgecolor="black",
                    linewidth=0.5,
                    label=condition,
                )

        for _, row in pca_df.iterrows():
            ax.text(
                row["PC1"],
                row["PC2"],
                row["Base_Sample"],
                fontsize=6,
                ha="right",
                va="bottom",
            )

        if scatter_for_colorbar is not None:
            cbar = fig.colorbar(scatter_for_colorbar, ax=ax)
            cbar.set_label(color_col)

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
            f"pca_before_metadata_covariate_correction_{grouping_col}_any_TPM_gt_{tpm_filter}.pdf",
        )

        after_pdf = os.path.join(
            output_dir,
            f"pca_after_metadata_covariate_correction_{grouping_col}_any_TPM_gt_{tpm_filter}.pdf",
        )
    else:
        before_pdf = None
        after_pdf = None

    plot_pca(
        pca_before,
        explained_before,
        title=f"Before metadata covariate correction: {grouping_col}",
        output_pdf=before_pdf,
    )

    plot_pca(
        pca_after,
        explained_after,
        title=f"After metadata covariate correction: {grouping_col}",
        output_pdf=after_pdf,
    )

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
        analysis_metadata[["Sample_Condition", "Base_Sample", "Condition"]],
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

    corrected_wide = corrected_wide.rename(
        columns={
            "Base_Sample": sample_col,
            "Cyclo": cyclo_tpm_col + output_suffix,
            "Noncyclo": noncyclo_tpm_col + output_suffix,
        }
    )

    corrected_long_df = long_df.merge(
        corrected_wide,
        on=[sample_col, grouping_col],
        how="left",
    )

    corrected_long_df = replace_original_with_corrected_tpm_and_scale_counts(
        corrected_long_df,
        cyclo_tpm_col=cyclo_tpm_col,
        noncyclo_tpm_col=noncyclo_tpm_col,
        cyclo_tpm_corrected_col=cyclo_tpm_col + output_suffix,
        noncyclo_tpm_corrected_col=noncyclo_tpm_col + output_suffix,
    )

    if output_dir is not None:
        corrected_long_df.to_csv(
            os.path.join(
                output_dir,
                f"long_format_for_isoranker_metadata_covariate_corrected_{grouping_col}.tsv.gz",
            ),
            sep="\t",
            index=False,
            compression="gzip",
        )

        analysis_metadata.to_csv(
            os.path.join(
                output_dir,
                f"sample_metadata_with_detected_covariates_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

        pca_before.to_csv(
            os.path.join(
                output_dir,
                f"pca_before_metadata_covariate_correction_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

        pca_after.to_csv(
            os.path.join(
                output_dir,
                f"pca_after_metadata_covariate_correction_{grouping_col}.tsv",
            ),
            sep="\t",
            index=False,
        )

    results = {
        "corrected_long_df": corrected_long_df,
        "metadata": analysis_metadata,
        "pca_before": pca_before,
        "pca_after": pca_after,
        "pca_model_before": pca_model_before,
        "pca_model_after": pca_model_after,
        "explained_before": explained_before,
        "explained_after": explained_after,
        "features_used_for_pca": list(keep_features),
        "continuous_correction_cols": continuous_correction_cols,
        "categorical_correction_cols": categorical_correction_cols,
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

