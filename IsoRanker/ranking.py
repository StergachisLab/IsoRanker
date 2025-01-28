import numpy as np

percentiles = [99.9, 99.5, 99, 98, 95]

def calculate_ranks_for_sample(df, group_col):
    """
    Calculate ranks for each sample at specified percentiles of z-scores.

    Parameters:
    - df: DataFrame containing columns 'Sample', group_col, 'z_score_of_test_stat', and 'test_statistic'.
    - group_col: The column name to group by (e.g., 'Isoform_PBid').
    - percentiles: List of percentiles to calculate ranks for (default: [99.5, 99.9]).

    Returns:
    - DataFrame with new rank columns added for each percentile.
    """
    # Filter out rows with NaN z-scores
    filtered_df = df[df['z_score_of_test_stat'].notna()].copy()

    # Initialize all rank columns globally
    for percentile in percentiles:
        rank_column_name = f'rank_top_{str(percentile).replace(".", "_")}_percentile'
        filtered_df[rank_column_name] = np.nan

    # Process each sample independently
    for sample in filtered_df['Sample'].unique():
        sample_df = filtered_df[filtered_df['Sample'] == sample]

        for percentile in percentiles:
            # Calculate percentile value
            percentile_value = np.percentile(sample_df['z_score_of_test_stat'], percentile)

            # Define rank column name
            rank_column_name = f'rank_top_{str(percentile).replace(".", "_")}_percentile'

            # Get top percentile subset and rank
            top_percentile_df = sample_df[sample_df['z_score_of_test_stat'] > percentile_value]
            sorted_df = top_percentile_df.sort_values(by='test_statistic', ascending=False).reset_index(drop=True)
            rank_mapping = {key: rank + 1 for rank, key in enumerate(sorted_df[group_col])}

            # Update the main DataFrame with ranks
            filtered_df.loc[
                (filtered_df['Sample'] == sample) & 
                (filtered_df[group_col].isin(rank_mapping.keys())),
                rank_column_name
            ] = filtered_df[group_col].map(rank_mapping)

    return filtered_df