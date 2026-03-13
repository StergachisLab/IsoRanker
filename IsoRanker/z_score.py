def calculate_z_score(df, group_col, stat_col):
    """
    Calculate z-scores within each group.
    
    Parameters:
    - df: Input DataFrame.
    - group_col: Column to group by (e.g., Isoform_PBid).
    - stat_col: Column containing the test statistic.
    """
    def z_score_func(group):
        # Work on a copy to avoid mutating the original group
        group = group.copy()
        for i, row in group.iterrows():
            others = group.loc[group.index != i, stat_col]
            median_others = others.median()
            sd_all = group[stat_col].std()
            group.at[i, 'z_score_of_test_stat'] = (
                0 if sd_all == 0 else (row[stat_col] - median_others) / sd_all
            )
        return group

    # Need to do this wrapper because using apply() on a group no longer enables the group column to be included in the output
    def wrapper(group):
        # Apply the inner z-score function
        out = z_score_func(group)

        # Ensure we don't accidentally mutate the original
        out = out.copy()

        # If the grouping column is missing (when include_groups=False),
        # add it back using the group key (group.name)
        if group_col not in out.columns:
            out[group_col] = group.name

        
        # Move grouping column to the front
        cols = [group_col] + [c for c in out.columns if c != group_col]
        out = out[cols]

        return out

    result = (
        df.groupby(group_col)
          .apply(wrapper, include_groups=False)
          .reset_index(drop=True)
    )

    return result


import numpy as np

# Using MAD instead of SD
def calculate_z_score_MAD(df, group_col, stat_col, variance_floor=None):
    """
    Calculate robust leave-one-out standardized scores within each group
    using the median absolute deviation (MAD).

    Parameters:
    - df: Input DataFrame.
    - group_col: Column to group by.
    - stat_col: Column containing the test statistic.
    - variance_floor: Minimum MAD value to use in denominator.
                      If None, defaults to 0.1.
    """
    def z_score_func(group):
        # Work on a copy to avoid mutating the original group
        group = group.copy()

        floor = 0.1 if variance_floor is None else variance_floor

        for i, row in group.iterrows():
            others = group.loc[group.index != i, stat_col]

            if len(others) == 0:
                group.at[i, 'z_score_of_test_stat'] = 0
                continue

            median_others = others.median()
            mad_others = np.median(np.abs(others - median_others))
            mad_others = max(mad_others, floor)

            group.at[i, 'z_score_of_test_stat'] = (
                0 if mad_others == 0 else (row[stat_col] - median_others) / (1.4826 * mad_others)
            )

        return group

    # Wrapper to restore grouping column, matching calculate_z_score()
    def wrapper(group):
        out = z_score_func(group)
        out = out.copy()

        if group_col not in out.columns:
            out[group_col] = group.name

        cols = [group_col] + [c for c in out.columns if c != group_col]
        out = out[cols]

        return out

    result = (
        df.groupby(group_col)
          .apply(wrapper, include_groups=False)
          .reset_index(drop=True)
    )

    return result