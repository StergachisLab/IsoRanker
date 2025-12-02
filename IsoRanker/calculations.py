def apply_hypothesis_test(df, group_col, test_statistic_func):
    """
    Apply a hypothesis test statistic function to grouped data.

    Parameters:
    - df: Input DataFrame.
    - group_col: Column to group by (e.g., Isoform_PBid).
    - test_statistic_func: Function to calculate the test statistic for each group.
    """

    # Need to do this wrapper because using apply() on a group no longer enables the group column to be included in the output
    def wrapper(group):
        # Apply the user-provided function
        out = test_statistic_func(group)

        # Ensure we don't accidentally mutate the original
        out = out.copy()

        # If the grouping column is missing (which it will be when include_groups=False),
        # add it back using the group key (group.name)
        if group_col not in out.columns:
            out[group_col] = group.name

        
        # Move grouping column to the front
        cols = [group_col] + [c for c in out.columns if c != group_col]
        out = out[cols]

        return out

    # Use include_groups=False to be future-proof:
    # - groups passed into wrapper DO NOT have group_col as a column
    # - group.name holds the group key (e.g., the Isoform_PBid)
    result = (
        df.groupby(group_col)
          .apply(wrapper, include_groups=False)
          .reset_index(drop=True)
    )

    return result
