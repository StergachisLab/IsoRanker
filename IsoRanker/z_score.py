def calculate_z_score(df, group_col, stat_col):
    """
    Calculate z-scores within each group.
    
    Parameters:
    - df: Input DataFrame.
    - group_col: Column to group by (e.g., Isoform_PBid).
    - stat_col: Column containing the test statistic.
    """
    def z_score_func(group):
        group = group.copy()
        for i, row in group.iterrows():
            others = group.loc[group.index != i, stat_col]
            median_others = others.median()
            sd_all = group[stat_col].std()
            group.at[i, 'z_score_of_test_stat'] = 0 if sd_all == 0 else (row[stat_col] - median_others) / sd_all
        return group

    return df.groupby(group_col).apply(z_score_func).reset_index(drop=True)
