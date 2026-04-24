import numpy as np


def _nanmedian_empty_safe(x):
    x = x[~np.isnan(x)]
    if x.size == 0:
        return np.nan
    return np.median(x)


def _leave_one_out_median(values):
    values = values.astype(float, copy=False)
    n = values.size
    out = np.empty(n, dtype=float)

    if n <= 1:
        out[:] = np.nan
        return out

    for i in range(n):
        out[i] = _nanmedian_empty_safe(np.delete(values, i))

    return out


def _leave_one_out_mad(values, floor):
    values = values.astype(float, copy=False)
    n = values.size

    medians = np.empty(n, dtype=float)
    mads = np.empty(n, dtype=float)

    if n <= 1:
        medians[:] = np.nan
        mads[:] = np.nan
        return medians, mads

    for i in range(n):
        others = np.delete(values, i)

        median_others = _nanmedian_empty_safe(others)

        if np.isnan(median_others):
            mad_others = np.nan
        else:
            mad_others = _nanmedian_empty_safe(np.abs(others - median_others))

        if np.isnan(mad_others):
            mad_others = floor
        else:
            mad_others = max(mad_others, floor)

        medians[i] = median_others
        mads[i] = mad_others

    return medians, mads


def calculate_z_score(df, group_col, stat_col):
    """
    Faster version of the original leave-one-out median z-score.

    Original formula:
        z = (current_value - median(other samples)) / std(all samples)
    """
    out = df.copy()
    values = out[stat_col].to_numpy(dtype=float)
    z = np.empty(len(out), dtype=float)
    z[:] = np.nan

    groups = out.groupby(group_col, sort=False).indices

    for _, idx in groups.items():
        idx = np.asarray(idx)

        group_values = values[idx]

        medians = _leave_one_out_median(group_values)

        # pandas std() uses ddof=1 by default
        sd_all = np.nanstd(group_values, ddof=1)

        if sd_all == 0:
            z[idx] = 0
        else:
            z[idx] = (group_values - medians) / sd_all

    out["z_score_of_test_stat"] = z
    return out


def calculate_z_score_MAD(df, group_col, stat_col, variance_floor=None):
    """
    Faster version of the original robust leave-one-out MAD z-score.

    Original formula:
        z = (current_value - median(other samples)) /
            (1.4826 * max(MAD(other samples), floor))
    """
    out = df.copy()
    values = out[stat_col].to_numpy(dtype=float)
    z = np.empty(len(out), dtype=float)
    z[:] = np.nan

    floor = 0.1 if variance_floor is None else variance_floor

    groups = out.groupby(group_col, sort=False).indices

    for _, idx in groups.items():
        idx = np.asarray(idx)

        group_values = values[idx]

        if group_values.size <= 1:
            z[idx] = 0
            continue

        medians, mads = _leave_one_out_mad(group_values, floor)

        denom = 1.4826 * mads

        z[idx] = np.where(
            denom == 0,
            0,
            (group_values - medians) / denom
        )

    out["z_score_of_test_stat"] = z
    return out