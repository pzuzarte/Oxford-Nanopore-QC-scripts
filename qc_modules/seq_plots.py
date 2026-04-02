"""
Core sequencing QC plot functions.

Each function accepts a DataFrame (and optional kwargs), produces a plot,
saves it to disk, and returns the output file path.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Apply the seaborn darkgrid theme globally for all plots in this module
sns.set_theme(style="darkgrid")

_PAL       = sns.color_palette()   # default "deep" palette
PASS_COLOR = _PAL[0]               # blue
FAIL_COLOR = _PAL[3]               # red
N50_COLOR  = _PAL[3]               # red (same as fail for N50 marker)

_FIGSIZE = (10, 8)


def _save(fig, path):
    fig.savefig(path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Yield over time
# ---------------------------------------------------------------------------

def plot_yield_over_time(df, outpath, max_points=50_000):
    """Cumulative pass-filter yield (Gb) vs run time (hr)."""
    pf = df[df['passes_filtering'] == True].sort_values('start_time')
    pf = pf.assign(cumulative_yield=pf['sequence_length_template'].cumsum())

    sample = pf if len(pf) <= max_points else pf.sample(max_points, random_state=42)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.scatter(
        sample['start_time'] / 3600,
        sample['cumulative_yield'] / 1e9,
        s=4, alpha=0.05, color=PASS_COLOR,
    )
    ax.set_title('Pass-Filter Cumulative Yield vs. Time', fontsize=16)
    ax.set_xlabel('Time (hr)', fontsize=13)
    ax.set_ylabel('Yield (Gb)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Reads per hour
# ---------------------------------------------------------------------------

def plot_reads_per_hour(df, outpath):
    """Histogram of read start times binned by hour."""
    run_time_hr = df['start_time'].max() / 3600
    n_bins = max(1, int(round(run_time_hr)))

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    sns.histplot(
        df['start_time'] / 3600,
        bins=n_bins,
        ax=ax,
    )
    ax.set_title('Reads per Hour', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Read Count', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length distribution
# ---------------------------------------------------------------------------

def plot_read_length_dist(df, n50, outpath, max_length=None, min_length=None):
    """Histogram of read lengths with N50 marked."""
    pf = df[df['passes_filtering'] == True]
    mean_len = pf['sequence_length_template'].mean()
    sd_len = pf['sequence_length_template'].std()
    if max_length is None:
        max_length = int(mean_len + 3 * sd_len)
    if min_length is None:
        min_length = 0

    plot_df = df[
        (df['sequence_length_template'] <= max_length) &
        (df['sequence_length_template'] >= min_length)
    ]

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    sns.histplot(
        plot_df['sequence_length_template'],
        bins=100,
        ax=ax,
    )
    ax.axvline(n50, color=N50_COLOR, linestyle='--', alpha=0.8, linewidth=1.5)
    ax.text(n50, ax.get_ylim()[1] * 0.95, f' N50 = {n50:,} bp',
            ha='left', va='top', color=N50_COLOR, fontsize=11)
    ax.set_xlim(min_length, max_length)
    ax.set_title('Read Length Distribution', fontsize=16)
    ax.set_xlabel('Read Length (bp)', fontsize=13)
    ax.set_ylabel('Read Count', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Proportion of bases by read length (N50 plot)
# ---------------------------------------------------------------------------

def plot_length_proportions(df, n50, outpath, max_points=100_000, min_length=None, max_length=None):
    """Proportion of PF bases sequenced by read length."""
    pf = df[df['passes_filtering'] == True].sort_values(
        'sequence_length_template', ascending=False
    ).copy()
    if min_length is not None:
        pf = pf[pf['sequence_length_template'] >= min_length]
    if max_length is not None:
        pf = pf[pf['sequence_length_template'] <= max_length]
    total_yield = pf['sequence_length_template'].sum()
    pf['cumulative_yield'] = pf['sequence_length_template'].cumsum()
    pf['proportion'] = pf['cumulative_yield'] / total_yield

    sample = pf if len(pf) <= max_points else pf.sample(max_points, random_state=42)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.scatter(
        sample['sequence_length_template'],
        sample['proportion'],
        s=4, alpha=0.4, color=PASS_COLOR,
    )
    ax.axvline(n50, color=N50_COLOR, linestyle='--', alpha=0.8, linewidth=1.5)
    ax.axhline(0.5, color=N50_COLOR, linestyle='--', alpha=0.8, linewidth=1.5)
    ax.text(n50, 1.0, f' N50 = {n50:,} bp', ha='left', va='top',
            color=N50_COLOR, fontsize=11)
    ax.set_title('PF Bases Sequenced by Read Length', fontsize=16)
    ax.set_xlabel('Read Length (bp)', fontsize=13)
    ax.set_ylabel('Proportion of Bases Sequenced', fontsize=13)
    ax.set_ylim(0, 1.05)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Q-score distribution
# ---------------------------------------------------------------------------

def plot_qscore_dist(df, outpath):
    """KDE of Q-scores, split by pass/fail."""
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    sns.kdeplot(
        df.loc[df['passes_filtering'] == True, 'mean_qscore_template'],
        color=PASS_COLOR, fill=True, alpha=0.4, label='Pass', ax=ax,
    )
    sns.kdeplot(
        df.loc[df['passes_filtering'] == False, 'mean_qscore_template'],
        color=FAIL_COLOR, fill=True, alpha=0.4, label='Fail', ax=ax,
    )
    ax.set_title('Q-Score Distribution', fontsize=16)
    ax.set_xlabel('Mean Q-Score', fontsize=13)
    ax.set_ylabel('Density', fontsize=13)
    ax.legend(fontsize=12)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length vs Q-score joint plot
# ---------------------------------------------------------------------------

def plot_length_vs_qscore(df, outpath, max_points=20_000, min_length=None, max_length=None):
    """Hex joint plot: read length vs Q-score (pass-filter reads)."""
    pf = df[df['passes_filtering'] == True]
    mean_len = pf['sequence_length_template'].mean()
    sd_len = pf['sequence_length_template'].std()
    if max_length is None:
        max_length = int(mean_len + 3 * sd_len)
    if min_length is None:
        min_length = 0

    sample = pf if len(pf) <= max_points else pf.sample(max_points, random_state=42)
    sample = sample[sample['sequence_length_template'] >= min_length]

    grid = sns.jointplot(
        x='sequence_length_template',
        y='mean_qscore_template',
        data=sample,
        kind='hex',
        height=8,
        xlim=(min_length, max_length),
        ylim=(0, 40),
    )
    grid.set_axis_labels('Read Length (bp)', 'Mean Q-Score', fontsize=13)
    grid.figure.suptitle('Read Length vs. Q-Score (PF reads)', y=1.01, fontsize=16)
    grid.figure.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(grid.figure)
    return outpath


# ---------------------------------------------------------------------------
# Sequencing speed over time
# ---------------------------------------------------------------------------

def plot_speed_over_time(df, outpath, max_points=50_000):
    """Scatter: sequencing speed (bp/s) vs run time, coloured by pass/fail."""
    sample = df if len(df) <= max_points else df.sample(max_points, random_state=42)
    sample = sample.dropna(subset=['sequencing_speed'])

    fig, ax = plt.subplots(figsize=(12, 8))
    palette = {True: PASS_COLOR, False: FAIL_COLOR}
    for passed, grp in sample.groupby('passes_filtering'):
        ax.scatter(
            grp['start_time'] / 3600,
            grp['sequencing_speed'],
            s=4, alpha=0.3,
            color=palette[passed],
            label='Pass' if passed else 'Fail',
        )
    ax.set_title('Sequencing Speed over Time', fontsize=16)
    ax.set_xlabel('Start Time (hr)', fontsize=13)
    ax.set_ylabel('Sequencing Speed (bp/s)', fontsize=13)
    ax.legend(fontsize=12)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Median current over time
# ---------------------------------------------------------------------------

def plot_median_current(df, outpath, max_points=50_000):
    """Scatter: median template current vs run time."""
    if 'median_template' not in df.columns:
        return None

    sample = df if len(df) <= max_points else df.sample(max_points, random_state=42)
    palette = {True: PASS_COLOR, False: FAIL_COLOR}

    fig, ax = plt.subplots(figsize=(12, 8))
    for passed, grp in sample.groupby('passes_filtering'):
        ax.scatter(
            grp['start_time'] / 3600,
            grp['median_template'],
            s=4, alpha=0.3,
            color=palette[passed],
            label='Pass' if passed else 'Fail',
        )
    ax.set_title('Median Current over Time', fontsize=16)
    ax.set_xlabel('Start Time (hr)', fontsize=13)
    ax.set_ylabel('Median Current (pA)', fontsize=13)
    ax.legend(fontsize=12)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Flowcell channel map  (3-panel: Q-score / read count / yield)
# ---------------------------------------------------------------------------

def _detect_flowcell_layout(max_channel):
    """
    Return (nrows, ncols) for the physical flowcell grid.
      MinION / GridION : ≤  512 channels  →  16 rows ×  32 cols  (sequential)
      PromethION       : >  512 channels  →  25 rows × 120 cols  (block-based)

    PromethION block layout (matches MinKNOW UI):
      12 blocks of 250 channels; each block occupies 10 cols × 25 rows.
        block  = (ch - 1) // 250      → which 10-column block (0-11)
        rem    = (ch - 1) % 250
        row    = rem // 10            → 0-24
        col    = rem % 10 + block*10  → 0-119
    Channels 1-10 fill row 0 of block 0 (cols 0-9),
    channels 11-20 fill row 1 of block 0, …,
    channels 251-260 fill row 0 of block 1 (cols 10-19), etc.
    """
    if max_channel <= 512:
        return 16, 32
    else:
        return 25, 120


def plot_channel_map(df, outpath):
    """
    3-panel heatmap of the flowcell showing per-channel mean Q-score,
    read count, and yield (Mb).  Grid layout is auto-detected:
      MinION/GridION → 32 cols × 16 rows
      PromethION     → 120 cols × 25 rows
    Channel 1 is always top-left; cells are rendered square.
    """
    pf = df[df['passes_filtering'] == True]

    grp = pf.groupby('channel').agg(
        mean_q=('mean_qscore_template', 'mean'),
        read_count=('read_id', 'count'),
        yield_mb=('sequence_length_template', lambda x: x.sum() / 1e6),
    ).reset_index()

    max_channel = int(df['channel'].max())
    nrows, ncols = _detect_flowcell_layout(max_channel)

    def make_grid(value_col):
        grid = np.full((nrows, ncols), np.nan)
        for _, row in grp.iterrows():
            ch = int(row['channel'])
            if max_channel <= 512:
                # MinION/GridION: simple sequential row-major layout
                r = (ch - 1) // ncols
                c = (ch - 1) % ncols
            else:
                # PromethION: 12 blocks of 250 channels, each block occupies
                # 10 columns × 25 rows — matches the MinKNOW channel map UI.
                block = (ch - 1) // 250
                rem   = (ch - 1) % 250
                r     = rem // 10
                c     = rem % 10 + block * 10
            if 0 <= r < nrows and 0 <= c < ncols:
                grid[r, c] = row[value_col]
        return grid

    q_grid     = make_grid('mean_q')
    count_grid = make_grid('read_count')
    yield_grid = make_grid('yield_mb')

    # Panels are stacked vertically so each spans the full report width.
    # Cell size is chosen per flowcell type so figures are appropriately scaled.
    #   PromethION (120 × 25): 0.18 in/cell → ~21.6 in wide data area
    #   MinION/GridION (32 × 16): 0.30 in/cell → ~9.6 in wide data area
    if max_channel <= 512:
        cell_inch = 0.30          # MinION / GridION
    else:
        cell_inch = 0.18          # PromethION
    panel_w = ncols * cell_inch
    panel_h = panel_w * (nrows / ncols)   # keeps cells square
    fig_w   = panel_w + 2.2               # + colorbar + margins
    fig_h   = 3 * panel_h + 3.0          # 3 stacked panels + spacing

    panels = [
        (q_grid,     'Mean Q-Score', 'Blues'),
        (count_grid, 'Read Count',   'Greens'),
        (yield_grid, 'Yield (Mb)',   'Oranges'),
    ]

    with sns.axes_style("white"):
        fig, axes = plt.subplots(3, 1, figsize=(fig_w, fig_h),
                                 constrained_layout=True)

    for ax, (data, title, cmap) in zip(axes, panels):
        masked = np.ma.masked_invalid(data)
        im = ax.imshow(masked, cmap=cmap,
                       aspect='equal',       # square cells
                       origin='upper',       # channel 1 at top-left
                       interpolation='nearest')
        plt.colorbar(im, ax=ax, shrink=0.8, pad=0.01)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel('Column', fontsize=11)
        ax.set_ylabel('Row', fontsize=11)

    fig.suptitle(
        f'Flowcell Channel Map — Pass-Filter Reads  ({ncols} cols × {nrows} rows)',
        fontsize=16,
    )
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Q-score over time
# ---------------------------------------------------------------------------

def plot_qscore_over_time(df, outpath, bin_minutes=10):
    """Mean Q-score (± 1 SD shaded) over run time, binned in <bin_minutes>-minute windows."""
    df2 = df.dropna(subset=['start_time']).copy()
    df2['time_bin_hr'] = (df2['start_time'] / 60 // bin_minutes * bin_minutes) / 60
    grp = (
        df2.groupby('time_bin_hr')['mean_qscore_template']
        .agg(['mean', 'std'])
        .reset_index()
    )
    grp.columns = ['time_hr', 'mean_q', 'std_q']
    grp['std_q'] = grp['std_q'].fillna(0)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(grp['time_hr'], grp['mean_q'], color=PASS_COLOR, linewidth=1.5)
    ax.fill_between(
        grp['time_hr'],
        (grp['mean_q'] - grp['std_q']).clip(lower=0),
        grp['mean_q'] + grp['std_q'],
        alpha=0.2, color=PASS_COLOR,
    )
    ax.set_title(f'Mean Q-Score over Time ({bin_minutes}-min bins)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Mean Q-Score', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Active pores over time
# ---------------------------------------------------------------------------

def plot_active_pores_over_time(df, outpath):
    """Bar chart of unique active channels per hour of run time."""
    df2 = df.dropna(subset=['start_time']).copy()
    df2['hour'] = (df2['start_time'] / 3600).astype(int)
    active = df2.groupby('hour')['channel'].nunique().reset_index()
    active.columns = ['hour', 'active_pores']

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.bar(active['hour'], active['active_pores'],
           color=PASS_COLOR, alpha=0.85, width=0.8)
    ax.set_title('Active Pores per Hour', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Unique Active Channels', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length CDF (base-weighted)
# ---------------------------------------------------------------------------

def plot_read_length_cdf(df, n50, outpath, max_length=None, min_length=None):
    """Base-weighted CDF of read lengths with N50 crosshair marked."""
    pf = df[df['passes_filtering'] == True].sort_values('sequence_length_template').copy()
    mean_len = pf['sequence_length_template'].mean()
    sd_len   = pf['sequence_length_template'].std()
    if max_length is None:
        max_length = int(mean_len + 3 * sd_len)
    if min_length is None:
        min_length = 0

    # Normalise to ALL pass-filter bases so that CDF(N50) == 0.5 by definition.
    # Filtering by length before computing total_bases would inflate CDF values
    # and shift the N50 crosshair off the curve.
    total_bases = pf['sequence_length_template'].sum()
    pf['cdf'] = pf['sequence_length_template'].cumsum() / total_bases

    # Restrict display to [min_length, max_length]; the curve may not start at 0
    # or reach 1.0 if reads exist outside the range — this is intentional.
    pf_display = pf[
        (pf['sequence_length_template'] <= max_length) &
        (pf['sequence_length_template'] >= min_length)
    ]

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(pf_display['sequence_length_template'], pf_display['cdf'],
            color=PASS_COLOR, linewidth=1.5)
    ax.axvline(n50, color=N50_COLOR, linestyle='--', linewidth=1.5, alpha=0.8)
    ax.axhline(0.5, color=N50_COLOR, linestyle='--', linewidth=1.5, alpha=0.8)
    ax.text(n50, 0.53, f' N50 = {n50:,} bp',
            color=N50_COLOR, fontsize=11, va='bottom')
    ax.set_xlim(min_length, max_length)
    ax.set_ylim(0, 1.02)
    ax.set_title('Read Length CDF — Fraction of Total Bases (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Read Length (bp)', fontsize=13)
    ax.set_ylabel('Fraction of Total Bases', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Q-score tier stacked area over time
# ---------------------------------------------------------------------------

def plot_qscore_bins_over_time(df, outpath):
    """Stacked area of read proportions in Q-score tiers (<10, 10–15, 15–20, ≥20) per hour."""
    pf = df[df['passes_filtering'] == True].dropna(subset=['start_time']).copy()
    pf['hour'] = (pf['start_time'] / 3600).astype(int)

    bins   = [0, 10, 15, 20, float('inf')]
    labels = ['Q < 10', 'Q10–15', 'Q15–20', 'Q ≥ 20']
    pf['q_tier'] = pd.cut(
        pf['mean_qscore_template'], bins=bins, labels=labels, right=False
    )

    pivot = (
        pf.groupby(['hour', 'q_tier'], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for lbl in labels:
        if lbl not in pivot.columns:
            pivot[lbl] = 0
    pivot = pivot[labels]
    pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    pal = sns.color_palette('RdYlGn', len(labels))
    pivot_pct.plot(kind='area', stacked=True, ax=ax, color=pal, alpha=0.85, linewidth=0)
    ax.set_ylim(0, 100)
    ax.set_title('Q-Score Tiers over Time (Pass-Filter Reads)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Proportion of Reads (%)', fontsize=13)
    ax.legend(title='Q-Score Tier', fontsize=11, loc='lower left')
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length vs run time (hexbin)
# ---------------------------------------------------------------------------

def plot_read_length_vs_time(df, outpath, max_length=None, min_length=None):
    """Hexbin density plot of read length vs run time (pass-filter reads)."""
    pf = df[df['passes_filtering'] == True].copy()
    mean_len = pf['sequence_length_template'].mean()
    sd_len   = pf['sequence_length_template'].std()
    if max_length is None:
        max_length = int(mean_len + 3 * sd_len)
    if min_length is None:
        min_length = 0
    pf = pf[
        (pf['sequence_length_template'] <= max_length) &
        (pf['sequence_length_template'] >= min_length)
    ]

    fig, ax = plt.subplots(figsize=(12, 7))
    hb = ax.hexbin(
        pf['start_time'] / 3600,
        pf['sequence_length_template'],
        gridsize=60,
        cmap='Blues',
        mincnt=1,
    )
    plt.colorbar(hb, ax=ax, label='Read Count')
    ax.set_title('Read Length vs. Run Time (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Read Length (bp)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Cumulative N50 over time
# ---------------------------------------------------------------------------

def plot_cumulative_n50_over_time(df, outpath):
    """
    Line plot of N50 recomputed cumulatively at each hour of run time.
    Shows whether fragment length improves, degrades, or stabilises as the
    run progresses.
    """
    pf = (
        df[df['passes_filtering'] == True]
        .dropna(subset=['start_time', 'sequence_length_template'])
        .sort_values('start_time')[['start_time', 'sequence_length_template']]
        .copy()
        .reset_index(drop=True)
    )
    pf['hour'] = (pf['start_time'] / 3600).astype(int)
    hours = sorted(pf['hour'].unique())

    hours_out, n50_out = [], []
    cumulative_lengths = []
    for h in hours:
        new = pf.loc[pf['hour'] == h, 'sequence_length_template'].values
        cumulative_lengths.extend(new)
        arr = np.sort(cumulative_lengths)[::-1]
        cumsum = np.cumsum(arr)
        idx = np.searchsorted(cumsum, cumsum[-1] * 0.5)
        n50_out.append(int(arr[min(idx, len(arr) - 1)]))
        hours_out.append(h)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(hours_out, n50_out, color=N50_COLOR, linewidth=1.8)
    ax.fill_between(hours_out, n50_out, alpha=0.15, color=N50_COLOR)
    ax.set_title('Cumulative N50 over Run Time', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('N50 (bp)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Translocation speed over time
# ---------------------------------------------------------------------------

def plot_translocation_speed_over_time(df, outpath, bin_minutes=10):
    """
    Mean translocation speed (bp/s) ± 1 SD over run time, binned in
    <bin_minutes>-minute windows.  Pass-filter reads only.
    """
    pf = (
        df[df['passes_filtering'] == True]
        .dropna(subset=['start_time', 'sequencing_speed'])
        .copy()
    )
    pf = pf[pf['sequencing_speed'] > 0]
    pf['time_bin_hr'] = (pf['start_time'] / 60 // bin_minutes * bin_minutes) / 60

    grp = (
        pf.groupby('time_bin_hr')['sequencing_speed']
        .agg(['mean', 'std'])
        .reset_index()
    )
    grp.columns = ['time_hr', 'mean_speed', 'std_speed']
    grp['std_speed'] = grp['std_speed'].fillna(0)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(grp['time_hr'], grp['mean_speed'], color=PASS_COLOR, linewidth=1.5)
    ax.fill_between(
        grp['time_hr'],
        (grp['mean_speed'] - grp['std_speed']).clip(lower=0),
        grp['mean_speed'] + grp['std_speed'],
        alpha=0.2, color=PASS_COLOR,
    )
    ax.set_title(f'Translocation Speed over Time ({bin_minutes}-min bins)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Translocation Speed (bp/s)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Pore survival curve
# ---------------------------------------------------------------------------

def plot_pore_survival(df, outpath):
    """
    Pore survival: fraction of channels active at T=0 still producing reads
    each hour.  Channels active in the first 10 min define the initial pool.
    """
    df2 = df.dropna(subset=['start_time', 'channel']).copy()
    initial_channels = set(df2[df2['start_time'] <= 600]['channel'].unique())
    if not initial_channels:
        return None

    df2['hour'] = (df2['start_time'] / 3600).astype(int)
    hours = sorted(df2['hour'].unique())

    hrs_out, survival = [], []
    for h in hours:
        active = set(df2[df2['hour'] == h]['channel'].unique())
        survival.append(len(initial_channels & active) / len(initial_channels) * 100)
        hrs_out.append(h)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(hrs_out, survival, color=PASS_COLOR, linewidth=2)
    ax.fill_between(hrs_out, survival, alpha=0.15, color=PASS_COLOR)
    ax.set_ylim(0, 105)
    ax.set_title('Pore Survival Curve', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Pores Still Active (%)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length distribution by Q-score tier (violin)
# ---------------------------------------------------------------------------

def plot_length_by_qscore_tier(df, outpath, max_length=None, min_length=None):
    """
    Violin plot of pass-filter read length split into Q-score tiers.
    Reveals whether long reads have systematically lower quality.
    """
    pf = df[df['passes_filtering'] == True].copy()
    if pf.empty:
        return None

    mean_len = pf['sequence_length_template'].mean()
    sd_len   = pf['sequence_length_template'].std()
    if max_length is None:
        max_length = int(mean_len + 3 * sd_len)
    pf = pf[pf['sequence_length_template'] <= max_length]
    if min_length is not None:
        pf = pf[pf['sequence_length_template'] >= min_length]

    bins   = [0, 10, 15, 20, float('inf')]
    labels = ['Q7–10', 'Q10–15', 'Q15–20', 'Q≥20']
    pf['q_tier'] = pd.cut(
        pf['mean_qscore_template'], bins=bins, labels=labels, right=False
    )

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    sns.violinplot(
        data=pf,
        x='q_tier', y='sequence_length_template',
        hue='q_tier',
        order=labels,
        hue_order=labels,
        palette='RdYlGn',
        inner='box',
        cut=0,
        legend=False,
        ax=ax,
    )
    ax.set_title('Read Length Distribution by Q-Score Tier (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Q-Score Tier', fontsize=13)
    ax.set_ylabel('Read Length (bp)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Yield by flowcell region / MUX group
# ---------------------------------------------------------------------------

def plot_yield_by_mux_group(df, outpath):
    """
    Bar chart of pass-filter yield per flowcell region.
      MinION / GridION  (≤512 ch)  → 4 groups of 128 channels each
      PromethION        (>512 ch)  → 12 blocks of 250 channels each
    Uneven bars reveal loading or pore-quality issues in specific regions.
    """
    pf = df[df['passes_filtering'] == True].copy()
    if pf.empty:
        return None

    max_channel = int(df['channel'].max())
    if max_channel <= 512:
        group_size = 128
        label_prefix = 'Group'
    else:
        group_size = 250
        label_prefix = 'Block'

    pf['region'] = label_prefix + ' ' + (
        ((pf['channel'] - 1) // group_size + 1).astype(int).astype(str)
    )

    grp = (
        pf.groupby('region')['sequence_length_template']
        .sum() / 1e9
    ).reset_index()
    grp.columns = ['region', 'yield_gb']

    # Sort numerically by the integer at the end of the label (e.g. "Group 2" < "Group 10")
    grp['_sort_key'] = grp['region'].str.extract(r'(\d+)$').astype(int)
    grp = grp.sort_values('_sort_key').reset_index(drop=True)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.bar(grp['region'], grp['yield_gb'], color=PASS_COLOR, alpha=0.85)
    ax.set_title('Pass-Filter Yield by Flowcell Region', fontsize=16)
    ax.set_xlabel('Region', fontsize=13)
    ax.set_ylabel('Yield (Gb)', fontsize=13)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# End-reason breakdown
# ---------------------------------------------------------------------------

def plot_end_reason(df, outpath):
    """
    Horizontal bar chart of read end reasons.
    Returns None silently if the end_reason column is absent (old files).
    """
    if 'end_reason' not in df.columns:
        return None

    counts = df['end_reason'].value_counts().sort_values()
    pal = sns.color_palette('muted', len(counts))

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    bars = ax.barh(counts.index, counts.values, color=pal[::-1], alpha=0.85)

    # Annotate each bar with the count
    for bar, val in zip(bars, counts.values):
        ax.text(
            bar.get_width() + counts.values.max() * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f'{val:,}',
            va='center', ha='left', fontsize=10,
        )

    ax.set_title('Read End Reason', fontsize=16)
    ax.set_xlabel('Read Count', fontsize=13)
    ax.set_ylabel('')
    ax.margins(x=0.15)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Animated channel strand-time video
# ---------------------------------------------------------------------------

def plot_channel_strand_video(df, outpath, fps=1, bin_hours=1):
    """
    Animated heatmap of per-channel strand time (seconds/bin) over the run.

    Each frame covers `bin_hours` of sequencing time and is played at `fps`
    frames per second, so by default (fps=1, bin_hours=1) one hour of run
    time equals one second of video — a 72-hour run → 72-second video.

    Output is MP4 (requires ffmpeg). Falls back to GIF automatically if
    ffmpeg is not available.

    Parameters
    ----------
    df         : sequencing summary DataFrame (must have start_time, duration, channel)
    outpath    : output path, e.g. 'run_strand_video.mp4'
    fps        : frames per second in the output video (default 1)
    bin_hours  : hours of run time collapsed into each frame (default 1)
    """
    import matplotlib.animation as animation

    pf = df.dropna(subset=['start_time', 'duration', 'channel']).copy()
    pf['hour_bin'] = (pf['start_time'] / 3600 // bin_hours).astype(int)

    max_channel = int(pf['channel'].max())
    nrows, ncols = _detect_flowcell_layout(max_channel)

    # Total strand time per channel per hour bin
    grp = (
        pf.groupby(['hour_bin', 'channel'])['duration']
        .sum()
        .reset_index(name='strand_time')
    )
    pivot = grp.pivot(index='hour_bin', columns='channel', values='strand_time')
    hours = sorted(pivot.index.tolist())

    # Pre-compute channel → (row, col) once
    ch_to_rc = {}
    for ch in [int(c) for c in pivot.columns]:
        if max_channel <= 512:
            r, c = (ch - 1) // ncols, (ch - 1) % ncols
        else:
            block = (ch - 1) // 250
            rem   = (ch - 1) % 250
            r     = rem // 10
            c     = rem % 10 + block * 10
        if 0 <= r < nrows and 0 <= c < ncols:
            ch_to_rc[ch] = (r, c)

    def build_frame(hour):
        grid = np.full((nrows, ncols), np.nan)
        if hour in pivot.index:
            for ch, val in pivot.loc[hour].items():
                ch = int(ch)
                if ch in ch_to_rc and not np.isnan(val):
                    r, c = ch_to_rc[ch]
                    grid[r, c] = val
        return grid

    # Robust colour scale capped at 99th percentile
    vmax = float(grp['strand_time'].quantile(0.99))

    # Figure sizing — same cell-inch rules as the static channel map
    cell_inch = 0.30 if max_channel <= 512 else 0.18
    panel_w   = ncols * cell_inch
    panel_h   = panel_w * (nrows / ncols)
    fig_w     = panel_w + 2.2
    fig_h     = panel_h + 1.8          # single panel + title/labels/colorbar

    with sns.axes_style("white"):
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    white_green = matplotlib.colors.LinearSegmentedColormap.from_list(
        'white_green', ['white', '#008000']
    )

    first = np.ma.masked_invalid(build_frame(hours[0]))
    im = ax.imshow(first, cmap=white_green, aspect='equal', origin='upper',
                   interpolation='nearest', vmin=0, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.01)
    cbar.set_label('Strand Time (s)', fontsize=11)
    h0 = hours[0] * bin_hours
    title_obj = ax.set_title(
        f'Channel Strand Time  [Hour {h0}–{h0 + bin_hours}]', fontsize=14
    )
    ax.set_xlabel('Column', fontsize=11)
    ax.set_ylabel('Row', fontsize=11)

    def update(i):
        hour    = hours[i]
        grid    = np.ma.masked_invalid(build_frame(hour))
        im.set_data(grid)
        h_start = hour * bin_hours
        title_obj.set_text(
            f'Channel Strand Time  [Hour {h_start}–{h_start + bin_hours}]'
        )
        return [im, title_obj]

    anim = animation.FuncAnimation(
        fig, update, frames=len(hours),
        interval=int(1000 / fps), blit=False,
    )

    # Try MP4 (H.264, widely compatible); fall back to GIF
    saved_path = outpath
    if outpath.endswith('.mp4'):
        try:
            writer = animation.FFMpegWriter(
                fps=fps, bitrate=2000,
                extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'],
            )
            anim.save(outpath, writer=writer, dpi=80)
        except Exception as exc:
            print(f'[ont_qc] ffmpeg unavailable ({exc}); saving as GIF instead.')
            saved_path = outpath.replace('.mp4', '.gif')
            anim.save(saved_path, writer='pillow', fps=fps)
    else:
        anim.save(outpath, writer='pillow', fps=fps)

    plt.close(fig)
    return saved_path
