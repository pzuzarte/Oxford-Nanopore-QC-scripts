"""
Throughput CSV plot functions.

Reads a throughput_*.csv file produced by MinKNOW (minute-resolution
cumulative counters) and generates run-health plots by differencing
adjacent rows to recover per-minute rates.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set_theme(style="darkgrid")

_PAL      = sns.color_palette()
_BLUE     = _PAL[0]
_ORANGE   = _PAL[1]
_GREEN    = _PAL[2]
_RED      = _PAL[3]

_FIGSIZE  = (13, 5)


def _save(fig, path):
    fig.savefig(path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    return path


def _load_rates(filepath):
    """
    Load throughput CSV and compute per-minute rate columns via differencing.
    Returns a DataFrame indexed by Experiment Time (hours).
    """
    df = pd.read_csv(filepath)
    df = df.sort_values('Experiment Time (minutes)').copy()
    df['time_hr'] = df['Experiment Time (minutes)'] / 60

    # Per-minute deltas (forward-fill NaN at row 0)
    df['reads_per_min']   = df['Reads'].diff().fillna(0).clip(lower=0)
    df['passed_per_min']  = df['Basecalled Reads Passed'].diff().fillna(0).clip(lower=0)
    df['failed_per_min']  = df['Basecalled Reads Failed'].diff().fillna(0).clip(lower=0)
    df['bases_per_min']   = df['Basecalled Bases'].diff().fillna(0).clip(lower=0)
    df['est_bases_per_min'] = df['Estimated Bases'].diff().fillna(0).clip(lower=0)

    # Gb/hr  =  bases/min × 60 / 1e9
    df['gbhr'] = df['bases_per_min'] * 60 / 1e9

    # Rolling pass rate (30-min window)
    total_per_min = (df['passed_per_min'] + df['failed_per_min']).replace(0, float('nan'))
    df['pass_rate'] = (df['passed_per_min'] / total_per_min * 100).rolling(30, min_periods=1).mean()

    return df


# ---------------------------------------------------------------------------
# Plot 1 — Sequencing rate (Gb/hr) over time
# ---------------------------------------------------------------------------

def plot_sequencing_rate(filepath, outpath):
    """Line plot of sequencing output rate (Gb/hr) vs run time."""
    df = _load_rates(filepath)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(df['time_hr'], df['gbhr'], color=_BLUE, linewidth=0.8, alpha=0.9)
    ax.fill_between(df['time_hr'], df['gbhr'], alpha=0.15, color=_BLUE)
    ax.set_title('Sequencing Rate over Time', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Output Rate (Gb/hr)', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Plot 2 — Read acquisition rate (reads/min) over time
# ---------------------------------------------------------------------------

def plot_read_rate(filepath, outpath):
    """Line plot of read acquisition rate (reads/min) vs run time."""
    df = _load_rates(filepath)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(df['time_hr'], df['reads_per_min'], color=_GREEN, linewidth=0.8, alpha=0.9)
    ax.fill_between(df['time_hr'], df['reads_per_min'], alpha=0.15, color=_GREEN)
    ax.set_title('Read Acquisition Rate over Time', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Reads per Minute', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Plot 3 — Rolling pass rate (%) over time
# ---------------------------------------------------------------------------

def plot_pass_rate(filepath, outpath):
    """Line plot of rolling 30-min pass rate (%) vs run time."""
    df = _load_rates(filepath)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(df['time_hr'], df['pass_rate'], color=_ORANGE, linewidth=1.2)
    ax.axhline(90, color=_RED, linestyle='--', linewidth=1, alpha=0.6, label='90 %')
    ax.set_ylim(0, 105)
    ax.set_title('Rolling Pass Rate over Time (30-min window)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Pass Rate (%)', fontsize=13)
    ax.legend(fontsize=11)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Plot 4 — Estimated vs Basecalled cumulative bases
# ---------------------------------------------------------------------------

def plot_estimated_vs_basecalled(filepath, outpath):
    """
    Dual-line plot comparing estimated bases (from signal) vs basecalled bases
    over run time. Divergence indicates skipped or backlogged reads.
    """
    df = pd.read_csv(filepath).sort_values('Experiment Time (minutes)').copy()
    df['time_hr'] = df['Experiment Time (minutes)'] / 60

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(df['time_hr'], df['Estimated Bases'] / 1e9,
            color=_ORANGE, linewidth=1.5, label='Estimated Bases')
    ax.plot(df['time_hr'], df['Basecalled Bases'] / 1e9,
            color=_BLUE, linewidth=1.5, label='Basecalled Bases')
    ax.set_title('Estimated vs. Basecalled Bases (Cumulative)', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Bases (Gb)', fontsize=13)
    ax.legend(fontsize=12)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Plot 5 — Yield per hour (bar chart)
# ---------------------------------------------------------------------------

def plot_yield_per_hour(filepath, outpath):
    """Bar chart of basecalled yield (Gb) produced in each clock hour of the run."""
    df = _load_rates(filepath)
    df['hour'] = df['time_hr'].astype(int)
    yield_per_hr = (
        df.groupby('hour')['bases_per_min']
        .sum()
        .reset_index()
    )
    yield_per_hr['yield_gb'] = yield_per_hr['bases_per_min'] / 1e9

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.bar(yield_per_hr['hour'], yield_per_hr['yield_gb'],
           color=_BLUE, alpha=0.85, width=0.8)
    ax.set_title('Basecalled Yield per Hour', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Yield (Gb)', fontsize=13)
    return _save(fig, outpath)
