"""
Barcode-specific QC plot functions.

Called only when the summary file contains a 'barcode_arrangement' column.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set_theme(style="darkgrid")

_FIGSIZE_BAR  = (12, 7)
_FIGSIZE_VIO  = (14, 7)


def _save(fig, path):
    fig.savefig(path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    return path


def _pf_barcodes(df, exclude_unclassified=True, barcodes=None):
    """Return pass-filter rows with a clean barcode column."""
    pf = df[df['passes_filtering'] == True].copy()
    if exclude_unclassified:
        pf = pf[~pf['barcode_arrangement'].str.lower().str.contains('unclassified', na=True)]
    pf = pf.dropna(subset=['barcode_arrangement'])
    if barcodes is not None:
        pf = pf[pf['barcode_arrangement'].isin(barcodes)]
    return pf


# ---------------------------------------------------------------------------
# Read counts per barcode
# ---------------------------------------------------------------------------

def plot_barcode_read_counts(df, outpath, exclude_unclassified=True, barcodes=None):
    """Horizontal bar chart: read counts per barcode (pass-filter)."""
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes)
    counts = (
        pf['barcode_arrangement']
        .value_counts()
        .reset_index()
    )
    counts.columns = ['barcode', 'read_count']
    counts = counts.sort_values('barcode')

    fig, ax = plt.subplots(figsize=_FIGSIZE_BAR)
    sns.barplot(data=counts, x='read_count', y='barcode', ax=ax, color=sns.color_palette()[0])
    ax.set_title('Read Counts per Barcode (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Read Count', fontsize=13)
    ax.set_ylabel('Barcode', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Yield per barcode
# ---------------------------------------------------------------------------

def plot_barcode_yield(df, outpath, exclude_unclassified=True, barcodes=None):
    """Horizontal bar chart: sequencing yield (Mb) per barcode (pass-filter)."""
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes)
    yield_df = (
        pf.groupby('barcode_arrangement')['sequence_length_template']
        .sum()
        .reset_index()
    )
    yield_df.columns = ['barcode', 'yield_mb']
    yield_df['yield_mb'] = yield_df['yield_mb'] / 1e6
    yield_df = yield_df.sort_values('barcode')

    fig, ax = plt.subplots(figsize=_FIGSIZE_BAR)
    sns.barplot(data=yield_df, x='yield_mb', y='barcode', ax=ax, color=sns.color_palette()[2])
    ax.set_title('Yield per Barcode (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Yield (Mb)', fontsize=13)
    ax.set_ylabel('Barcode', fontsize=13)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Q-score distribution per barcode
# ---------------------------------------------------------------------------

def plot_barcode_qscore(df, outpath, exclude_unclassified=True, barcodes=None):
    """Violin plot: Q-score distribution per barcode (pass-filter)."""
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes)

    # Order barcodes alphabetically
    order = sorted(pf['barcode_arrangement'].unique())

    fig, ax = plt.subplots(figsize=_FIGSIZE_VIO)
    sns.violinplot(
        data=pf,
        x='barcode_arrangement',
        y='mean_qscore_template',
        hue='barcode_arrangement',
        order=order,
        ax=ax,
        palette='Blues',
        legend=False,
        inner='box',
        density_norm='width',
    )
    ax.set_title('Q-Score Distribution per Barcode (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Barcode', fontsize=13)
    ax.set_ylabel('Mean Q-Score', fontsize=13)
    ax.tick_params(axis='x', rotation=45)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length distribution per barcode
# ---------------------------------------------------------------------------

def plot_barcode_read_length(df, outpath, exclude_unclassified=True, min_length=None, max_length=None, barcodes=None):
    """Violin plot: read length distribution per barcode (pass-filter)."""
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes)

    # Cap at mean + 3 SD to avoid extreme outliers, unless max_length is provided
    mean_len = pf['sequence_length_template'].mean()
    sd_len   = pf['sequence_length_template'].std()
    cap      = max_length if max_length is not None else mean_len + 3 * sd_len
    pf = pf[pf['sequence_length_template'] <= cap].copy()
    if min_length is not None:
        pf = pf[pf['sequence_length_template'] >= min_length]

    order = sorted(pf['barcode_arrangement'].unique())

    fig, ax = plt.subplots(figsize=_FIGSIZE_VIO)
    sns.violinplot(
        data=pf,
        x='barcode_arrangement',
        y='sequence_length_template',
        hue='barcode_arrangement',
        order=order,
        ax=ax,
        palette='Greens',
        legend=False,
        inner='box',
        density_norm='width',
    )
    ax.set_title('Read Length Distribution per Barcode (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Barcode', fontsize=13)
    ax.set_ylabel('Read Length (bp)', fontsize=13)
    ax.tick_params(axis='x', rotation=45)
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Read length CDF per barcode
# ---------------------------------------------------------------------------

def plot_barcode_read_length_cdf(df, outpath, exclude_unclassified=True, min_length=None, max_length=None, barcodes=None):
    """Base-weighted read length CDF — one line per barcode (pass-filter)."""
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes)
    mean_len = pf['sequence_length_template'].mean()
    sd_len   = pf['sequence_length_template'].std()
    max_len  = max_length if max_length is not None else int(mean_len + 3 * sd_len)
    min_len  = min_length if min_length is not None else 0
    pf = pf[
        (pf['sequence_length_template'] <= max_len) &
        (pf['sequence_length_template'] >= min_len)
    ]

    bc_names = sorted(pf['barcode_arrangement'].unique())
    pal = sns.color_palette('tab20', len(bc_names))

    fig, ax = plt.subplots(figsize=(12, 7))
    for bc, color in zip(bc_names, pal):
        sub = (
            pf[pf['barcode_arrangement'] == bc]
            .sort_values('sequence_length_template')
        )
        total = sub['sequence_length_template'].sum()
        if total == 0:
            continue
        cdf = sub['sequence_length_template'].cumsum() / total
        ax.plot(sub['sequence_length_template'], cdf,
                label=bc, color=color, linewidth=1.5, alpha=0.85)

    ax.axhline(0.5, color='grey', linestyle='--', linewidth=1, alpha=0.6)
    ax.set_xlim(min_len, max_len)
    ax.set_ylim(0, 1.02)
    ax.set_title('Read Length CDF per Barcode (Pass-Filter)', fontsize=16)
    ax.set_xlabel('Read Length (bp)', fontsize=13)
    ax.set_ylabel('Fraction of Total Bases', fontsize=13)
    ax.legend(title='Barcode', fontsize=9, bbox_to_anchor=(1.01, 1), loc='upper left')
    return _save(fig, outpath)


# ---------------------------------------------------------------------------
# Barcode accumulation over time
# ---------------------------------------------------------------------------

def plot_barcode_accumulation(df, outpath, exclude_unclassified=True, barcodes=None):
    """
    Stacked area chart of cumulative read counts per barcode over run time.
    Shows when each sample was sequenced and whether yields were consistent.
    """
    pf = _pf_barcodes(df, exclude_unclassified, barcodes=barcodes).dropna(subset=['start_time']).copy()
    pf['hour'] = (pf['start_time'] / 3600).astype(int)

    pivot = (
        pf.groupby(['hour', 'barcode_arrangement'])
        .size()
        .unstack(fill_value=0)
    )
    pivot = pivot[sorted(pivot.columns)]
    pivot_cum = pivot.cumsum()

    pal = sns.color_palette('tab20', len(pivot_cum.columns))
    fig, ax = plt.subplots(figsize=(13, 6))
    pivot_cum.plot(kind='area', stacked=True, ax=ax, color=pal, alpha=0.8, linewidth=0)
    ax.set_title('Cumulative Read Accumulation per Barcode', fontsize=16)
    ax.set_xlabel('Run Time (hr)', fontsize=13)
    ax.set_ylabel('Cumulative Read Count', fontsize=13)
    ax.legend(title='Barcode', fontsize=9, bbox_to_anchor=(1.01, 1), loc='upper left')
    return _save(fig, outpath)
