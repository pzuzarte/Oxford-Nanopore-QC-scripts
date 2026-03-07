"""
Pore activity / duty-time plot.

Reads a pore_activity*.csv file produced by MinKNOW and renders a stacked
bar chart showing channel state over experiment time.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

try:
    import plotly.graph_objects as go
    _PLOTLY_AVAILABLE = True
except ImportError:
    _PLOTLY_AVAILABLE = False

sns.set_theme(style="darkgrid")


# State order and colour palette (matches MinKNOW colours where possible)
STATE_COLORS = {
    'strand':                     '#2196F3',   # blue
    'pore':                       '#4CAF50',   # green
    'adapter':                    '#8BC34A',   # light green
    'unblocking':                 '#FF9800',   # orange
    'no_pore':                    '#9E9E9E',   # grey
    'multiple':                   '#9C27B0',   # purple
    'saturated':                  '#F44336',   # red
    'locked':                     '#795548',   # brown
    'unavailable':                '#607D8B',   # blue-grey
    'zero':                       '#212121',   # near-black
    'pending_manual_reset':       '#FF5722',   # deep orange
    'pending_mux_change':         '#FFC107',   # amber
    'unclassified':               '#E0E0E0',   # light grey
    'unclassified_following_reset': '#BDBDBD', # grey
    'unknown_negative':           '#B71C1C',   # dark red
    'unknown_positive':           '#1B5E20',   # dark green
    'disabled':                   '#000000',   # black
}

STATE_ORDER = list(STATE_COLORS.keys())


def plot_duty_time(filepath, outpath, sample_every=10):
    """
    Stacked bar chart of pore channel states over experiment time.

    Parameters
    ----------
    filepath : str
        Path to pore_activity*.csv file.
    outpath : str
        Destination PNG path.
    sample_every : int
        Take every Nth time point to reduce plot density.
    """
    df = pd.read_csv(filepath)

    # Pivot to wide format: index = time, columns = state, values = state time
    wide = pd.pivot_table(
        df,
        values='State Time (samples)',
        index='Experiment Time (minutes)',
        columns='Channel State',
        aggfunc='sum',
        fill_value=0,
    ).reset_index()

    # Downsample time axis
    wide = wide.iloc[::sample_every, :].copy()

    # Keep only states that are present in this file
    available_states = [s for s in STATE_ORDER if s in wide.columns]
    other_states = [c for c in wide.columns if c != 'Experiment Time (minutes)' and c not in STATE_ORDER]
    plot_states = available_states + other_states

    colors = [STATE_COLORS.get(s, '#CCCCCC') for s in plot_states]

    fig, ax = plt.subplots(figsize=(16, 6))

    bottom = pd.Series([0.0] * len(wide), index=wide.index)
    for state, color in zip(plot_states, colors):
        if state in wide.columns:
            ax.bar(
                wide['Experiment Time (minutes)'],
                wide[state],
                bottom=bottom,
                color=color,
                label=state,
                width=(wide['Experiment Time (minutes)'].diff().median() or 1),
                align='edge',
            )
            bottom = bottom + wide[state].values

    ax.set_title('Pore Activity / Duty Time', fontsize=16)
    ax.set_xlabel('Experiment Time (minutes)', fontsize=13)
    ax.set_ylabel('State Time (samples)', fontsize=13)
    ax.legend(
        loc='upper left',
        bbox_to_anchor=(1.01, 1),
        fontsize=9,
        frameon=True,
    )
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(fig)
    return outpath


def plot_occupancy_over_time(filepath, outpath, sample_every=10):
    """
    Line plot of pore occupancy over run time.

    Occupancy (%) = strand / (strand + pore + adapter + unblocking) × 100

    Only states where a pore is actively available for sequencing are included
    in the denominator, so occupancy reflects how productively available pores
    are being used.

    Parameters
    ----------
    filepath : str
        Path to pore_activity*.csv file.
    outpath : str
        Destination PNG path.
    sample_every : int
        Take every Nth time point to reduce density.
    """
    df = pd.read_csv(filepath)

    wide = pd.pivot_table(
        df,
        values='State Time (samples)',
        index='Experiment Time (minutes)',
        columns='Channel State',
        aggfunc='sum',
        fill_value=0,
    ).reset_index().sort_values('Experiment Time (minutes)')

    wide = wide.iloc[::sample_every, :].copy()

    numerator_states   = ['strand']
    denominator_states = ['strand', 'pore', 'adapter', 'unblocking']

    numerator   = sum(wide[s] for s in numerator_states   if s in wide.columns)
    denominator = sum(wide[s] for s in denominator_states if s in wide.columns)
    occupancy   = (numerator / denominator.replace(0, float('nan'))) * 100

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(wide['Experiment Time (minutes)'], occupancy,
            color='#2196F3', linewidth=1.5)
    ax.fill_between(wide['Experiment Time (minutes)'], occupancy,
                    alpha=0.15, color='#2196F3')
    ax.set_ylim(0, 100)
    ax.set_title('Pore Occupancy over Time', fontsize=16)
    ax.set_xlabel('Experiment Time (minutes)', fontsize=13)
    ax.set_ylabel('Occupancy (%)', fontsize=13)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(fig)
    return outpath


def plot_duty_time_plotly(filepath, sample_every=10):
    """
    Interactive Plotly stacked bar chart of pore channel states over time.

    Returns an HTML string (self-contained, plotly.js embedded) suitable for
    direct insertion into the QC report.  Returns None if plotly is not installed.

    Parameters
    ----------
    filepath : str
        Path to pore_activity*.csv file.
    sample_every : int
        Take every Nth time point to keep the chart responsive.
    """
    if not _PLOTLY_AVAILABLE:
        print('[ont_qc] plotly not installed — skipping interactive duty-time chart. '
              'Install with: pip install plotly')
        return None

    df = pd.read_csv(filepath)

    wide = pd.pivot_table(
        df,
        values='State Time (samples)',
        index='Experiment Time (minutes)',
        columns='Channel State',
        aggfunc='sum',
        fill_value=0,
    ).reset_index().sort_values('Experiment Time (minutes)')

    wide = wide.iloc[::sample_every, :].copy()

    available_states = [s for s in STATE_ORDER if s in wide.columns]
    other_states     = [c for c in wide.columns
                        if c != 'Experiment Time (minutes)' and c not in STATE_ORDER]
    plot_states = available_states + other_states

    fig = go.Figure()
    for state in plot_states:
        if state in wide.columns:
            fig.add_trace(go.Bar(
                x=wide['Experiment Time (minutes)'],
                y=wide[state],
                name=state,
                marker_color=STATE_COLORS.get(state, '#CCCCCC'),
                hovertemplate='%{x} min<br>' + state + ': %{y:,.0f}<extra></extra>',
            ))

    fig.update_layout(
        barmode='stack',
        title=dict(text='Pore Activity / Duty Time', font=dict(size=18)),
        xaxis=dict(title='Experiment Time (minutes)', tickfont=dict(size=12)),
        yaxis=dict(title='State Time (samples)',      tickfont=dict(size=12)),
        legend=dict(title='Channel State', x=1.01, y=1, bgcolor='rgba(0,0,0,0)'),
        height=520,
        margin=dict(l=60, r=180, t=60, b=60),
        plot_bgcolor='white',
        paper_bgcolor='white',
    )
    fig.update_xaxes(showgrid=True, gridcolor='#eeeeee')
    fig.update_yaxes(showgrid=True, gridcolor='#eeeeee')

    return fig.to_html(full_html=False, include_plotlyjs=True)
