import pandas as pd


def compute_n50(pf_df):
    """
    Compute N50 from a pass-filter DataFrame sorted by length (desc).
    Returns None if the DataFrame is empty.
    """
    if pf_df.empty:
        return None
    sorted_df = pf_df.sort_values('sequence_length_template', ascending=False).copy()
    sorted_df['cumsum'] = sorted_df['sequence_length_template'].cumsum()
    total = sorted_df['sequence_length_template'].sum()
    row = sorted_df[sorted_df['cumsum'] <= total / 2]
    if row.empty:
        return int(sorted_df['sequence_length_template'].iloc[0])
    return int(row['sequence_length_template'].iloc[-1])


def compute_stats(df, cutoff_length=20000):
    """
    Compute summary statistics from a sequencing summary DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full dataframe (pass + fail reads).
    cutoff_length : int
        Read length threshold for proportion-above metric.

    Returns
    -------
    dict
    """
    pf = df[df['passes_filtering'] == True]

    total_reads = len(df)
    pf_reads = len(pf)
    pct_pf = round((pf_reads / total_reads * 100), 2) if total_reads > 0 else 0.0

    yield_gb = round(df['sequence_length_template'].sum() / 1e9, 3)
    pf_yield_gb = round(pf['sequence_length_template'].sum() / 1e9, 3)

    run_time_hr = round(df['start_time'].max() / 3600, 2)

    # Channels that produced a read in the first 10 minutes — proxy for
    # pores available for sequencing at the start of the run.
    pores_t0 = int(df[df['start_time'] <= 600]['channel'].nunique())

    n50 = compute_n50(pf)

    mean_read_length = round(pf['sequence_length_template'].mean(), 0) if pf_reads > 0 else 0
    sd_read_length = round(pf['sequence_length_template'].std(), 0) if pf_reads > 0 else 0
    max_read_length = int(pf['sequence_length_template'].max()) if pf_reads > 0 else 0

    mean_qscore = round(pf['mean_qscore_template'].mean(), 2) if pf_reads > 0 else 0

    pf_above = pf[pf['sequence_length_template'] >= cutoff_length]
    total_pf_bases = pf['sequence_length_template'].sum()
    bases_above = pf_above['sequence_length_template'].sum()
    prop_above = round(bases_above / total_pf_bases, 4) if total_pf_bases > 0 else 0.0

    top10 = (
        pf.sort_values('sequence_length_template', ascending=False)
        .head(10)[['read_id', 'sequence_length_template']]
        .reset_index(drop=True)
    )
    top10.index += 1  # 1-based rank

    return {
        'run_time_hr': run_time_hr,
        'pores_t0': pores_t0,
        'total_reads': total_reads,
        'pf_reads': pf_reads,
        'pct_pf': pct_pf,
        'yield_gb': yield_gb,
        'pf_yield_gb': pf_yield_gb,
        'n50': n50,
        'mean_read_length': int(mean_read_length),
        'sd_read_length': int(sd_read_length),
        'max_read_length': max_read_length,
        'mean_qscore': mean_qscore,
        'cutoff_length': cutoff_length,
        'prop_above_cutoff': prop_above,
        'top10_reads': top10,
    }


def write_summary_txt(stats, outpath):
    """Write a human-readable summary text file."""
    with open(outpath, 'w') as f:
        print(f"Run time:                         {stats['run_time_hr']} hrs", file=f)
        print(f"Pores available @ T=0:            {stats['pores_t0']:,}", file=f)
        print(f"Total reads:                      {stats['total_reads']:,}", file=f)
        print(f"Pass-filter reads:                {stats['pf_reads']:,}", file=f)
        print(f"% pass-filter:                    {stats['pct_pf']} %", file=f)
        print(f"Total yield:                      {stats['yield_gb']} Gb", file=f)
        print(f"Pass-filter yield:                {stats['pf_yield_gb']} Gb", file=f)
        print(f"N50:                              {stats['n50']:,} bp", file=f)
        print(f"Mean PF read length:              {stats['mean_read_length']:,} bp", file=f)
        print(f"SD PF read length:                {stats['sd_read_length']:,} bp", file=f)
        print(f"Longest PF read:                  {stats['max_read_length']:,} bp", file=f)
        print(f"Mean PF Q-score:                  {stats['mean_qscore']}", file=f)
        print(
            f"Proportion PF bases >= {stats['cutoff_length']:,} bp:  {stats['prop_above_cutoff']}",
            file=f,
        )
        print("\nTop 10 longest PF reads:", file=f)
        print(stats['top10_reads'].to_string(index=True), file=f)
