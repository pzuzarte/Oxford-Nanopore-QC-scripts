#!/usr/bin/env python3
"""
ont_qc.py — Oxford Nanopore Sequencing QC Suite
================================================
Generates quality-control plots and a self-contained HTML report from the
output files of an Oxford Nanopore Technology (ONT) sequencing run.

The script reads one or more of the following MinKNOW output files:

  sequencing_summary*.txt   Per-read statistics (required). Contains read
                            lengths, Q-scores, channel assignments, timing,
                            and — for barcoded runs — barcode calls.
                            Barcoded runs are detected automatically.

  pore_activity*.csv        Per-minute channel state counts (optional).
                            Used to plot pore duty-time over the run.

  throughput_*.csv          Per-minute cumulative counters (optional).
                            Used to plot sequencing rate, read rate, pass
                            rate, and estimated vs basecalled base comparison.

All three files are auto-detected in the current directory if not specified
explicitly. Output is written to <runName>_qc/ by default.

Outputs
-------
  <runName>_summary.txt         Plain-text summary statistics
  <runName>_*.png               Individual QC plot images
  <runName>_report.html         Self-contained HTML report with all plots
                                and statistics embedded (no external deps)

Usage
-----
  # Minimal — auto-detect all input files:
  python ont_qc.py --runName MyRun

  # Explicit input files with subsampling (useful for very large runs):
  python ont_qc.py --file sequencing_summary.txt \\
                   --poreActivity pore_activity.csv \\
                   --throughput throughput_data.csv \\
                   --runName MyRun \\
                   --subsample 0.5

  # Barcoded run (barcode plots generated automatically):
  python ont_qc.py --file sequencing_summary_barcoded.txt --runName MyRun

  # Custom output directory and read-length axis cap:
  python ont_qc.py --runName MyRun --outdir /results/qc --maxLength 50000

Email Delivery
--------------
Add --email <address> to zip all outputs and send them as an attachment
once the report is complete. Three delivery methods are available via
--emailMethod:

  smtp (default)
    Standard SMTP. Supports open relays (no credentials) and authenticated
    servers. For authenticated servers pass --smtpUser and --smtpPass on
    the command line, or set the environment variables ONT_QC_SMTP_USER
    and ONT_QC_SMTP_PASS instead.

    Port guide:
      25   unauthenticated relay (default)
      587  STARTTLS (most institutional / Gmail SMTP servers)
      465  SSL/TLS

    Gmail note: use an App Password, not your account password.
    Generate one at: myaccount.google.com/apppasswords

    Example (Gmail):
      python ont_qc.py --runName MyRun --email you@gmail.com \\
        --smtpHost smtp.gmail.com --smtpPort 587 \\
        --smtpUser you@gmail.com --smtpPass <app-password>

    Example (institutional relay, no auth):
      python ont_qc.py --runName MyRun --email you@lab.org \\
        --smtpHost smtp.lab.org --smtpPort 25

  sendmail
    Pipes the message to the local sendmail/postfix binary. No credentials
    are needed — the system MTA handles delivery. Works on most Linux
    servers and macOS systems with postfix configured.

    Example:
      python ont_qc.py --runName MyRun --email you@lab.org \\
        --emailMethod sendmail

  gmail-oauth
    Sends via the Gmail API using OAuth2. No password is stored.
    On the first run a browser window opens for a one-time Google sign-in;
    the access token is cached at ~/.config/ont_qc/gmail_token.json and
    refreshed automatically on subsequent runs.

    Requirements:
      pip install google-auth-oauthlib google-api-python-client
    Download credentials.json from:
      Google Cloud Console -> APIs & Services -> Credentials
      -> Create Credentials -> OAuth client ID -> Desktop app

    Example:
      python ont_qc.py --runName MyRun --email you@gmail.com \\
        --emailMethod gmail-oauth \\
        --gmailCredentials ~/Downloads/credentials.json
"""

import argparse
import os
import sys
import warnings

# Seaborn ≤ 0.13 triggers a pandas FutureWarning about use_inf_as_na.
# This is an internal seaborn issue; suppress it to keep output clean.
warnings.filterwarnings(
    'ignore',
    message='use_inf_as_na option is deprecated',
    category=FutureWarning,
)

from qc_modules import loader, stats, seq_plots, barcode_plots, duty_plots, throughput_plots, report, emailer


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Generate QC plots and an HTML report from an ONT sequencing summary file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        '--file', type=str, default=None,
        help='Path to sequencing_summary*.txt (auto-detected if omitted).',
    )
    p.add_argument(
        '--runName', type=str, default='',
        help='Run name prefix used in output filenames and the HTML report.',
    )
    p.add_argument(
        '--outdir', type=str, default=None,
        help='Output directory (default: <runName>_qc/ or ont_qc_output/).',
    )
    p.add_argument(
        '--poreActivity', type=str, default=None,
        help='Path to pore_activity*.csv (auto-detected if omitted).',
    )
    p.add_argument(
        '--throughput', type=str, default=None,
        help='Path to throughput_*.csv (auto-detected if omitted).',
    )
    p.add_argument(
        '--subsample', type=float, default=1.0,
        help='Fraction of reads to use (0 < value <= 1). Useful for very large files.',
    )
    p.add_argument(
        '--maxLength', type=int, default=None,
        help='X-axis cap for read length plots (default: mean + 3×SD).',
    )
    p.add_argument(
        '--prop', type=int, default=20_000,
        help='Read length cutoff for proportion-above-cutoff metric (default: 20000).',
    )
    p.add_argument(
        '--video', action='store_true', default=False,
        help='Generate an animated MP4 (or GIF) showing per-channel strand time over the run. '
             '1 hour of sequencing = 1 second of video. Requires ffmpeg for MP4 output.',
    )

    # Email options
    email_group = p.add_argument_group(
        'email',
        'Send the zipped report by email. Three methods are available:\n'
        '  smtp         Standard SMTP (default). Pass --smtpUser / --smtpPass\n'
        '               for authenticated servers, or omit for open relays.\n'
        '  sendmail     Pipe to local sendmail/postfix binary. No credentials\n'
        '               needed if the system MTA is configured.\n'
        '  gmail-oauth  Gmail API with OAuth2. No password stored. Opens a\n'
        '               browser once for sign-in; token cached afterwards.\n'
        '               Requires --gmailCredentials and:\n'
        '               pip install google-auth-oauthlib google-api-python-client',
    )
    email_group.add_argument(
        '--email', type=str, default=None, metavar='ADDRESS',
        help='Recipient email address. Triggers sending a zipped report.',
    )
    email_group.add_argument(
        '--emailMethod', type=str, default='smtp',
        choices=['smtp', 'sendmail', 'gmail-oauth'],
        help='Delivery method: smtp | sendmail | gmail-oauth (default: smtp).',
    )
    email_group.add_argument(
        '--smtpHost', type=str, default='localhost',
        help='[smtp] SMTP server hostname (default: localhost).',
    )
    email_group.add_argument(
        '--smtpPort', type=int, default=25,
        help='[smtp] Port: 25 = unauthenticated, 587 = STARTTLS, 465 = SSL (default: 25).',
    )
    email_group.add_argument(
        '--smtpFrom', type=str, default='ont_qc@localhost',
        help='Sender address in the From: header (default: ont_qc@localhost).',
    )
    email_group.add_argument(
        '--smtpUser', type=str, default=None, metavar='USER',
        help='[smtp] SMTP username (overrides ONT_QC_SMTP_USER env var).',
    )
    email_group.add_argument(
        '--smtpPass', type=str, default=None, metavar='PASS',
        help='[smtp] SMTP password (overrides ONT_QC_SMTP_PASS env var).',
    )
    email_group.add_argument(
        '--gmailCredentials', type=str, default=None, metavar='FILE',
        help='[gmail-oauth] Path to credentials.json from Google Cloud Console.',
    )

    return p.parse_args()


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------

def make_outdir(args):
    if args.outdir:
        outdir = args.outdir
    elif args.runName:
        outdir = f'./{args.runName}_qc'
    else:
        outdir = './ont_qc_output'
    os.makedirs(outdir, exist_ok=True)
    return outdir


def p(outdir, run_name, suffix):
    """Build an output file path: <outdir>/<run_name>_<suffix>"""
    prefix = f'{run_name}_' if run_name else ''
    return os.path.join(outdir, f'{prefix}{suffix}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # --- Resolve input files ------------------------------------------------
    summary_file = args.file or loader.find_summary_file()
    # Search for companion files in the same directory as the summary file
    search_dir   = os.path.dirname(os.path.abspath(summary_file))
    pore_file       = args.poreActivity or loader.find_pore_activity_file(search_dir)
    throughput_file = args.throughput   or loader.find_throughput_file(search_dir)

    print(f'[ont_qc] Summary file   : {summary_file}')
    if pore_file:
        print(f'[ont_qc] Pore activity  : {pore_file}')
    if throughput_file:
        print(f'[ont_qc] Throughput CSV : {throughput_file}')

    # --- Load data ----------------------------------------------------------
    print('[ont_qc] Loading sequencing summary...')
    df, is_barcoded = loader.load_summary(summary_file, subsample=args.subsample)

    if is_barcoded:
        print('[ont_qc] Barcoded run detected.')
    else:
        print('[ont_qc] Non-barcoded run.')

    print(f'[ont_qc] Loaded {len(df):,} reads.')

    # --- Output directory ---------------------------------------------------
    outdir = make_outdir(args)
    run_name = args.runName
    print(f'[ont_qc] Output directory: {outdir}')

    # --- Statistics ---------------------------------------------------------
    print('[ont_qc] Computing statistics...')
    run_stats = stats.compute_stats(df, cutoff_length=args.prop)
    stats.write_summary_txt(run_stats, p(outdir, run_name, 'summary.txt'))
    print('[ont_qc] Summary statistics written.')

    n50 = run_stats['n50']
    max_length = args.maxLength  # may be None → functions will auto-compute

    # --- Core plots ---------------------------------------------------------
    plot_registry = {
        'throughput':   [],
        'run_health':   [],
        'read_quality': [],
        'channel_map':  [],
        'barcodes':     [],
        'duty_time':    [],
    }

    print('[ont_qc] Generating core plots...')

    path = seq_plots.plot_yield_over_time(df, p(outdir, run_name, 'yield_over_time.png'))
    plot_registry['throughput'].append((path, 'Pass-Filter Cumulative Yield vs. Time'))

    path = seq_plots.plot_reads_per_hour(df, p(outdir, run_name, 'reads_per_hour.png'))
    plot_registry['throughput'].append((path, 'Reads per Hour'))

    path = seq_plots.plot_active_pores_over_time(df, p(outdir, run_name, 'active_pores.png'))
    plot_registry['run_health'].append((path, 'Active Pores per Hour'))

    path = seq_plots.plot_qscore_over_time(df, p(outdir, run_name, 'qscore_over_time.png'))
    plot_registry['run_health'].append((path, 'Mean Q-Score over Time'))

    path = seq_plots.plot_cumulative_n50_over_time(df, p(outdir, run_name, 'cumulative_n50.png'))
    plot_registry['run_health'].append((path, 'Cumulative N50 over Time'))

    path = seq_plots.plot_read_length_dist(
        df, n50, p(outdir, run_name, 'read_length_dist.png'),
        max_length=max_length,
    )
    plot_registry['read_quality'].append((path, 'Read Length Distribution'))

    path = seq_plots.plot_read_length_cdf(
        df, n50, p(outdir, run_name, 'read_length_cdf.png'),
        max_length=max_length,
    )
    plot_registry['read_quality'].append((path, 'Read Length CDF'))

    path = seq_plots.plot_length_proportions(
        df, n50, p(outdir, run_name, 'length_proportions.png')
    )
    plot_registry['read_quality'].append((path, 'Proportion of Bases by Read Length'))

    path = seq_plots.plot_qscore_dist(df, p(outdir, run_name, 'qscore_dist.png'))
    plot_registry['read_quality'].append((path, 'Q-Score Distribution (Pass vs. Fail)'))

    path = seq_plots.plot_qscore_bins_over_time(df, p(outdir, run_name, 'qscore_bins_over_time.png'))
    plot_registry['read_quality'].append((path, 'Q-Score Tiers over Time'))

    path = seq_plots.plot_length_vs_qscore(
        df, p(outdir, run_name, 'length_vs_qscore.png')
    )
    plot_registry['read_quality'].append((path, 'Read Length vs. Q-Score'))

    path = seq_plots.plot_length_by_qscore_tier(
        df, p(outdir, run_name, 'length_by_qscore_tier.png'),
        max_length=max_length,
    )
    if path:
        plot_registry['read_quality'].append((path, 'Read Length by Q-Score Tier'))

    path = seq_plots.plot_end_reason(df, p(outdir, run_name, 'end_reason.png'))
    if path:
        plot_registry['read_quality'].append((path, 'Read End Reason'))

    path = seq_plots.plot_read_length_vs_time(
        df, p(outdir, run_name, 'read_length_vs_time.png'),
        max_length=max_length,
    )
    plot_registry['read_quality'].append((path, 'Read Length vs. Run Time'))

    path = seq_plots.plot_speed_over_time(df, p(outdir, run_name, 'speed_over_time.png'))
    plot_registry['read_quality'].append((path, 'Sequencing Speed over Time'))

    path = seq_plots.plot_translocation_speed_over_time(
        df, p(outdir, run_name, 'translocation_speed.png')
    )
    plot_registry['run_health'].append((path, 'Translocation Speed over Time'))

    path = seq_plots.plot_pore_survival(df, p(outdir, run_name, 'pore_survival.png'))
    if path:
        plot_registry['run_health'].append((path, 'Pore Survival Curve'))

    path = seq_plots.plot_median_current(df, p(outdir, run_name, 'median_current.png'))
    if path:
        plot_registry['read_quality'].append((path, 'Median Current over Time'))

    path = seq_plots.plot_channel_map(df, p(outdir, run_name, 'channel_map.png'))
    plot_registry['channel_map'].append((path, 'Flowcell Channel Map (Q-score / Reads / Yield)'))

    path = seq_plots.plot_yield_by_mux_group(df, p(outdir, run_name, 'yield_by_region.png'))
    if path:
        plot_registry['channel_map'].append((path, 'PF Yield by Flowcell Region'))

    # --- Barcode plots ------------------------------------------------------
    if is_barcoded:
        print('[ont_qc] Generating barcode plots...')

        path = barcode_plots.plot_barcode_read_counts(
            df, p(outdir, run_name, 'barcode_read_counts.png')
        )
        plot_registry['barcodes'].append((path, 'Read Counts per Barcode'))

        path = barcode_plots.plot_barcode_yield(
            df, p(outdir, run_name, 'barcode_yield.png')
        )
        plot_registry['barcodes'].append((path, 'Yield per Barcode'))

        path = barcode_plots.plot_barcode_qscore(
            df, p(outdir, run_name, 'barcode_qscore.png')
        )
        plot_registry['barcodes'].append((path, 'Q-Score Distribution per Barcode'))

        path = barcode_plots.plot_barcode_read_length(
            df, p(outdir, run_name, 'barcode_read_length.png')
        )
        plot_registry['barcodes'].append((path, 'Read Length Distribution per Barcode'))

        path = barcode_plots.plot_barcode_read_length_cdf(
            df, p(outdir, run_name, 'barcode_read_length_cdf.png')
        )
        plot_registry['barcodes'].append((path, 'Read Length CDF per Barcode'))

        path = barcode_plots.plot_barcode_accumulation(
            df, p(outdir, run_name, 'barcode_accumulation.png')
        )
        plot_registry['barcodes'].append((path, 'Cumulative Read Accumulation per Barcode'))

    # --- Pore activity plot -------------------------------------------------
    plotly_sections = {}
    if pore_file:
        print('[ont_qc] Generating pore activity plot...')
        path = duty_plots.plot_duty_time(pore_file, p(outdir, run_name, 'duty_time.png'))
        plot_registry['duty_time'].append((path, 'Pore Activity / Duty Time'))

        path = duty_plots.plot_occupancy_over_time(
            pore_file, p(outdir, run_name, 'occupancy_over_time.png')
        )
        plot_registry['duty_time'].append((path, 'Pore Occupancy over Time'))

        plotly_html = duty_plots.plot_duty_time_plotly(pore_file)
        if plotly_html:
            plotly_sections['duty_time'] = plotly_html

    # --- Throughput CSV plots -----------------------------------------------
    if throughput_file:
        print('[ont_qc] Generating throughput plots...')

        path = throughput_plots.plot_sequencing_rate(
            throughput_file, p(outdir, run_name, 'sequencing_rate.png')
        )
        plot_registry['run_health'].append((path, 'Sequencing Rate (Gb/hr)'))

        path = throughput_plots.plot_read_rate(
            throughput_file, p(outdir, run_name, 'read_rate.png')
        )
        plot_registry['run_health'].append((path, 'Read Acquisition Rate (reads/min)'))

        path = throughput_plots.plot_pass_rate(
            throughput_file, p(outdir, run_name, 'pass_rate_over_time.png')
        )
        plot_registry['run_health'].append((path, 'Rolling Pass Rate over Time'))

        path = throughput_plots.plot_estimated_vs_basecalled(
            throughput_file, p(outdir, run_name, 'estimated_vs_basecalled.png')
        )
        plot_registry['run_health'].append((path, 'Estimated vs. Basecalled Bases'))

        path = throughput_plots.plot_yield_per_hour(
            throughput_file, p(outdir, run_name, 'yield_per_hour.png')
        )
        plot_registry['throughput'].append((path, 'Yield per Hour'))

    # --- Strand-time video --------------------------------------------------
    vid_path = None
    if args.video:
        print('[ont_qc] Generating channel strand-time video (this may take a moment)...')
        vid_path = seq_plots.plot_channel_strand_video(
            df, p(outdir, run_name, 'strand_video.mp4'),
            fps=2, bin_hours=0.5,
        )
        print(f'[ont_qc] Video saved: {vid_path}')

    # --- HTML report --------------------------------------------------------
    print('[ont_qc] Building HTML report...')
    report_path = p(outdir, run_name, 'report.html')
    report.build_report(run_stats, plot_registry, report_path, run_name,
                        video_path=vid_path, plotly_sections=plotly_sections)
    print(f'[ont_qc] HTML report written: {report_path}')

    # --- Email --------------------------------------------------------------
    if args.email:
        print(f'[ont_qc] Zipping outputs for email...')
        zip_path = emailer.zip_outputs(outdir, run_name)
        emailer.send_email(
            recipient         = args.email,
            zip_path          = zip_path,
            run_name          = run_name,
            method            = args.emailMethod,
            smtp_host         = args.smtpHost,
            smtp_port         = args.smtpPort,
            smtp_from         = args.smtpFrom,
            smtp_user         = args.smtpUser,
            smtp_pass         = args.smtpPass,
            gmail_credentials = args.gmailCredentials,
        )

    print('[ont_qc] Done.')


if __name__ == '__main__':
    main()
