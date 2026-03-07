#!/usr/bin/env python3
"""
plot_guide.py — ONT QC Suite Plot Reference Guide
===================================================
Generates a self-contained PDF describing every plot produced by ont_qc.py:
what it shows, how to read it, and what to look for.

Usage
-----
  python plot_guide.py                      # writes ONT_QC_Plot_Guide.pdf
  python plot_guide.py --out my_guide.pdf   # custom output path

Requirements
------------
  pip install reportlab
"""

import argparse
import sys

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, HRFlowable,
        PageBreak, ListFlowable, ListItem, Table, TableStyle,
        KeepTogether,
    )
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
except ImportError:
    print(
        '\n[plot_guide] reportlab is required to generate the PDF guide.\n'
        '  Install it with:  pip install reportlab\n'
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Colour palette (matches the HTML report header)
# ---------------------------------------------------------------------------
BLUE_DARK  = colors.HexColor('#1565C0')
BLUE_LIGHT = colors.HexColor('#E3F0FF')
GREY_RULE  = colors.HexColor('#CCCCCC')
ORANGE     = colors.HexColor('#E65100')


# ---------------------------------------------------------------------------
# Custom paragraph styles
# ---------------------------------------------------------------------------
def _build_styles():
    base = getSampleStyleSheet()

    styles = {}

    styles['cover_title'] = ParagraphStyle(
        'cover_title',
        parent=base['Title'],
        fontSize=30,
        leading=36,
        textColor=BLUE_DARK,
        spaceAfter=6,
        alignment=TA_CENTER,
    )
    styles['cover_sub'] = ParagraphStyle(
        'cover_sub',
        parent=base['Normal'],
        fontSize=13,
        textColor=colors.HexColor('#555555'),
        alignment=TA_CENTER,
        spaceAfter=4,
    )
    styles['section_heading'] = ParagraphStyle(
        'section_heading',
        parent=base['Heading1'],
        fontSize=16,
        leading=20,
        textColor=colors.white,
        backColor=BLUE_DARK,
        borderPad=(6, 8, 6, 8),
        spaceBefore=18,
        spaceAfter=4,
    )
    styles['plot_heading'] = ParagraphStyle(
        'plot_heading',
        parent=base['Heading2'],
        fontSize=12,
        leading=15,
        textColor=BLUE_DARK,
        spaceBefore=12,
        spaceAfter=3,
        leftIndent=0,
    )
    styles['body'] = ParagraphStyle(
        'body',
        parent=base['Normal'],
        fontSize=10,
        leading=14,
        alignment=TA_JUSTIFY,
        spaceAfter=4,
    )
    styles['bullet'] = ParagraphStyle(
        'bullet',
        parent=base['Normal'],
        fontSize=10,
        leading=14,
        leftIndent=14,
        firstLineIndent=-10,
        spaceAfter=2,
    )
    styles['label'] = ParagraphStyle(
        'label',
        parent=base['Normal'],
        fontSize=10,
        leading=13,
        textColor=ORANGE,
        fontName='Helvetica-Bold',
        spaceAfter=1,
    )
    styles['note'] = ParagraphStyle(
        'note',
        parent=base['Normal'],
        fontSize=9,
        leading=12,
        textColor=colors.HexColor('#666666'),
        leftIndent=12,
        spaceAfter=6,
    )
    styles['toc_section'] = ParagraphStyle(
        'toc_section',
        parent=base['Normal'],
        fontSize=11,
        leading=16,
        textColor=BLUE_DARK,
        fontName='Helvetica-Bold',
        spaceBefore=6,
    )
    styles['toc_item'] = ParagraphStyle(
        'toc_item',
        parent=base['Normal'],
        fontSize=10,
        leading=14,
        leftIndent=16,
    )

    return styles


# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------
def _section(title, styles):
    """Full-width dark-blue section banner."""
    return [
        Spacer(1, 0.3 * cm),
        Paragraph(f'  {title}', styles['section_heading']),
        Spacer(1, 0.1 * cm),
    ]


def _plot_block(title, data_source, description, look_for, interpret, styles,
                optional=False):
    """
    Returns a KeepTogether block for one plot entry.

    Parameters
    ----------
    title       : str  — plot name (as it appears in the report)
    data_source : str  — which input file it comes from
    description : str  — one-paragraph description of what the plot shows
    look_for    : list[str] — bullet points: what to examine
    interpret   : list[str] — bullet points: how to interpret findings
    styles      : dict
    optional    : bool — show "Optional" badge if True
    """
    items = []

    badge = ' <font color="#E65100"><b>[optional]</b></font>' if optional else ''
    items.append(Paragraph(f'{title}{badge}', styles['plot_heading']))

    items.append(Paragraph(
        f'<font color="#666666"><i>Source: {data_source}</i></font>',
        styles['note'],
    ))

    items.append(Paragraph(description, styles['body']))

    items.append(Paragraph('What to look for:', styles['label']))
    for pt in look_for:
        items.append(Paragraph(f'• {pt}', styles['bullet']))

    items.append(Spacer(1, 0.15 * cm))
    items.append(Paragraph('How to interpret:', styles['label']))
    for pt in interpret:
        items.append(Paragraph(f'• {pt}', styles['bullet']))

    items.append(HRFlowable(width='100%', thickness=0.5,
                             color=GREY_RULE, spaceAfter=4))

    return KeepTogether(items)


# ---------------------------------------------------------------------------
# Plot catalogue
# ---------------------------------------------------------------------------
def _build_plots(styles):
    """Return a flat list of ReportLab flowables for the full plot catalogue."""
    story = []

    # ================================================================
    # SECTION 1 — THROUGHPUT
    # ================================================================
    story += _section('1  Throughput', styles)

    story.append(_plot_block(
        title='Pass-Filter Cumulative Yield vs. Time',
        data_source='sequencing_summary*.txt',
        description=(
            'A scatter plot of pass-filter reads plotted cumulatively over run time. '
            'Each point represents a read; its x-position is the time it was sequenced '
            '(hours) and its y-position is the total yield accumulated up to that point '
            '(gigabases). The overall shape of the curve reflects how productively the '
            'flowcell performed throughout the run.'
        ),
        look_for=[
            'Steepness of the early slope — a steep initial rise indicates a well-loaded, '
            'active flowcell.',
            'Changes in slope — shoulders or plateaux indicate periods of reduced output.',
            'Whether the curve levels off well before the run end.',
        ],
        interpret=[
            'A smooth S-curve (fast start, gradual plateau) is normal: pores deplete over time.',
            'A very shallow early slope suggests poor loading, buffer issues, or a '
            'cold start — check active pores at T=0.',
            'Sudden flat sections followed by recovery are often MUX scan events where '
            'channels are temporarily taken offline.',
            'A curve that plateaus very early and never recovers may indicate a failed '
            'flowcell or a library prep issue.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Reads per Hour',
        data_source='sequencing_summary*.txt',
        description=(
            'A histogram binned by clock hour showing how many reads were produced '
            'in each hour of the run. Unlike the cumulative yield curve, this '
            'non-cumulative view makes productivity stalls immediately visible as '
            'short or absent bars.'
        ),
        look_for=[
            'Hours with dramatically lower read counts — these are stalls.',
            'Whether read rate declines smoothly or drops abruptly.',
            'Any hours with zero or near-zero reads.',
        ],
        interpret=[
            'A gradual decline from left to right is expected as pores deplete.',
            'Sharp drops followed by partial recovery are usually MUX scans.',
            'A single empty hour early in the run may indicate the sequencer '
            'paused for a MUX scan or pore equilibration.',
            'Sustained low output after a drop (no recovery) suggests pore death, '
            'a blockage event, or the flowcell running out of active pores.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Yield per Hour',
        data_source='throughput_*.csv',
        description=(
            'A bar chart of total basecalled yield (Gb) produced in each clock hour, '
            'derived from the throughput CSV. This complements "Reads per Hour" by '
            'accounting for read length — an hour with few but very long reads will '
            'still show high yield here.'
        ),
        look_for=[
            'Consistency of yield across hours.',
            'Hours where yield drops without a corresponding drop in read count '
            '(indicates shorter reads).',
        ],
        interpret=[
            'Yield per hour naturally declines as pores deplete — this is normal.',
            'A sudden drop in yield without a matching drop in read count means '
            'reads became shorter — check translocation speed and N50 over time.',
            'Very high early yield that falls rapidly may indicate a large DNA '
            'fragment library being sequenced through quickly.',
        ],
        styles=styles,
        optional=True,
    ))

    # ================================================================
    # SECTION 2 — RUN HEALTH OVER TIME
    # ================================================================
    story += _section('2  Run Health over Time', styles)

    story.append(_plot_block(
        title='Active Pores per Hour',
        data_source='sequencing_summary*.txt',
        description=(
            'A bar chart showing the number of unique channels that produced at least '
            'one read in each clock hour. This is the primary indicator of flowcell '
            'health: as pores deplete, die, or become blocked, this number declines.'
        ),
        look_for=[
            'The count at hour 0 — this is your effective flowcell loading.',
            'Rate of decline — slow and smooth vs. sudden drops.',
            'Any partial recoveries after drops (MUX scan fingerprint).',
        ],
        interpret=[
            'A high starting count (close to the flowcell maximum: 512 for MinION, '
            '2675 for PromethION) indicates good loading.',
            'Starting counts below ~50% of the flowcell maximum suggest poor loading, '
            'an aged flowcell, or a library too dilute or too concentrated.',
            'Step-like drops then partial recoveries are MUX scans switching to '
            'fresh pore wells — this is normal and healthy.',
            'A monotonic, steep decline without recoveries suggests irreversible '
            'pore fouling (common with dirty libraries or aged flowcells).',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Mean Q-Score over Time',
        data_source='sequencing_summary*.txt',
        description=(
            'A line plot of mean Q-score per 10-minute bin (all reads, pass + fail) '
            'with ±1 standard deviation shaded. Tracks whether basecall accuracy '
            'changes as the run progresses.'
        ),
        look_for=[
            'Whether Q-score is stable, improving, or declining over time.',
            'Width of the shaded band — narrow = consistent quality.',
            'Early vs. late run quality difference.',
        ],
        interpret=[
            'A stable flat line at or above your threshold (Q10 for kit14, Q7 for '
            'older kits) is ideal.',
            'A gradual decline is normal and expected as pores age; a steep early '
            'decline suggests reagent or temperature issues.',
            'Sudden improvement mid-run can follow a MUX scan loading fresh pores '
            'with cleaner signal.',
            'Very wide SD bands indicate heterogeneous quality — common with '
            'mixed-length or degraded libraries.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Cumulative N50 over Time',
        data_source='sequencing_summary*.txt',
        description=(
            'The N50 is recomputed cumulatively at the end of each clock hour and '
            'plotted as a line. N50 is the length such that 50% of all sequenced '
            'bases come from reads at least this long. Watching it evolve over time '
            'reveals whether fragment length is improving, degrading, or stable.'
        ),
        look_for=[
            'Whether N50 rises, falls, or stabilises as the run proceeds.',
            'How quickly N50 stabilises — small runs may show large early fluctuations.',
        ],
        interpret=[
            'N50 typically rises in the first few hours as more reads accumulate '
            'and then plateaus once the dataset is statistically stable.',
            'A steadily falling N50 late in the run means shorter fragments are '
            'being sequenced — long molecules may have been consumed or '
            'pore-blocking events are truncating reads.',
            'A very low N50 that never rises suggests a degraded library; '
            'consider running a new prep.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Translocation Speed over Time',
        data_source='sequencing_summary*.txt',
        description=(
            'Mean translocation speed (bp/s) ± 1 SD per 10-minute bin, restricted '
            'to pass-filter reads. Translocation speed reflects how fast DNA moves '
            'through the nanopore and is sensitive to temperature, buffer conditions, '
            'and motor protein function.'
        ),
        look_for=[
            'Baseline speed (typically 400–500 bp/s for R10.4 Kit14 at 400 bps mode).',
            'Drift in speed over time — upward or downward.',
            'Sudden shifts that correlate with changes in other metrics.',
        ],
        interpret=[
            'Speed is set by the motor protein and should be relatively stable for '
            'a given kit — minor drift of ±10% is normal.',
            'A gradual speed increase often indicates a warming sequencer or '
            'evaporation concentrating the buffer.',
            'Sudden speed drops may indicate the sequencer is struggling to keep '
            'up with data throughput (adaptive sampling backpressure) or a '
            'temperature drop.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Pore Survival Curve',
        data_source='sequencing_summary*.txt',
        description=(
            'Channels that produced a read in the first 10 minutes define the '
            'starting pool (100%). Each subsequent hour shows what fraction of '
            'those original pores are still active. This is the classic ONT pore '
            'survival metric — directly analogous to a cell survival curve in biology.'
        ),
        look_for=[
            'Half-life of the pore population — the hour at which ~50% survive.',
            'Shape of the curve: exponential decay vs. abrupt steps vs. plateau.',
            'Final survival fraction at run end.',
        ],
        interpret=[
            'Smooth exponential decline is normal; individual pores burn out '
            'stochastically over the run.',
            'Step-like drops aligned with MUX scans are expected — those pores '
            'were switched to new wells, so they may reappear briefly.',
            'A very steep early drop (>50% loss in the first few hours) suggests '
            'pore fouling from dirty DNA, excessive adapter, or a library that '
            'is too concentrated.',
            'A flat curve with >80% survival through the whole run indicates '
            'minimal pore loss — great for long-read sequencing.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Sequencing Rate (Gb/hr)',
        data_source='throughput_*.csv',
        description=(
            'Instantaneous output rate in Gb/hr derived by differencing adjacent '
            'rows in the throughput CSV. Shows the moment-to-moment productivity '
            'of the sequencer rather than cumulative totals.'
        ),
        look_for=[
            'Peak rate — reached typically in the first 1–2 hours.',
            'Shape of decline — smooth or irregular.',
            'Spikes or troughs that correlate with events in other plots.',
        ],
        interpret=[
            'Rate declines naturally as pores deplete; this is expected.',
            'Short spikes followed by drops are often MUX scans: a brief pause '
            'then burst of output from freshly selected pores.',
            'Extended flat sections near zero before run end suggest the flowcell '
            'should have been stopped earlier.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Read Acquisition Rate (reads/min)',
        data_source='throughput_*.csv',
        description=(
            'Number of reads captured per minute, derived from the throughput CSV. '
            'Complements yield-based metrics: read rate is high when many short reads '
            'are captured, even if yield per read is low.'
        ),
        look_for=[
            'Sustained high rates indicate good pore activity.',
            'Discordance between read rate and yield rate (could indicate short reads).',
        ],
        interpret=[
            'A high read rate with low yield suggests many short reads — '
            'check N50 and read length distribution.',
            'A very low read rate despite many active pores may indicate '
            'the library is at low concentration and pores are idle between reads.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Rolling Pass Rate over Time',
        data_source='throughput_*.csv',
        description=(
            '30-minute rolling average of the fraction of reads that pass the '
            'quality filter (Q-score threshold), expressed as a percentage. '
            'A dashed reference line is drawn at 90%.'
        ),
        look_for=[
            'Whether pass rate stays above the reference line.',
            'Trends: rising, falling, or stable.',
            'Whether pass rate correlates with changes in active pore count.',
        ],
        interpret=[
            'Pass rates consistently above 90% indicate high-quality signal; '
            'below 80% is a concern for most applications.',
            'A declining pass rate over time may indicate degrading signal quality '
            'from ageing pores.',
            'A very low pass rate early in the run (first 30 min) often normalises '
            'once pores equilibrate — watch for recovery.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Estimated vs. Basecalled Bases (Cumulative)',
        data_source='throughput_*.csv',
        description=(
            'Two cumulative lines plotted together: estimated bases (from raw signal '
            'length) and basecalled bases (reads that have been processed by the '
            'basecaller). When real-time basecalling keeps up with acquisition, the '
            'lines overlap. A growing gap indicates the basecaller is lagging behind.'
        ),
        look_for=[
            'How closely the two lines track each other.',
            'Whether any gap closes later in the run (basecaller catching up).',
        ],
        interpret=[
            'Lines that overlap throughout indicate real-time basecalling running '
            'at full speed — ideal.',
            'A widening gap means the basecaller cannot keep up with acquisition; '
            'the remaining reads will be basecalled post-run.',
            'A gap that closes near the end of a run suggests the basecaller '
            'caught up after sequencing slowed down.',
        ],
        styles=styles,
        optional=True,
    ))

    # ================================================================
    # SECTION 3 — READ QUALITY
    # ================================================================
    story += _section('3  Read Quality', styles)

    story.append(_plot_block(
        title='Read Length Distribution',
        data_source='sequencing_summary*.txt',
        description=(
            'A histogram of read lengths (all reads). The x-axis is capped at '
            'mean + 3 standard deviations to prevent extreme outliers from '
            'compressing the bulk of the distribution. A vertical dashed red line '
            'marks the N50.'
        ),
        look_for=[
            'Shape: unimodal (one prep), bimodal (two size populations), or '
            'heavy right tail (ultra-long reads present).',
            'Position of the N50 relative to the distribution peak.',
            'Proportion of very short reads (< 1 kb) — often adapter dimers.',
        ],
        interpret=[
            'A roughly log-normal distribution peaked around 5–15 kb is typical '
            'for standard high-molecular-weight gDNA.',
            'A large spike at very short lengths (<500 bp) may indicate carry-over '
            'of short DNA fragments, adapter dimers, or incomplete adapter ligation.',
            'A bimodal distribution may indicate mixing of two library pools or '
            'incomplete size selection.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Read Length CDF (Base-Weighted)',
        data_source='sequencing_summary*.txt',
        description=(
            'The cumulative distribution function of read lengths weighted by base '
            'count. The y-axis shows the fraction of total sequenced bases that '
            'come from reads up to a given length. N50 crosshairs are drawn: the '
            'curve crosses y=0.5 at exactly the N50 length by definition.'
        ),
        look_for=[
            'Steepness: a steep curve means most bases are in long reads.',
            'Whether the curve reaches 1.0 within the displayed range '
            '(if not, very long reads extend beyond the axis cap).',
        ],
        interpret=[
            'A curve that rises slowly and reaches 0.5 at a high length value '
            'indicates a long-read-enriched library — desirable for most long-read applications.',
            'A curve that reaches 0.5 at a low length value (low N50) suggests '
            'a short-fragment-dominated library.',
            'The upper tail of the curve shows what fraction of total bases '
            'come from ultra-long reads.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Proportion of Bases by Read Length',
        data_source='sequencing_summary*.txt',
        description=(
            'A scatter plot of read length (x) vs. cumulative proportion of total '
            'pass-filter bases (y), with reads sorted longest-first. '
            'This is equivalent to the inverse-sorted CDF and directly answers: '
            '"What fraction of my total yield comes from reads longer than X?"'
        ),
        look_for=[
            'The x-value at y=0.5 is the N50.',
            'How flat the curve is at large x — a long flat tail means a useful '
            'population of ultra-long reads.',
        ],
        interpret=[
            'If the curve reaches 0.9 (90% of bases) at a short read length, '
            'the library is short-fragment dominated.',
            'A curve with a long, slowly-declining tail beyond the N50 indicates '
            'a broad, well-distributed size range — typical of good HMW DNA extraction.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Q-Score Distribution (Pass vs. Fail)',
        data_source='sequencing_summary*.txt',
        description=(
            'Kernel density estimate (KDE) curves of per-read mean Q-score, '
            'with pass-filter reads in blue and fail reads in red, both overlaid. '
            'The pass/fail split is defined by the Q-score threshold set in MinKNOW '
            '(Q10 for Kit14, Q7 for older kits).'
        ),
        look_for=[
            'Separation between pass and fail peaks.',
            'Width of the pass-filter peak — narrow is good (consistent quality).',
            'Whether any pass-filter reads have surprisingly low Q-scores.',
        ],
        interpret=[
            'Well-separated, narrow peaks indicate clean, high-quality data.',
            'A broad pass peak suggests heterogeneous library quality — '
            'possibly degraded DNA or a mixture of fragment sizes.',
            'A large fail peak overlapping the pass peak near the threshold '
            'is common and expected; reads near the threshold are borderline.',
            'For kit14 data, the pass peak should be centred around Q15–Q20.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Q-Score Tiers over Time',
        data_source='sequencing_summary*.txt',
        description=(
            'Stacked area chart showing the proportion of pass-filter reads '
            'in each of four Q-score tiers (Q7–10, Q10–15, Q15–20, Q≥20) per '
            'clock hour. Colours grade from red (low Q) to green (high Q). '
            'Each column always sums to 100%.'
        ),
        look_for=[
            'Whether the green (high Q) fraction increases or decreases over time.',
            'Any sudden shifts in the tier composition mid-run.',
        ],
        interpret=[
            'A stable tier distribution throughout the run indicates '
            'consistent basecall accuracy.',
            'An increasing red fraction late in the run suggests pore ageing '
            'is degrading signal quality.',
            'A sudden shift to higher Q-score tiers after a MUX scan reflects '
            'fresh, clean pores producing better signal.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Read Length vs. Q-Score (Hex Joint Plot)',
        data_source='sequencing_summary*.txt',
        description=(
            'A 2-D hexbin density plot of read length (x) against mean Q-score (y) '
            'for a random sample of up to 20,000 pass-filter reads. Hex colour '
            'encodes read density. Marginal distributions are shown on each axis.'
        ),
        look_for=[
            'The location of the densest hexbin cluster — this is the modal read.',
            'Correlation between length and Q-score: does quality drop for longer reads?',
            'Any unexpected clusters (e.g., a separate population of low-Q, long reads).',
        ],
        interpret=[
            'A slight negative correlation between length and Q-score is normal: '
            'very long reads have more accumulated pore noise.',
            'A strong negative correlation may indicate the motor protein is '
            'struggling with very long molecules.',
            'Two distinct clusters may suggest a mixed library (e.g., '
            'two PCR products of different sizes).',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Read Length Distribution by Q-Score Tier (Violin Plot)',
        data_source='sequencing_summary*.txt',
        description=(
            'Violin plots — one per Q-score tier — showing the full distribution '
            'of pass-filter read lengths within each quality band. An inner '
            'box-and-whisker plot marks the median and IQR. Colours grade '
            'red-to-green from low to high Q.'
        ),
        look_for=[
            'Whether the median read length decreases as Q-score tier increases.',
            'Width of the violin — wide = many reads; narrow = few.',
            'Any tier that is entirely absent or extremely narrow.',
        ],
        interpret=[
            'A consistent read length across all tiers means quality does not '
            'depend on fragment length for this library.',
            'A clear pattern where high-Q reads are short indicates the library '
            'has excellent short-read basecalling but long reads are noisier.',
            'An empty or near-empty Q≥20 tier is expected for R9.4.1 data; '
            'for R10.4 Kit14 data this tier should be well populated.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Read End Reason',
        data_source='sequencing_summary*.txt',
        description=(
            'A horizontal bar chart showing the count of reads grouped by the '
            'reason the sequencer stopped reading them. Each bar is labelled '
            'with its exact count. Common end reasons are: signal_positive '
            '(read completed normally), signal_negative (pore blocked), '
            'mux_change (taken off by MUX scan), unblock_mux_change '
            '(rejected by adaptive sampling), and data_service_unblock_mux_change '
            '(rejected by a data service rule). This plot is only generated when '
            'the end_reason column is present.'
        ),
        look_for=[
            'Which end reason dominates.',
            'The proportion of unblock_mux_change events — these are adaptive '
            'sampling rejections.',
            'The proportion of signal_negative events — these are blockages.',
        ],
        interpret=[
            'signal_positive should be the most common end reason for a normal run.',
            'High unblock_mux_change is expected and desirable on adaptive '
            'sampling runs; on non-adaptive runs it may indicate the library '
            'contains many short fragments that are being forcibly ejected.',
            'High signal_negative (blockages) can indicate dirty DNA, excess '
            'adapter, or pores that are becoming clogged.',
            'High mux_change counts correlate with MUX scan events — normal.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Read Length vs. Run Time (Hexbin)',
        data_source='sequencing_summary*.txt',
        description=(
            'A 2-D hexbin density plot showing how read length (y) is distributed '
            'across run time (x), restricted to pass-filter reads. Reveals whether '
            'the fragment length composition changes as the run progresses.'
        ),
        look_for=[
            'Whether the length distribution shifts over time.',
            'Diagonal banding — often indicates a correlated relationship.',
            'Any time windows where the distribution collapses to very short reads.',
        ],
        interpret=[
            'A stable, vertically symmetric blob throughout the run indicates '
            'no time-dependent length bias.',
            'A gradual compression toward shorter reads late in the run is common: '
            'long fragments are consumed first by the most active pores.',
            'Very short reads late in the run with no long reads may indicate '
            'that long fragments have been fully sequenced and only short '
            'DNA remnants remain.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Sequencing Speed over Time (Scatter)',
        data_source='sequencing_summary*.txt',
        description=(
            'All reads (pass and fail) plotted as points with start time on x '
            'and sequencing speed (bp/s) on y. Pass-filter reads are blue, '
            'fail reads are red. A random subsample of up to 50,000 reads '
            'is shown to keep the plot tractable for large runs.'
        ),
        look_for=[
            'The modal speed band — is it where you expect for your kit?',
            'Outlier reads with extreme high or low speed.',
            'Any time-dependent drift in the speed distribution.',
        ],
        interpret=[
            'R10.4 Kit14 400 bps mode → expected speed ~400 bp/s.',
            'R9.4.1 450 bps mode → expected speed ~450 bp/s.',
            'Very fast reads (>1000 bp/s) are often noise or very short fragments '
            'that zip through the pore without the motor protein slowing them down.',
            'Very slow reads (<100 bp/s) may indicate a stalling motor protein '
            'or a difficult DNA secondary structure.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Median Current over Time',
        data_source='sequencing_summary*.txt',
        description=(
            'Scatter plot of per-read median ion current (pA) vs. start time for '
            'all reads (pass and fail). The median template current reflects the '
            'electrical environment of the nanopore. This plot is only generated '
            'when the median_template column is present (older files).'
        ),
        look_for=[
            'Drift in the median current over the run.',
            'Whether fail reads cluster at systematically different current values.',
        ],
        interpret=[
            'A stable median current around 60–80 pA is typical for most kits '
            'and flowcells.',
            'A gradual upward drift may indicate buffer evaporation concentrating '
            'the salt solution and increasing conductivity.',
            'A large spread between pass and fail median currents suggests the '
            'current threshold is important for quality classification.',
        ],
        styles=styles,
        optional=True,
    ))

    # ================================================================
    # SECTION 4 — FLOWCELL CHANNEL MAP
    # ================================================================
    story += _section('4  Flowcell Channel Map', styles)

    story.append(_plot_block(
        title='Flowcell Channel Map (3-Panel Heatmap)',
        data_source='sequencing_summary*.txt',
        description=(
            'Three stacked heatmaps, each showing a property per physical channel '
            'on the flowcell membrane. The grid matches the MinKNOW channel layout: '
            '16×32 for MinION/GridION (512 channels) or 25×120 for PromethION '
            '(3,000 channels). Only pass-filter reads contribute. '
            'Panel 1: mean Q-score. Panel 2: read count. Panel 3: yield (Mb).'
        ),
        look_for=[
            'Spatial uniformity — should be roughly even across the flowcell.',
            'Any dark (dead) patches — rows, columns, or quadrants with no activity.',
            'Rows or columns with unusually high or low Q-score relative to neighbours.',
        ],
        interpret=[
            'Uniform colour across all three panels indicates even loading and '
            'uniform membrane quality — the ideal outcome.',
            'An entire row of dark channels often indicates a bad connection or '
            'a failed electrode group.',
            'A quadrant of dead channels corresponds to one of the four MUX groups '
            'failing — check the Pore Activity plot for that group.',
            'High read count but low yield in a region indicates short reads '
            'from that part of the membrane.',
            'Low Q-score in one corner may indicate uneven temperature distribution '
            'or physical contamination on the membrane surface.',
        ],
        styles=styles,
    ))

    story.append(_plot_block(
        title='Pass-Filter Yield by Flowcell Region',
        data_source='sequencing_summary*.txt',
        description=(
            'A bar chart comparing total pass-filter yield (Gb) across flowcell '
            'regions. For MinION/GridION, channels are divided into 4 groups of '
            '128 channels each. For PromethION, 12 blocks of 250 channels. '
            'Bars are ordered numerically by region number and coloured blue.'
        ),
        look_for=[
            'Uniformity of bar heights — large disparities indicate regional problems.',
            'One or more regions with zero or near-zero yield.',
        ],
        interpret=[
            'Roughly equal bars across all regions indicate uniform flowcell '
            'loading and performance.',
            'One region significantly below others may indicate that MUX group '
            'had fewer active pores, possibly due to a loading artefact or '
            'a section of membrane that was damaged or air-exposed.',
            'For PromethION, large inter-block variation is common due to the '
            'physical distance between the blocks on the chip.',
        ],
        styles=styles,
    ))

    # ================================================================
    # SECTION 5 — BARCODE ANALYSIS
    # ================================================================
    story += _section('5  Barcode Analysis  [barcoded runs only]', styles)

    story.append(Paragraph(
        'These plots are generated automatically when the sequencing summary file '
        'contains a barcode_arrangement column. Unclassified reads are excluded '
        'from all barcode plots by default.',
        styles['body'],
    ))
    story.append(Spacer(1, 0.2 * cm))

    story.append(_plot_block(
        title='Read Counts per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Horizontal bar chart of pass-filter read counts for each classified '
            'barcode, sorted alphabetically. Provides a quick visual balance '
            'check across samples.'
        ),
        look_for=[
            'Balance between barcodes — are any samples dramatically under- or over-represented?',
            'Any barcodes with zero reads.',
        ],
        interpret=[
            'For balanced multiplexing, bars should be roughly equal height.',
            'A factor of 2–3× variation is common and acceptable; 10× variation '
            'indicates a pipetting error or very different DNA concentrations '
            'in the pooled library.',
            'Zero reads for a barcode means that sample was not detected at all — '
            'check the library preparation for that sample.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Yield per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Horizontal bar chart of total pass-filter yield (Mb) per barcode. '
            'Combines read count and read length — a barcode with fewer but '
            'longer reads may have the same yield as one with many short reads.'
        ),
        look_for=[
            'Which barcodes achieved coverage targets.',
            'Discordance between read count and yield bars for the same barcode.',
        ],
        interpret=[
            'A barcode with low read count but normal yield has longer reads — '
            'check the read length violin plot.',
            'A barcode with high read count but low yield has unusually short reads '
            '— check for sample degradation or adapter contamination.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Q-Score Distribution per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Violin plot with one violin per barcode (width-normalised) showing '
            'the distribution of mean Q-scores for pass-filter reads. An inner '
            'box-and-whisker plot marks the median and IQR.'
        ),
        look_for=[
            'Any barcode with consistently lower Q-score than the others.',
            'Unusually wide violins indicating high Q-score variance in one sample.',
        ],
        interpret=[
            'Most barcodes should have similar Q-score distributions if '
            'the multiplexed libraries had similar fragment size and purity.',
            'A single barcode with low Q may indicate that sample had PCR '
            'inhibitors, damaged DNA, or a different GC content.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Read Length Distribution per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Violin plot with one violin per barcode showing the distribution '
            'of pass-filter read lengths. Lengths are capped at mean + 3 SD to '
            'prevent extreme outliers from collapsing all violins.'
        ),
        look_for=[
            'Whether all barcodes have similar size distributions.',
            'Any barcode with a very different modal length.',
        ],
        interpret=[
            'If barcodes were prepared from samples of similar size, violins '
            'should look nearly identical.',
            'A barcode with a much lower median length may represent a degraded '
            'DNA sample or a library with excess short-fragment contamination.',
            'A much longer median for one barcode may indicate that sample '
            'had a larger insert size or was not size-selected.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Read Length CDF per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Base-weighted cumulative distribution function (CDF) of read length '
            'for each barcode — one coloured line per barcode. A dashed grey line '
            'marks y=0.5 (the N50 level). The x-axis is capped at mean + 3 SD.'
        ),
        look_for=[
            'How close together the lines are — separation indicates '
            'different size distributions.',
            'Where each line crosses y=0.5 — this is the per-barcode N50.',
        ],
        interpret=[
            'Tightly clustered lines indicate reproducible size selection across '
            'samples in the same pool.',
            'A line shifted far left (early y=0.5 crossing) indicates that '
            'barcode has shorter reads — check the violin plot and Q-score plot '
            'for corroborating evidence.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Cumulative Read Accumulation per Barcode',
        data_source='sequencing_summary*.txt (barcoded)',
        description=(
            'Stacked area chart of cumulative read count per barcode over run time '
            '(one colour per barcode). Shows when each sample was sequenced and '
            'whether yield accumulation was consistent across the run.'
        ),
        look_for=[
            'Whether all barcodes accumulate reads at similar rates.',
            'Any barcode that stops accumulating reads mid-run.',
            'Whether proportions stay stable or shift over the run.',
        ],
        interpret=[
            'Parallel, evenly-spaced accumulation lines indicate balanced '
            'sequencing of all samples throughout the run.',
            'A barcode that plateaus early relative to others may represent '
            'a sample with fewer molecules in the pool.',
            'Changing proportions over time may reflect different fragment length '
            'distributions competing for pore occupancy.',
        ],
        styles=styles,
        optional=True,
    ))

    # ================================================================
    # SECTION 6 — PORE ACTIVITY / DUTY TIME
    # ================================================================
    story += _section('6  Pore Activity / Duty Time  [requires pore_activity*.csv]', styles)

    story.append(Paragraph(
        'These plots require the pore_activity*.csv file produced by MinKNOW. '
        'This file records the count of channels in each state at each minute '
        'of the run.',
        styles['body'],
    ))
    story.append(Spacer(1, 0.2 * cm))

    story.append(_plot_block(
        title='Pore Activity / Duty Time (Interactive)',
        data_source='pore_activity*.csv',
        description=(
            'An interactive Plotly stacked bar chart showing how many channels '
            'are in each state at each minute of the run. States include: '
            'strand (actively sequencing), pore (open, available), adapter '
            '(reading adapter sequence), unblocking (ejecting a molecule), '
            'no_pore, locked, saturated, and others. Hover over any bar to '
            'see exact counts. Toggle states on/off by clicking the legend. '
            'A static PNG version is also generated as a fallback.'
        ),
        look_for=[
            'The "strand" state (blue) — this is time spent actively sequencing.',
            'The "no_pore" state (grey) — channels with no functional pore.',
            'The "locked" state (brown) — pores stuck and unable to sequence.',
            'Sudden increases in "saturated" (red) state.',
        ],
        interpret=[
            'High "strand" fraction relative to total channels = high duty time = '
            'efficient sequencing.',
            'Increasing "no_pore" over time is expected as pores deplete — '
            'these channels will not recover without a fresh flowcell.',
            'A surge in "locked" state may indicate DNA or protein contamination '
            'blocking the pore entrance.',
            'High "saturated" counts (positive current overload) can occur with '
            'very concentrated libraries or if the flowcell is too hot.',
            'Periodic spikes in "unblocking" state are MUX scan events ejecting '
            'molecules from all active channels simultaneously.',
        ],
        styles=styles,
        optional=True,
    ))

    story.append(_plot_block(
        title='Pore Occupancy over Time',
        data_source='pore_activity*.csv',
        description=(
            'Line plot (with shaded area) of pore occupancy (%) over run time. '
            'Occupancy is defined as: strand ÷ (strand + pore + adapter + '
            'unblocking) × 100. The denominator counts only channels that have '
            'an active pore — so occupancy reflects how productively the '
            'available pores are being used, independent of how many pores '
            'have died.'
        ),
        look_for=[
            'Baseline occupancy — is it consistently high?',
            'Whether occupancy declines, rises, or fluctuates over the run.',
            'Periodic dips corresponding to MUX scans.',
        ],
        interpret=[
            'Occupancy above 70% indicates a library at good concentration: '
            'available pores spend most of their time on strands.',
            'Occupancy below 30% suggests the library is too dilute — pores '
            'spend most of their time open and waiting for a molecule to arrive.',
            'Occupancy above 95% consistently may indicate the library is too '
            'concentrated, causing strand queuing and premature stalls.',
            'Occupancy that declines over time may indicate library depletion '
            'as molecules are consumed.',
        ],
        styles=styles,
        optional=True,
    ))

    # ================================================================
    # SECTION 7 — SUMMARY TABLE
    # ================================================================
    story += _section('7  Top 10 Longest Pass-Filter Reads', styles)

    story.append(_plot_block(
        title='Top 10 Longest Pass-Filter Reads Table',
        data_source='sequencing_summary*.txt',
        description=(
            'A table listing the 10 longest pass-filter reads by read length '
            '(bp), including each read\'s unique identifier (read_id) and '
            'length in base pairs. Reads are ranked 1 (longest) to 10.'
        ),
        look_for=[
            'The length of the longest read — is it consistent with the expected '
            'maximum fragment size for your extraction method?',
            'Large gaps between consecutive entries in the list.',
        ],
        interpret=[
            'Read IDs can be used to extract specific reads from the FASTQ/BAM '
            'for manual inspection: samtools view -h file.bam <read_id>.',
            'A longest read of 500 kb or more indicates excellent ultra-long DNA '
            'extraction; above 1 Mb is exceptional.',
            'If the longest read is only 20–30 kb, the library lacks very long '
            'fragments — consider switching to an ultra-long extraction protocol.',
            'Very similar lengths among the top 10 (e.g., all around 200 kb) '
            'may indicate a gel plug extraction with consistent fragment size.',
        ],
        styles=styles,
    ))

    return story


# ---------------------------------------------------------------------------
# Cover page and table of contents
# ---------------------------------------------------------------------------
def _build_cover(styles):
    story = []
    story.append(Spacer(1, 3 * cm))
    story.append(Paragraph('ONT QC Suite', styles['cover_title']))
    story.append(Paragraph('Plot Reference Guide', styles['cover_title']))
    story.append(Spacer(1, 0.6 * cm))
    story.append(HRFlowable(width='60%', thickness=2, color=BLUE_DARK,
                             spaceAfter=12, hAlign='CENTER'))
    story.append(Paragraph(
        'A complete description of every plot generated by <b>ont_qc.py</b>.',
        styles['cover_sub'],
    ))
    story.append(Paragraph(
        'Includes what each visualisation shows, what to look for,',
        styles['cover_sub'],
    ))
    story.append(Paragraph(
        'and how to interpret normal vs. abnormal findings.',
        styles['cover_sub'],
    ))
    story.append(Spacer(1, 2 * cm))

    # Usage table
    usage_data = [
        ['Command', 'Description'],
        ['python ont_qc.py --runName MyRun',
         'Auto-detect all input files in the current directory'],
        ['python ont_qc.py --file summary.txt --runName MyRun',
         'Explicit summary file path'],
        ['... --poreActivity pore_activity.csv',
         'Explicit pore activity file (otherwise auto-detected)'],
        ['... --throughput throughput_data.csv',
         'Explicit throughput file (otherwise auto-detected)'],
        ['... --subsample 0.5',
         'Use 50% of reads — faster for very large files'],
        ['... --outdir /path/to/results',
         'Custom output directory (default: <runName>_qc/)'],
        ['... --maxLength 50000',
         'Cap read-length axis at 50 kb instead of mean + 3 SD'],
        ['... --prop 50000',
         'Change proportion-above-cutoff threshold (default: 20,000 bp)'],
        ['... --video',
         'Generate animated channel strand-time video (requires ffmpeg)'],
        ['python plot_guide.py',
         'Generate this PDF (requires: pip install reportlab)'],
    ]
    usage_tbl = Table(usage_data, colWidths=[8.0 * cm, 8.5 * cm])
    usage_tbl.setStyle(TableStyle([
        ('BACKGROUND',   (0, 0), (-1, 0), BLUE_DARK),
        ('TEXTCOLOR',    (0, 0), (-1, 0), colors.white),
        ('FONTNAME',     (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE',     (0, 0), (-1, -1), 8),
        ('FONTNAME',     (0, 1), (0, -1), 'Courier'),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [BLUE_LIGHT, colors.white]),
        ('GRID',         (0, 0), (-1, -1), 0.5, GREY_RULE),
        ('TOPPADDING',   (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ('LEFTPADDING',  (0, 0), (-1, -1), 7),
        ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('VALIGN',       (0, 0), (-1, -1), 'MIDDLE'),
    ]))
    story.append(usage_tbl)
    story.append(Spacer(1, 0.6 * cm))

    # Input file legend table
    legend_data = [
        ['Input File',            'Description',                                    'Required?'],
        ['sequencing_summary*.txt', 'Per-read TSV from Albacore, Guppy, or Dorado', 'Yes'],
        ['pore_activity*.csv',      'Per-minute channel state counts',              'Optional'],
        ['throughput_*.csv',        'Per-minute cumulative counters',               'Optional'],
    ]
    tbl = Table(legend_data, colWidths=[5.5 * cm, 8.5 * cm, 2.5 * cm])
    tbl.setStyle(TableStyle([
        ('BACKGROUND',   (0, 0), (-1, 0), BLUE_DARK),
        ('TEXTCOLOR',    (0, 0), (-1, 0), colors.white),
        ('FONTNAME',     (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE',     (0, 0), (-1, -1), 9),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [BLUE_LIGHT, colors.white]),
        ('GRID',         (0, 0), (-1, -1), 0.5, GREY_RULE),
        ('TOPPADDING',   (0, 0), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('LEFTPADDING',  (0, 0), (-1, -1), 8),
        ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ('VALIGN',       (0, 0), (-1, -1), 'MIDDLE'),
    ]))
    story.append(tbl)
    story.append(PageBreak())
    return story


def _build_toc(styles):
    story = []
    story.append(Paragraph('Contents', styles['section_heading']))
    story.append(Spacer(1, 0.3 * cm))

    toc = [
        ('1  Throughput', [
            'Pass-Filter Cumulative Yield vs. Time',
            'Reads per Hour',
            'Yield per Hour  [optional]',
        ]),
        ('2  Run Health over Time', [
            'Active Pores per Hour',
            'Mean Q-Score over Time',
            'Cumulative N50 over Time',
            'Translocation Speed over Time',
            'Pore Survival Curve',
            'Sequencing Rate (Gb/hr)  [optional]',
            'Read Acquisition Rate  [optional]',
            'Rolling Pass Rate over Time  [optional]',
            'Estimated vs. Basecalled Bases  [optional]',
        ]),
        ('3  Read Quality', [
            'Read Length Distribution',
            'Read Length CDF (Base-Weighted)',
            'Proportion of Bases by Read Length',
            'Q-Score Distribution (Pass vs. Fail)',
            'Q-Score Tiers over Time',
            'Read Length vs. Q-Score (Hex Joint Plot)',
            'Read Length Distribution by Q-Score Tier (Violin)',
            'Read End Reason',
            'Read Length vs. Run Time (Hexbin)',
            'Sequencing Speed over Time',
            'Median Current over Time  [optional]',
        ]),
        ('4  Flowcell Channel Map', [
            'Flowcell Channel Map (3-Panel Heatmap)',
            'Pass-Filter Yield by Flowcell Region',
        ]),
        ('5  Barcode Analysis  [barcoded runs only]', [
            'Read Counts per Barcode',
            'Yield per Barcode',
            'Q-Score Distribution per Barcode',
            'Read Length Distribution per Barcode',
            'Read Length CDF per Barcode',
            'Cumulative Read Accumulation per Barcode',
        ]),
        ('6  Pore Activity / Duty Time  [requires pore_activity*.csv]', [
            'Pore Activity / Duty Time (Interactive)',
            'Pore Occupancy over Time',
        ]),
        ('7  Top 10 Longest Pass-Filter Reads', [
            'Top 10 Longest Pass-Filter Reads Table',
        ]),
    ]

    for section_title, items in toc:
        story.append(Paragraph(section_title, styles['toc_section']))
        for item in items:
            story.append(Paragraph(f'  {item}', styles['toc_item']))
        story.append(Spacer(1, 0.1 * cm))

    story.append(PageBreak())
    return story


# ---------------------------------------------------------------------------
# Page numbering
# ---------------------------------------------------------------------------
def _add_page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(colors.HexColor('#999999'))
    canvas.drawRightString(
        doc.pagesize[0] - 1.5 * cm,
        1.0 * cm,
        f'Page {doc.page}',
    )
    canvas.drawString(
        1.5 * cm,
        1.0 * cm,
        'ONT QC Suite — Plot Reference Guide',
    )
    canvas.restoreState()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description='Generate the ONT QC Plot Reference Guide PDF.'
    )
    ap.add_argument(
        '--out', default='ONT_QC_Plot_Guide.pdf',
        help='Output PDF path (default: ONT_QC_Plot_Guide.pdf)',
    )
    args = ap.parse_args()

    styles = _build_styles()

    doc = SimpleDocTemplate(
        args.out,
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.8 * cm,
        bottomMargin=2.0 * cm,
        title='ONT QC Suite — Plot Reference Guide',
        author='ont_qc.py',
    )

    story = []
    story += _build_cover(styles)
    story += _build_toc(styles)
    story += _build_plots(styles)

    doc.build(story, onFirstPage=_add_page_number, onLaterPages=_add_page_number)
    print(f'[plot_guide] PDF written: {args.out}')


if __name__ == '__main__':
    main()
