"""
HTML report builder.

Generates a self-contained HTML file with embedded base64 images and a
summary statistics table. No external CSS/JS dependencies.
"""

import base64
import os
from datetime import datetime


# ---------------------------------------------------------------------------
# HTML template fragments
# ---------------------------------------------------------------------------

_PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ONT QC Report – {run_name}</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f5f5f5; color: #333; }}
    header {{ background: #1565C0; color: white; padding: 24px 32px; }}
    header h1 {{ font-size: 1.8em; }}
    header p  {{ font-size: 0.9em; opacity: 0.85; margin-top: 4px; }}
    main {{ max-width: 1400px; margin: 24px auto; padding: 0 24px 48px; }}
    section {{ margin-bottom: 40px; }}
    h2 {{ font-size: 1.3em; color: #1565C0; border-bottom: 2px solid #1565C0;
          padding-bottom: 6px; margin-bottom: 16px; }}
    .stats-table {{ border-collapse: collapse; width: 100%; max-width: 640px;
                    background: white; border-radius: 6px; overflow: hidden;
                    box-shadow: 0 1px 4px rgba(0,0,0,.12); }}
    .stats-table th, .stats-table td {{ padding: 10px 16px; text-align: left;
                                         font-size: 0.95em; }}
    .stats-table tr:nth-child(even) {{ background: #f0f4ff; }}
    .stats-table th {{ background: #1565C0; color: white; }}
    .plot-grid {{ display: flex; flex-wrap: wrap; gap: 20px; }}
    .plot-card {{ background: white; border-radius: 8px;
                  box-shadow: 0 1px 4px rgba(0,0,0,.12); overflow: hidden;
                  flex: 1 1 580px; max-width: 780px; }}
    .plot-card img {{ width: 100%; display: block; }}
    .plot-card p {{ padding: 8px 14px; font-size: 0.85em; color: #666; }}
    .video-card {{ background: white; border-radius: 8px;
                   box-shadow: 0 1px 4px rgba(0,0,0,.12); overflow: hidden;
                   width: 100%; }}
    .video-card video, .video-card img {{ width: 100%; display: block; }}
    .video-card p {{ padding: 8px 14px; font-size: 0.85em; color: #666; }}
    .plotly-card {{ background: white; border-radius: 8px;
                    box-shadow: 0 1px 4px rgba(0,0,0,.12); overflow: hidden;
                    width: 100%; padding: 8px; }}
    .reads-table {{ border-collapse: collapse; width: 100%; max-width: 860px;
                    background: white; border-radius: 6px; overflow: hidden;
                    box-shadow: 0 1px 4px rgba(0,0,0,.12); font-family: monospace; }}
    .reads-table th, .reads-table td {{ padding: 9px 16px; text-align: left; font-size: 0.9em; }}
    .reads-table tr:nth-child(even) {{ background: #f0f4ff; }}
    .reads-table th {{ background: #1565C0; color: white; font-family: 'Segoe UI', Arial, sans-serif; }}
    footer {{ text-align: center; font-size: 0.8em; color: #999; padding: 16px; }}
    .plot-caption {{ position: relative; display: flex; align-items: center; gap: 6px;
                     padding: 8px 14px; font-size: 0.85em; color: #666; }}
    .info-btn {{ display: inline-flex; align-items: center; justify-content: center;
                 width: 16px; height: 16px; border-radius: 50%;
                 background: #1565C0; color: white; font-size: 11px; font-weight: bold;
                 cursor: pointer; border: none; flex-shrink: 0; font-style: normal;
                 line-height: 1; padding: 0; }}
    .info-btn:hover {{ background: #0d47a1; }}
    .info-popup {{ position: fixed; display: none;
                   background: white; border: 1px solid #ccc; border-radius: 6px;
                   box-shadow: 0 4px 16px rgba(0,0,0,.25);
                   z-index: 1000; min-width: 440px; max-width: 700px;
                   font-size: 0.82em; color: #333;
                   font-family: 'Segoe UI', Arial, sans-serif; overflow: hidden; }}
    .info-popup.visible {{ display: block; }}
    .info-popup-header {{ display: flex; justify-content: space-between; align-items: center;
                          padding: 7px 10px; background: #1565C0; color: white;
                          cursor: grab; user-select: none; font-size: 1em; font-weight: 600; }}
    .info-popup-header:active {{ cursor: grabbing; }}
    .info-popup-body {{ padding: 12px 14px; max-height: 70vh; overflow-y: auto; }}
    .info-popup-close {{ background: none; border: none; color: white; font-size: 15px;
                         cursor: pointer; padding: 0 2px; line-height: 1; opacity: 0.85; }}
    .info-popup-close:hover {{ opacity: 1; }}
    .info-popup table {{ border-collapse: collapse; width: 100%; }}
    .info-popup th, .info-popup td {{ padding: 5px 8px; text-align: left;
                                      border-bottom: 1px solid #eee; font-size: 1em; }}
    .info-popup th {{ background: #f0f4ff; font-weight: 600; }}
    .info-popup td:first-child {{ font-family: monospace; white-space: nowrap;
                                  color: #1565C0; }}
  </style>
</head>
<body>
<header>
  <h1>ONT QC Report</h1>
  <p>Run: {run_name} &nbsp;|&nbsp; Generated: {generated}</p>
</header>
<main>
{sections}
</main>
<footer>Generated by ont_qc.py</footer>
<script>
document.addEventListener('click', function(e) {{
  if (e.target.classList.contains('info-btn')) {{
    var popup = e.target.nextElementSibling;
    var isVisible = popup.classList.contains('visible');
    document.querySelectorAll('.info-popup.visible').forEach(function(p) {{ p.classList.remove('visible'); }});
    if (!isVisible) {{
      popup.style.visibility = 'hidden';
      popup.style.display = 'block';
      var rect = e.target.getBoundingClientRect();
      var pw = popup.offsetWidth, ph = popup.offsetHeight;
      var top = rect.top - ph - 8;
      if (top < 10) top = rect.bottom + 8;
      var left = rect.left;
      if (left + pw > window.innerWidth - 10) left = window.innerWidth - pw - 10;
      if (left < 10) left = 10;
      popup.style.left = left + 'px';
      popup.style.top  = top  + 'px';
      popup.style.display = '';
      popup.style.visibility = '';
      popup.classList.add('visible');
    }}
    e.stopPropagation();
  }} else if (e.target.classList.contains('info-popup-close')) {{
    e.target.closest('.info-popup').classList.remove('visible');
  }} else if (e.target.closest && e.target.closest('.info-popup')) {{
    return;
  }} else {{
    document.querySelectorAll('.info-popup.visible').forEach(function(p) {{ p.classList.remove('visible'); }});
  }}
}});
var _drag = {{ el: null, sx: 0, sy: 0, ox: 0, oy: 0 }};
document.addEventListener('mousedown', function(e) {{
  var hdr = e.target.closest && e.target.closest('.info-popup-header');
  if (hdr && !e.target.classList.contains('info-popup-close')) {{
    _drag.el = hdr.parentElement;
    _drag.sx = e.clientX; _drag.sy = e.clientY;
    _drag.ox = parseInt(_drag.el.style.left) || 0;
    _drag.oy = parseInt(_drag.el.style.top)  || 0;
    e.preventDefault();
  }}
}});
document.addEventListener('mousemove', function(e) {{
  if (_drag.el) {{
    _drag.el.style.left = (_drag.ox + e.clientX - _drag.sx) + 'px';
    _drag.el.style.top  = (_drag.oy + e.clientY - _drag.sy) + 'px';
  }}
}});
document.addEventListener('mouseup', function() {{ _drag.el = null; }});
</script>
</body>
</html>
"""

_SECTION_TEMPLATE = """\
<section>
  <h2>{title}</h2>
  {content}
</section>
"""

_PLOT_GRID_TEMPLATE = '<div class="plot-grid">{cards}</div>'

_PLOT_CARD_TEMPLATE = """\
<div class="plot-card">
  <img src="data:image/png;base64,{b64}" alt="{caption}">
  <p>{caption}</p>
</div>
"""

_PLOT_CARD_INFO_TEMPLATE = """\
<div class="plot-card">
  <img src="data:image/png;base64,{b64}" alt="{caption}">
  <div class="plot-caption">
    <span>{caption}</span>
    <button class="info-btn" aria-label="More info">i</button>
    <div class="info-popup">
      <div class="info-popup-header">
        <span>{caption}</span>
        <button class="info-popup-close" aria-label="Close">&#x2715;</button>
      </div>
      <div class="info-popup-body">{info_html}</div>
    </div>
  </div>
</div>
"""

_VIDEO_CARD_MP4_TEMPLATE = """\
<div class="video-card">
  <video controls loop autoplay muted onloadedmetadata="this.playbackRate=1.75;">
    <source src="data:video/mp4;base64,{b64}" type="video/mp4">
    Your browser does not support embedded video.
  </video>
  <p>{caption}</p>
</div>
"""

_VIDEO_CARD_GIF_TEMPLATE = """\
<div class="video-card">
  <img src="data:image/gif;base64,{b64}" alt="{caption}">
  <p>{caption}</p>
</div>
"""


# ---------------------------------------------------------------------------
# Per-plot tooltip content
# ---------------------------------------------------------------------------

_END_REASON_INFO = (
    '<table>'
    '<thead><tr><th>End Reason</th><th>Meaning</th></tr></thead>'
    '<tbody>'
    '<tr><td>signal_positive</td><td>Normal end — strand finished translocating through the pore. The most desirable outcome.</td></tr>'
    '<tr><td>signal_negative</td><td>Read terminated by a negative ionic current event; often a stalled or damaged strand.</td></tr>'
    '<tr><td>mux_change</td><td>Read cut short by a routine MinKNOW mux scan (~every 1.5 hr) to assess pore health.</td></tr>'
    '<tr><td>unblock_mux_change</td><td>Pore actively unblocked (voltage reversed to eject a stuck strand), then a mux scan followed.</td></tr>'
    '<tr><td>data_service_unblock_mux_change</td><td>Unblock triggered by adaptive sampling / data service (e.g. ReadFish or Dorado real-time rejection).</td></tr>'
    '<tr><td>unknown</td><td>End reason undetermined; common in older Guppy or Albacore files.</td></tr>'
    '<tr><td>forced</td><td>MinKNOW forcibly ended the read, typically due to a pore quality threshold being crossed.</td></tr>'
    '<tr><td>user_regulated_pause</td><td>The run was paused manually by the user mid-read.</td></tr>'
    '</tbody></table>'
)

_DUTY_TIME_INFO = (
    '<table>'
    '<thead><tr><th>State</th><th>Colour</th><th>Meaning</th></tr></thead>'
    '<tbody>'
    '<tr><td>strand</td><td><span style="color:#2196F3">&#9632;</span> blue</td><td>DNA strand actively translocating — productive sequencing time.</td></tr>'
    '<tr><td>pore</td><td><span style="color:#4CAF50">&#9632;</span> green</td><td>Pore open and available but no strand present. Healthy idle state.</td></tr>'
    '<tr><td>adapter</td><td><span style="color:#8BC34A">&#9632;</span> light green</td><td>Adapter threading into pore — transitional state just before a read begins.</td></tr>'
    '<tr><td>unblocking</td><td><span style="color:#FF9800">&#9632;</span> orange</td><td>Reverse voltage applied to eject a stuck molecule from a blocked pore.</td></tr>'
    '<tr><td>no_pore</td><td><span style="color:#9E9E9E">&#9632;</span> grey</td><td>No functional pore in this channel — well is empty or pore was lost.</td></tr>'
    '<tr><td>multiple</td><td><span style="color:#9C27B0">&#9632;</span> purple</td><td>Multiple pores detected in one channel — excluded from sequencing.</td></tr>'
    '<tr><td>saturated</td><td><span style="color:#F44336">&#9632;</span> red</td><td>Signal saturated (clipped) — usually a bubble or debris in the well.</td></tr>'
    '<tr><td>locked</td><td><span style="color:#795548">&#9632;</span> brown</td><td>Pore blocked by a non-DNA molecule; cannot be automatically unblocked.</td></tr>'
    '<tr><td>unavailable</td><td><span style="color:#607D8B">&#9632;</span> blue-grey</td><td>Temporarily unavailable — typically during mux scans or device operations.</td></tr>'
    '<tr><td>zero</td><td><span style="color:#212121">&#9632;</span> near-black</td><td>Zero current detected; pore present but inactive.</td></tr>'
    '<tr><td>pending_manual_reset</td><td><span style="color:#FF5722">&#9632;</span> deep orange</td><td>Requires a manual voltage reset before the channel can be reused.</td></tr>'
    '<tr><td>pending_mux_change</td><td><span style="color:#FFC107">&#9632;</span> amber</td><td>Waiting for the next mux scan before MinKNOW reassigns the channel.</td></tr>'
    '<tr><td>unclassified</td><td><span style="color:#9E9E9E">&#9632;</span> light grey</td><td>State could not be classified by MinKNOW.</td></tr>'
    '<tr><td>unclassified_following_reset</td><td><span style="color:#BDBDBD">&#9632;</span> grey</td><td>Unclassified immediately after a voltage reset event.</td></tr>'
    '<tr><td>unknown_negative</td><td><span style="color:#B71C1C">&#9632;</span> dark red</td><td>Unknown state with negative current bias — usually transient.</td></tr>'
    '<tr><td>unknown_positive</td><td><span style="color:#1B5E20">&#9632;</span> dark green</td><td>Unknown state with positive current bias — usually transient.</td></tr>'
    '<tr><td>disabled</td><td><span style="color:#000000">&#9632;</span> black</td><td>Channel manually or automatically disabled.</td></tr>'
    '</tbody></table>'
)

_INFO_TOOLTIPS = {
    'Read End Reason':          _END_REASON_INFO,
    'Pore Activity / Duty Time': _DUTY_TIME_INFO,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _img_to_b64(path):
    with open(path, 'rb') as fh:
        return base64.b64encode(fh.read()).decode('ascii')


def _video_section(video_path):
    """Return an HTML section string embedding the video as base64."""
    if not video_path or not os.path.isfile(video_path):
        return ''
    b64 = _img_to_b64(video_path)
    caption = 'Per-channel strand time animated over run time (1 frame = 1 hour)'
    if video_path.endswith('.mp4'):
        card = _VIDEO_CARD_MP4_TEMPLATE.format(b64=b64, caption=caption)
    else:
        card = _VIDEO_CARD_GIF_TEMPLATE.format(b64=b64, caption=caption)
    return _SECTION_TEMPLATE.format(title='Channel Strand-Time Animation', content=card)


def _stats_table(stats):
    rows = [
        ('Run time',              f"{stats['run_time_hr']} hr"),
        ('Pores available @ T=0', f"{stats.get('pores_t0', 'N/A'):,}"),
        ('Total reads',           f"{stats['total_reads']:,}"),
        ('Pass-filter reads',     f"{stats['pf_reads']:,}"),
        ('% pass-filter',         f"{stats['pct_pf']} %"),
        ('Total yield',           f"{stats['yield_gb']} Gb"),
        ('Pass-filter yield',     f"{stats['pf_yield_gb']} Gb"),
        ('N50',                   f"{stats['n50']:,} bp" if stats['n50'] else 'N/A'),
        ('Mean PF read length',   f"{stats['mean_read_length']:,} bp"),
        ('SD PF read length',     f"{stats['sd_read_length']:,} bp"),
        ('Longest PF read',       f"{stats['max_read_length']:,} bp"),
        ('Mean PF Q-score',       str(stats['mean_qscore'])),
        (f"Prop. bases ≥ {stats['cutoff_length']:,} bp",
                                  str(stats['prop_above_cutoff'])),
    ]
    tbody = ''.join(
        f'<tr><td>{k}</td><td>{v}</td></tr>' for k, v in rows
    )
    return (
        '<table class="stats-table">'
        '<thead><tr><th>Metric</th><th>Value</th></tr></thead>'
        f'<tbody>{tbody}</tbody>'
        '</table>'
    )


def _top_reads_table(top10_df):
    """Render the top-10 longest reads as a simple HTML table."""
    if top10_df is None or top10_df.empty:
        return '<p>No pass-filter reads found.</p>'
    rows = ''
    for rank, row in top10_df.iterrows():
        length = f"{int(row['sequence_length_template']):,}"
        rows += f'<tr><td>{rank}</td><td>{row["read_id"]}</td><td>{length}</td></tr>'
    return (
        '<table class="reads-table">'
        '<thead><tr><th>#</th><th>Read ID</th><th>Length (bp)</th></tr></thead>'
        f'<tbody>{rows}</tbody>'
        '</table>'
    )


def _plot_cards(paths_captions, info_tooltips=None):
    """paths_captions: list of (path, caption) tuples. Skip None paths."""
    info_tooltips = info_tooltips or {}
    cards = []
    for path, caption in paths_captions:
        if path and os.path.isfile(path):
            b64 = _img_to_b64(path)
            if caption in info_tooltips:
                cards.append(_PLOT_CARD_INFO_TEMPLATE.format(
                    b64=b64, caption=caption, info_html=info_tooltips[caption]
                ))
            else:
                cards.append(_PLOT_CARD_TEMPLATE.format(b64=b64, caption=caption))
    return _PLOT_GRID_TEMPLATE.format(cards='\n'.join(cards))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_report(stats, plot_registry, outpath, run_name, video_path=None,
                 plotly_sections=None):
    """
    Build a self-contained HTML QC report.

    Parameters
    ----------
    stats : dict
        Output of stats.compute_stats().
    plot_registry : dict
        Keys are section names; values are lists of (filepath, caption) tuples.
        Recognised section keys (in order):
          'throughput', 'read_quality', 'channel_map', 'barcodes', 'duty_time'
    outpath : str
        Destination HTML file path.
    run_name : str
        Run identifier shown in the report header.
    """
    sections_html = []

    # 1. Summary statistics
    sections_html.append(
        _SECTION_TEMPLATE.format(
            title='Run Summary',
            content=_stats_table(stats),
        )
    )

    # 2. Ordered plot sections
    section_order = [
        ('throughput',   'Throughput'),
        ('run_health',   'Run Health over Time'),
        ('read_quality', 'Read Quality'),
        ('channel_map',  'Flowcell Channel Map'),
        ('barcodes',     'Barcode Analysis'),
        ('duty_time',    'Pore Activity / Duty Time'),
    ]

    plotly_sections = plotly_sections or {}

    for key, title in section_order:
        if key in plotly_sections and plotly_sections[key]:
            # Plotly chart first; any additional static PNGs (e.g. occupancy) follow
            parts = [f'<div class="plotly-card">{plotly_sections[key]}</div>']
            if key in plot_registry and plot_registry[key]:
                parts.append(
                    f'<div style="margin-top: 20px;">{_plot_cards(plot_registry[key], info_tooltips=_INFO_TOOLTIPS)}</div>'
                )
            sections_html.append(
                _SECTION_TEMPLATE.format(title=title, content=''.join(parts))
            )
        elif key in plot_registry and plot_registry[key]:
            content = _plot_cards(plot_registry[key], info_tooltips=_INFO_TOOLTIPS)
            if content.strip():
                sections_html.append(
                    _SECTION_TEMPLATE.format(title=title, content=content)
                )
        # Insert video section after channel map
        if key == 'channel_map':
            vid_html = _video_section(video_path)
            if vid_html:
                sections_html.append(vid_html)

    # Last section: top-10 longest reads
    sections_html.append(
        _SECTION_TEMPLATE.format(
            title='Top 10 Longest Pass-Filter Reads',
            content=_top_reads_table(stats.get('top10_reads')),
        )
    )

    html = _PAGE_TEMPLATE.format(
        run_name=run_name or 'Unknown',
        generated=datetime.now().strftime('%Y-%m-%d %H:%M'),
        sections='\n'.join(sections_html),
    )

    with open(outpath, 'w', encoding='utf-8') as fh:
        fh.write(html)

    return outpath
