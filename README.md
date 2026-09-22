# Oxford Nanopore QC Scripts

Quality-control tools for Oxford Nanopore Technology (ONT) sequencing runs. Generates plots and a self-contained HTML report from MinKNOW output files.

---

## Contents

| File | Description |
|---|---|
| `ont_qc.py` | Main QC script — produces plots and an HTML report |
| `qc_modules/` | Supporting modules imported by `ont_qc.py` |
| `plot_active_runs.py` | Quick-update plots for runs in progress (SSHes to PromethION) |
| `plot_guide.py` | Generates `ONT_QC_Plot_Guide.pdf` |
| `ONT_QC_Plot_Guide.pdf` | Reference guide describing every plot |

---

## Requirements

```
pip install matplotlib numpy pandas seaborn plotly
```

| Package | Minimum version | Used by |
|---|---|---|
| matplotlib | 3.8 | `ont_qc.py`, `plot_active_runs.py` |
| numpy | 1.26 | `ont_qc.py` |
| pandas | 2.2 | `ont_qc.py` |
| seaborn | 0.13 | `ont_qc.py` |
| plotly | 5.22 | `ont_qc.py` (interactive duty-time plot) |
| reportlab | 4.0 | `plot_guide.py` only (PDF generation) |

A `requirements.txt` is included for convenience:

```
pip install -r requirements.txt
```

> **Note:** `reportlab` is only needed if you want to regenerate `ONT_QC_Plot_Guide.pdf`. It is not required to run `ont_qc.py`.

---

## Usage

### ont_qc.py

```bash
# Minimal — auto-detect all input files in the current directory
python ont_qc.py --runName MyRun

# Explicit summary file
python ont_qc.py --file /path/to/sequencing_summary.txt --runName MyRun

# With optional companion files
python ont_qc.py --file sequencing_summary.txt \
                 --poreActivity pore_activity.csv \
                 --throughput throughput_data.csv \
                 --runName MyRun

# Custom output directory and read-length axis cap
python ont_qc.py --runName MyRun --outdir /results/qc --maxLength 50000

# Subsample to 50% of reads (useful for very large files)
python ont_qc.py --file sequencing_summary.txt --runName MyRun --subsample 0.5

# Generate an animated channel strand-time video (requires ffmpeg)
python ont_qc.py --runName MyRun --video
```

**Input files** (all auto-detected if in the same directory as the summary file):

| File | Required | Description |
|---|---|---|
| `sequencing_summary*.txt` | Yes | Per-read statistics (Albacore, Guppy, or Dorado) |
| `pore_activity*.csv` | No | Per-minute channel state counts |
| `throughput_*.csv` | No | Per-minute cumulative counters |

**Outputs** are written to `<runName>_qc/` by default:

- `<runName>_report.html` — self-contained HTML report (no external dependencies)
- `<runName>_*.png` — individual plot images
- `<runName>_summary.txt` — plain-text run statistics

### plot_active_runs.py

Generates quick-update plots (yield, read rate, length distribution, barcodes) for runs currently in progress on a PromethION. Requires SSH access to the instrument.

```bash
python plot_active_runs.py
python plot_active_runs.py --outdir my_plots/
```

---

## Plot Reference Guide

`ONT_QC_Plot_Guide.pdf` describes every plot produced by `ont_qc.py`: what it shows, what to look for, and how to interpret normal vs. abnormal findings.

To regenerate it after editing `plot_guide.py`:

```bash
pip install reportlab
python plot_guide.py
```
