#!/usr/bin/env python3
"""
plot_active_runs.py

Single SSH call to the PromethION: discovers active runs, aggregates
read counts and yield into hourly buckets, samples read lengths (piggy-
backing on the same awk pass at no extra I/O cost), and collects recent
barcode counts, then generates diagnostic plots per run:
  1. Cumulative yield (Gb) vs. time
  2. Reads per hour
  3. Read length distribution + N50 estimate  (every 500th read sampled)
  4. Barcode distribution from last 10k reads (barcoded runs only)

Usage:
    python plot_active_runs.py [--outdir <dir>] [--min-barcode-pct <float>]

Output PNGs are saved to ./active_run_plots/ by default.
"""

import argparse
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
PROM_HOST        = "prom@10.13.3.129"
DATA_DIR         = "/data/reads/tmp"
ACTIVE_MINS      = 120    # ignore files not modified within this many minutes
SAMPLE_STEP      = 500    # sample every Nth read for length distribution
MIN_BARCODE_PCT  = 1.0    # hide barcodes below this % of last-10k reads


# ---------------------------------------------------------------------------
# Single SSH call: hourly stats + length sample + barcode counts
# ---------------------------------------------------------------------------
def fetch_all_run_data():
    """
    One SSH connection.  For each active run the remote script:
      - runs a single awk pass over the full summary file, collecting
        hourly yield/read counts AND printing every 500th read length
        inline (same pass, no extra I/O)
      - runs awk over the last 10k reads for barcode distribution

    Returns list of (sample, run_id, hourly_rows, lengths, barcode_counts)
      hourly_rows    = [(hour, read_count, bases_float), ...]
      lengths        = [int, ...]  sampled read lengths
      barcode_counts = {barcode: count} or None
    """
    remote_script = (
        f"ACTIVE_MINS={ACTIVE_MINS}\n"
        f"DATA_DIR={DATA_DIR}\n"
        f"SAMPLE_STEP={SAMPLE_STEP}\n"
        "\n"
        "find \"$DATA_DIR\" -name 'sequencing_summary*.txt.tmp' "
            "-type f -mmin \"-$ACTIVE_MINS\" 2>/dev/null | sort | "
            "while IFS= read -r f; do\n"
        "    run_dir=$(dirname \"$f\")\n"
        "    run_id=$(basename \"$run_dir\")\n"
        "    sample=$(basename \"$(dirname \"$run_dir\")\")\n"
        "\n"
        "    # -- single awk pass: hourly stats + sampled lengths ----------\n"
        "    printf 'RUN\\t%s\\t%s\\n' \"$sample\" \"$run_id\"\n"
        "    awk -F'\\t' -v step=\"$SAMPLE_STEP\" '\n"
        "    NR==1 {\n"
        "        for (i=1; i<=NF; i++) {\n"
        "            if ($i == \"start_time\")               tcol=i\n"
        "            if ($i == \"sequence_length_template\") lcol=i\n"
        "        }\n"
        "        next\n"
        "    }\n"
        "    tcol {\n"
        "        n++\n"
        "        h = int($tcol / 3600)\n"
        "        reads[h]++\n"
        "        if (lcol) bases[h] += $lcol\n"
        "        if (lcol && n % step == 0 && $lcol+0 > 0)\n"
        "            printf \"L %d\\n\", $lcol\n"
        "    }\n"
        "    END {\n"
        "        for (h in reads)\n"
        "            printf \"%d %d %.0f\\n\", h, reads[h]+0, bases[h]+0\n"
        "    }\n"
        "    ' \"$f\" | sort -k1,1\n"
        "    printf 'END\\n'\n"
        "\n"
        "    # -- barcode counts from last 10k reads -----------------------\n"
        "    BCOL=$(head -1 \"$f\" | awk -F'\\t' '{\n"
        "        for(i=1;i<=NF;i++) if($i==\"barcode_arrangement\"){print i; exit}\n"
        "    }')\n"
        "    if [ -n \"$BCOL\" ]; then\n"
        "        printf 'BARCODES\\t%s\\t%s\\n' \"$sample\" \"$run_id\"\n"
        "        tail -n 10000 \"$f\" | awk -F'\\t' -v bcol=\"$BCOL\" \\\n"
        "            '{counts[$bcol]++}\n"
        "             END{for(bc in counts) printf \"%s\\t%d\\n\", bc, counts[bc]}' \\\n"
        "            | sort\n"
        "        printf 'BARCODES_END\\n'\n"
        "    fi\n"
        "done\n"
    )

    proc = subprocess.Popen(
        ["ssh", "-o", "ConnectTimeout=30", PROM_HOST, "bash"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    proc.stdin.write(remote_script)
    proc.stdin.close()

    run_data  = {}
    run_order = []

    current_run  = None
    current_rows = []
    current_lens = []
    in_barcodes  = False
    bc_target    = None
    bc_counts    = {}

    for raw in proc.stdout:
        line = raw.rstrip()

        if line.startswith("RUN\t"):
            parts = line.split("\t")
            if len(parts) == 3:
                current_run  = (parts[1], parts[2])
                current_rows = []
                current_lens = []
                print(f"  Processing {parts[1]} ...", flush=True)

        elif line == "END":
            if current_run is not None:
                run_data[current_run] = {
                    "rows": current_rows,
                    "lengths": current_lens,
                    "barcodes": None,
                }
                run_order.append(current_run)
                print(
                    f"    {len(current_rows)} hours, "
                    f"{len(current_lens)} sampled lengths.",
                    flush=True,
                )
            current_run  = None
            current_rows = []
            current_lens = []

        elif line.startswith("BARCODES\t"):
            parts = line.split("\t")
            if len(parts) == 3:
                bc_target   = (parts[1], parts[2])
                bc_counts   = {}
                in_barcodes = True

        elif line == "BARCODES_END":
            if bc_target and bc_target in run_data:
                run_data[bc_target]["barcodes"] = bc_counts
            in_barcodes = False
            bc_target   = None
            bc_counts   = {}

        elif in_barcodes:
            parts = line.split("\t")
            if len(parts) == 2:
                try:
                    bc_counts[parts[0]] = int(parts[1])
                except ValueError:
                    pass

        elif line.startswith("L "):
            try:
                current_lens.append(int(line[2:]))
            except ValueError:
                pass

        elif current_run is not None:
            parts = line.split()
            if len(parts) == 3:
                try:
                    current_rows.append(
                        (int(parts[0]), int(parts[1]), float(parts[2]))
                    )
                except ValueError:
                    pass

    proc.wait()

    return [
        (
            s, r,
            run_data[(s, r)]["rows"],
            run_data[(s, r)]["lengths"],
            run_data[(s, r)]["barcodes"],
        )
        for s, r in run_order
    ]


# ---------------------------------------------------------------------------
# N50 estimate from sampled lengths
# ---------------------------------------------------------------------------
def estimate_n50(lengths):
    """N50 of the sampled lengths -- a good proxy for the full dataset."""
    if not lengths:
        return 0
    sorted_l = sorted(lengths, reverse=True)
    target   = sum(sorted_l) / 2
    cumsum   = 0
    for l in sorted_l:
        cumsum += l
        if cumsum >= target:
            return l
    return 0


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_run(sample, run_id, rows, lengths, barcodes, outdir, min_barcode_pct):
    if not rows:
        print(f"  No data for {sample} -- skipping.")
        return

    # --- hourly data --------------------------------------------------------
    max_h      = max(h for h, _, _ in rows)
    reads_by_h = [0] * (max_h + 1)
    bases_by_h = [0.0] * (max_h + 1)
    for h, r, b in rows:
        reads_by_h[h] = r
        bases_by_h[h] = b

    all_hours = list(range(max_h + 1))
    cum_gb, total = [], 0.0
    for b in bases_by_h:
        total += b
        cum_gb.append(total / 1e9)

    # --- barcode data -------------------------------------------------------
    plot_barcodes = None
    if barcodes:
        bc_total  = sum(barcodes.values())
        threshold = bc_total * (min_barcode_pct / 100.0)
        filtered  = {
            bc: cnt for bc, cnt in barcodes.items()
            if cnt >= threshold and bc not in ("not_set", "")
        }
        if len({bc for bc in filtered if bc != "unclassified"}) >= 2:
            plot_barcodes = filtered

    # --- figure: always 2x2, hide unused panels ----------------------------
    for style_name in ("seaborn-v0_8", "seaborn"):
        try:
            plt.style.use(style_name)
            break
        except OSError:
            continue

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1, ax2  = axes[0]
    ax3, ax4  = axes[1]
    fig.suptitle(sample, fontsize=13, fontweight="bold")

    # Plot 1: cumulative yield
    ax1.plot(all_hours, cum_gb, linewidth=2, color="#4472C4")
    ax1.set_xlabel("Time (hr)")
    ax1.set_ylabel("Yield (Gb)")
    ax1.set_title("Cumulative Yield vs. Time")

    # Plot 2: reads per hour
    ax2.bar(all_hours, reads_by_h, width=0.8, color="#4472C4", alpha=0.8)
    ax2.set_xlabel("Run Time (hr)")
    ax2.set_ylabel("Read Count")
    ax2.set_title("Reads per Hour")

    # Plot 3: read length distribution
    if lengths:
        n50      = estimate_n50(lengths)
        sorted_l = sorted(lengths)
        cap      = sorted_l[int(len(sorted_l) * 0.95)]   # trim top 5% outliers
        plot_l   = [l for l in lengths if l <= cap]

        ax3.hist(plot_l, bins=100, color="#4472C4", alpha=0.8, edgecolor="none")
        ax3.axvline(n50, color="#C0392B", linewidth=1.5, linestyle="--",
                    label=f"N50 ~{n50:,} bp")
        ax3.set_xlabel("Read Length (bp)")
        ax3.set_ylabel("Count (sampled)")
        ax3.set_title(
            f"Read Length Distribution  (1 in {SAMPLE_STEP} reads sampled)"
        )
        ax3.legend(fontsize=9)
    else:
        ax3.set_visible(False)

    # Plot 4: barcode distribution (last 10k reads)
    if plot_barcodes:
        def bc_sort_key(bc):
            return (bc == "unclassified", bc)

        labels = sorted(plot_barcodes.keys(), key=bc_sort_key)
        counts = [plot_barcodes[bc] for bc in labels]
        colors = ["#A8C8E8" if bc == "unclassified" else "#4472C4"
                  for bc in labels]

        ax4.bar(labels, counts, color=colors)
        ax4.set_xlabel("Barcode")
        ax4.set_ylabel("Read Count")
        ax4.set_title(
            f"Barcode Distribution  (last 10k reads, "
            f">{min_barcode_pct:.0f}% threshold)"
        )
        ax4.tick_params(axis="x", rotation=45)
    else:
        ax4.set_visible(False)

    plt.tight_layout()
    fname = os.path.join(outdir, f"{sample}.png")
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir", default="active_run_plots",
        help="Directory for output PNGs (default: active_run_plots/)"
    )
    parser.add_argument(
        "--min-barcode-pct", type=float, default=MIN_BARCODE_PCT,
        help=(
            f"Hide barcodes below this %% of last-10k reads "
            f"(default: {MIN_BARCODE_PCT})"
        ),
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Connecting to {PROM_HOST} (one password prompt)...")
    all_runs = fetch_all_run_data()

    if not all_runs:
        print("No active runs found.")
        return

    print(f"\nFound {len(all_runs)} active run(s). Generating plots...\n")
    for sample, run_id, rows, lengths, barcodes in all_runs:
        print(f"  {sample}  ({run_id})")
        plot_run(
            sample, run_id, rows, lengths, barcodes,
            args.outdir, args.min_barcode_pct,
        )

    print(f"\nDone. Plots saved to: {os.path.abspath(args.outdir)}/")


if __name__ == "__main__":
    main()
