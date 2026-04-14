#!/bin/bash
#$ -N ont_qc
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=04:00:00
#$ -o ont_qc.$JOB_ID.out
#$ -e ont_qc.$JOB_ID.err

# ---------------------------------------------------------------------------
# SGE submission script for ont_qc.py
#
# Usage:
#   qsub submit_ont_qc.sh
#
# Customise the variables in the USER SETTINGS section below before
# submitting.  All ont_qc.py arguments are optional except --runName.
# ---------------------------------------------------------------------------

# ===========================================================================
# USER SETTINGS — edit these before submitting
# ===========================================================================

# Path to the ont_qc.py script
SCRIPT_DIR=~/Desktop/Oxford-Nanopore-QC-scripts

# Run name (used as filename prefix and in the HTML report)
RUN_NAME="MyRun"

# Path to sequencing_summary*.txt  (leave blank for auto-detection)
SUMMARY_FILE=""

# Path to pore_activity*.csv  (leave blank for auto-detection)
PORE_ACTIVITY=""

# Path to throughput_*.csv  (leave blank for auto-detection)
THROUGHPUT=""

# Output directory  (leave blank for default: <RUN_NAME>_qc/)
OUTDIR=""

# Subsample fraction (1.0 = use all reads)
SUBSAMPLE=1.0

# Read-length axis cap in bp  (leave blank for auto)
MAX_LENGTH=""

# Minimum read length shown on length-axis plots  (leave blank for auto)
MIN_LENGTH=""

# Proportion-above-cutoff threshold in bp  (default: 2000)
PROP=2000

# Barcode restriction: integer for top-N, or space-separated barcode names
# e.g. BARCODES="12"  or  BARCODES="barcode01 barcode02"
# Leave blank to include all barcodes
BARCODES=""

# Set to 1 to generate an animated channel video (requires ffmpeg on PATH)
VIDEO=0

# ===========================================================================
# BUILD COMMAND — no edits needed below this line
# ===========================================================================

CMD="python ${SCRIPT_DIR}/ont_qc.py --runName ${RUN_NAME}"

[ -n "$SUMMARY_FILE"  ] && CMD="$CMD --file $SUMMARY_FILE"
[ -n "$PORE_ACTIVITY" ] && CMD="$CMD --poreActivity $PORE_ACTIVITY"
[ -n "$THROUGHPUT"    ] && CMD="$CMD --throughput $THROUGHPUT"
[ -n "$OUTDIR"        ] && CMD="$CMD --outdir $OUTDIR"
[ -n "$MAX_LENGTH"    ] && CMD="$CMD --maxLength $MAX_LENGTH"
[ -n "$MIN_LENGTH"    ] && CMD="$CMD --minLength $MIN_LENGTH"
[ "$SUBSAMPLE" != "1.0" ] && CMD="$CMD --subsample $SUBSAMPLE"
[ "$PROP" != "2000"   ] && CMD="$CMD --prop $PROP"
[ -n "$BARCODES"      ] && CMD="$CMD --barcodes $BARCODES"
[ "$VIDEO" -eq 1      ] && CMD="$CMD --video"

echo "=============================="
echo "Job:    $JOB_NAME ($JOB_ID)"
echo "Host:   $(hostname)"
echo "Date:   $(date)"
echo "Command: $CMD"
echo "=============================="

eval "$CMD"
exit $?
