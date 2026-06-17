#!/bin/bash
# ---------------------------------------------------------------------------
# check_active_runs.sh
#
# SSH to the PromethION and print a quick summary of any runs currently
# in progress.  Active runs are identified by open (*.txt.tmp) sequencing
# summary files under /data/reads/tmp.
#
# Usage:  bash check_active_runs.sh
# ---------------------------------------------------------------------------

PROM_HOST="prom@10.13.3.129"
DATA_DIR="/data/reads/tmp"

echo "PromethION run monitor -- $(date)"
echo "Host: $PROM_HOST"
echo ""

ssh -o ConnectTimeout=10 "$PROM_HOST" DATA_DIR="$DATA_DIR" bash <<'REMOTE'

NOW=$(date +%s)

# Only consider files modified in the last ACTIVE_MINS minutes.
# This excludes stale .txt.tmp files left behind by old/crashed runs.
ACTIVE_MINS=120

mapfile -t RUN_DIRS < <(
    find "$DATA_DIR" -name "sequencing_summary*.txt.tmp" -type f \
         -mmin "-${ACTIVE_MINS}" 2>/dev/null \
    | xargs -I{} dirname {} \
    | sort -u
)

if [ "${#RUN_DIRS[@]}" -eq 0 ]; then
    echo "No active runs found."
    exit 0
fi

echo "Active runs: ${#RUN_DIRS[@]}"
echo "------------------------------------------------------------"

for run_dir in "${RUN_DIRS[@]}"; do

    RUN_ID=$(basename "$run_dir")
    SAMPLE=$(basename "$(dirname "$run_dir")")

    # Run start time -- directory birth time, fall back to ctime
    START_EPOCH=$(stat -c '%W' "$run_dir" 2>/dev/null)
    [ "${START_EPOCH:-0}" -eq 0 ] && START_EPOCH=$(stat -c '%Z' "$run_dir" 2>/dev/null)
    ELAPSED=$(( NOW - START_EPOCH ))
    ELAPSED_H=$(( ELAPSED / 3600 ))
    ELAPSED_M=$(( (ELAPSED % 3600) / 60 ))

    # Active-run indicator: the open .txt.tmp file
    f=$(find "$run_dir" -name "sequencing_summary*.txt.tmp" -type f | head -1)
    [ -z "$f" ] && continue

    # Last-write time (from the open file only)
    MOD_EPOCH=$(stat -c '%Y' "$f" 2>/dev/null)
    IDLE=$(( NOW - MOD_EPOCH ))
    if   [ "$IDLE" -lt 60 ];   then LAST="${IDLE}s ago"
    elif [ "$IDLE" -lt 3600 ]; then LAST="$(( IDLE/60 ))m ago"
    else                            LAST="$(( IDLE/3600 ))h $(( (IDLE%3600)/60 ))m ago"
    fi
    [ "$IDLE" -gt 7200 ] && LAST="$LAST  [STALLED?]"

    # Collect ALL summary files for this run: MinKNOW rotates completed chunks
    # to .txt and keeps only the current chunk open as .txt.tmp, so we must
    # sum across every file to get the true total yield.
    mapfile -t run_files < <(
        find "$run_dir" -name "sequencing_summary*.txt*" -type f | sort
    )

    # Full-file awk pass: yield, read count, pass rate, mean Q.
    read -r READS YIELD PASS FAIL MEANQ < <(
        awk 'FNR==1 {
                if (!lcol) {
                    for(i=1;i<=NF;i++) {
                        if($i=="sequence_length_template" && !lcol) lcol=i
                        if($i=="sequence_length") lcol=i
                        if($i=="passes_filtering")         pcol=i
                        if($i=="mean_qscore_template")     qcol=i
                    }
                }
                next
             }
             {
                reads++
                if(lcol) bases += $lcol
                if(pcol && (toupper($pcol)=="TRUE"  || $pcol=="1")) pass++
                if(pcol && (toupper($pcol)=="FALSE" || $pcol=="0")) fail++
                if(qcol) qsum += $qcol
             }
             END {
                avg_q = (reads>0) ? qsum/reads : 0
                printf "%d %.0f %d %d %.1f\n",
                       reads+0, bases+0, pass+0, fail+0, avg_q
             }
        ' FS='\t' "${run_files[@]}"
    )

    # Active pore estimates: unique channels in the first and last 10k reads.
    # head/tail make both operations fast regardless of file size.
    PORES_START=$(awk -F'\t' 'NR>1 && NR<=10001{ch[$8]++} NR>10001{exit} END{print length(ch)}' "$f")
    PORES_NOW=$(tail -n 10000 "$f" \
        | awk -F'\t' '{ch[$8]++} END{print length(ch)}')

    # Format yield
    if   [ "${YIELD:-0}" -ge 1000000000 ] 2>/dev/null; then
        YIELD_FMT=$(awk "BEGIN {printf \"%.2f Gb\", $YIELD/1e9}")
    elif [ "${YIELD:-0}" -ge 1000000 ] 2>/dev/null; then
        YIELD_FMT=$(awk "BEGIN {printf \"%.0f Mb\", $YIELD/1e6}")
    else
        YIELD_FMT="< 1 Mb"
    fi

    # Pass rate
    TOTAL_QC=$(( PASS + FAIL ))
    if [ "$TOTAL_QC" -gt 0 ]; then
        PASS_RATE=$(awk "BEGIN {printf \"%.0f\", 100*$PASS/$TOTAL_QC}")
    else
        PASS_RATE="n/a"
    fi

    echo "Sample       : $SAMPLE"
    echo "Run ID       : $RUN_ID"
    printf "Running      : %dh %02dm\n" "$ELAPSED_H" "$ELAPSED_M"
    echo "Yield        : $YIELD_FMT"
    echo "Reads        : $READS  (pass rate: ${PASS_RATE}%)"
    echo "Mean Q       : $MEANQ"
    echo "Active pores : ~$PORES_NOW now  /  ~$PORES_START at start  (est. from 10k reads)"
    echo "Updated      : $LAST"
    echo "------------------------------------------------------------"

done

REMOTE

[ $? -ne 0 ] && echo "ERROR: SSH to $PROM_HOST failed." && exit 1
