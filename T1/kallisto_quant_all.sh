#!/usr/bin/env bash
set -euo pipefail

# --- basic paths ---
PROJECT_ROOT="/scratch/bx96/ma9796/rnaseq_project"

TRIM_DIR="${PROJECT_ROOT}/trimmed"
OUT_DIR="${PROJECT_ROOT}/kallisto_results"
INDEX="${PROJECT_ROOT}/index/mouse_transcriptome_plus_ncRNA.k31.idx"

THREADS=6

mkdir -p "$OUT_DIR"

echo "Using trimmed reads from : $TRIM_DIR"
echo "Writing kallisto results : $OUT_DIR"
echo "Using index              : $INDEX"
echo "Threads                  : $THREADS"
echo

# loop over R1 files
for R1 in "${TRIM_DIR}"/*_R1.trimmed.fastq.gz; do
    base="$(basename "$R1" _R1.trimmed.fastq.gz)"
    R2="${TRIM_DIR}/${base}_R2.trimmed.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "Skipping ${base}: missing $R2" >&2
        continue
    fi

    sample_out="${OUT_DIR}/${base}"
    mkdir -p "$sample_out"

    echo "Running kallisto for sample: ${base}"
    echo "  R1: $R1"
    echo "  R2: $R2"
    echo "  out: $sample_out"
    echo

    kallisto quant \
        -i "$INDEX" \
        -o "$sample_out" \
        -t "$THREADS" \
         "$R1" "$R2"

    echo "Done: ${base}"
done

echo "All samples finished."
