#!/usr/bin/env bash
set -euo pipefail

# Config: paths + basic settings

# Project root
PROJECT_ROOT="/scratch/bx96/ma9796/rnaseq_project"

# concatenated raw FASTQs live
FASTQ_DIR="${PROJECT_ROOT}/raw"

# Trimmed reads
TRIM_DIR="${PROJECT_ROOT}/trimmed"

# QC outputs
QC_FASTP_DIR="${PROJECT_ROOT}/qc/fastp"
QC_FASTQC_DIR="${PROJECT_ROOT}/qc/fastqc_trimmed"
QC_MULTIQC_DIR="${PROJECT_ROOT}/qc/multiqc"

# Threads to use per sample
THREADS=6   # match this to ncpus in PBS job

# Trimming thresholds (fastp)
Q_CUTOFF=20        # minimum base quality
MIN_LENGTH=35      # discard reads shorter than this after trimming

# Setup

mkdir -p "$TRIM_DIR" "$QC_FASTP_DIR" "$QC_FASTQC_DIR" "$QC_MULTIQC_DIR"

# Loop over samples (paired-end)


for R1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
    base="$(basename "$R1" _R1.fastq.gz)"
    R2="${FASTQ_DIR}/${base}_R2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "Skipping ${base}: missing adj file $R2"
        continue
    fi

    echo "Processing sample: $base"
    echo "  R1: $R1"
    echo "  R2: $R2"

    # Output filenames
    R1_OUT="${TRIM_DIR}/${base}_R1.trimmed.fastq.gz"
    R2_OUT="${TRIM_DIR}/${base}_R2.trimmed.fastq.gz"
    HTML_REPORT="${QC_FASTP_DIR}/${base}_fastp.html"
    JSON_REPORT="${QC_FASTP_DIR}/${base}_fastp.json"

    #Trim with fastp (adapters + quality)

    echo "  [fastp] Trimming and QC..."
    fastp \
        -i "$R1" -I "$R2" \
        -o "$R1_OUT" -O "$R2_OUT" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --qualified_quality_phred "$Q_CUTOFF" \
        --length_required "$MIN_LENGTH" \
        --html "$HTML_REPORT" \
        --json "$JSON_REPORT"

    echo "  [fastp] Done for $base"
    echo "    Trimmed: $R1_OUT"
    echo "             $R2_OUT"
    echo "    Reports: $HTML_REPORT"
    echo "             $JSON_REPORT"

    # FastQC on trimmed reads
    echo "  [FastQC] Running on trimmed reads..."
    fastqc \
        --threads "$THREADS" \
        --outdir "$QC_FASTQC_DIR" \
        "$R1_OUT" "$R2_OUT"

    echo "  [FastQC] Done for $base"
    echo
done

# MultiQC summary over all QC reports

echo "[MultiQC] Summarising fastp + FastQC reports..."

multiqc \
    "$QC_FASTP_DIR" "$QC_FASTQC_DIR" \
    -o "$QC_MULTIQC_DIR"

echo "[MultiQC] Report: ${QC_MULTIQC_DIR}/multiqc_report.html"
echo "All trimming and QC completed."

