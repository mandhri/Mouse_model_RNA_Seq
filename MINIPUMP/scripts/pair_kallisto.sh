#!/usr/bin/env bash
set -euo pipefail

# Project layout
PROJ_DIR="/mnt/vol1/Mouse_model_RNA_Seq"
TRIM_DIR="${PROJ_DIR}/trimmed_fastq"
INDEX_DIR="${PROJ_DIR}/index"
RESULTS_DIR="${PROJ_DIR}/kallisto_results"
LOG_DIR="${RESULTS_DIR}/logs"

# Kallisto settings
THREADS=16
IDX="${INDEX_DIR}/mouse_transcriptome.k31.idx"   # cDNA-only index (Ensembl r115)

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"


# sanity checks 
command -v kallisto >/dev/null 2>&1 || { echo "ERROR: kallisto not found in PATH"; exit 1; }
[[ -d "$TRIM_DIR" ]] || { echo "ERROR: TRIM_DIR not found: $TRIM_DIR"; exit 1; }

# BUILD INDEX (cDNA only)

if [[ ! -f "$IDX" ]]; then
  echo "[INFO] Index not found. Building at: $IDX"
  cd "$INDEX_DIR"

  if [[ ! -f Mus_musculus.GRCm39.cdna.all.fa.gz ]]; then
    echo "[INFO] Downloading Ensembl r115 mouse cDNA FASTA…"
    wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
      "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
  fi

  # Unzip and build kallisto index
  gunzip -c Mus_musculus.GRCm39.cdna.all.fa.gz > Mus_musculus.GRCm39.cdna.all.fa
  kallisto index -i "$IDX" Mus_musculus.GRCm39.cdna.all.fa
  echo "[INFO] Index built."
fi


# QUANTIFICATION


echo "[INFO] Quantifying all samples in: $TRIM_DIR"
shopt -s nullglob

for R1 in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
  base="$(basename "$R1")"
  sample="${base%_R1_trimmed.fastq.gz}"
  R2="${TRIM_DIR}/${sample}_R2_trimmed.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "[WARN] Skipping ${sample}: missing mate file $R2"
    continue
  fi

  OUT_DIR="${RESULTS_DIR}/${sample}"
  mkdir -p "$OUT_DIR"

  # Skip if already done
  if [[ -s "${OUT_DIR}/abundance.tsv" ]]; then
    echo "[INFO] ${sample}: already quantified, skipping."
    continue
  fi

  echo "[INFO] Running kallisto for sample: ${sample}"
  echo "  R1: $R1"
  echo "  R2: $R2"
  echo "  Out: $OUT_DIR"

  kallisto quant \
    -i "$IDX" \
    -o "$OUT_DIR" \
    -t "$THREADS" \
    "$R1" "$R2"

  echo "[INFO] Finished ${sample}"
  echo
done

echo "[INFO] All samples processed."




























