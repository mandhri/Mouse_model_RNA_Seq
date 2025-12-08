#!/usr/bin/env bash
set -euo pipefail


# -------------- CONFIG -------------
PROJ_DIR=/mnt/vol1/Mouse_model_RNA_Seq
TRIM_DIR="$PROJ_DIR/trimmed_fastq"
INDEX_DIR="$PROJ_DIR/index"
RESULTS_DIR="$PROJ_DIR/kallisto_results_plus_ncRNA"
THREADS=8

# Combined (cDNA + ncRNA) transcriptome index (Ensembl r115, GRCm39)
IDX="$INDEX_DIR/mouse_transcriptome_plus_ncRNA.k31.idx"

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" 

if [[ ! -s "$IDX" ]]; then
  echo "[INFO] Index not found. Building at: $IDX"
  cd "$INDEX_DIR"

  # Download if fasta files are missing
  [[ -f Mus_musculus.GRCm39.cdna.all.fa.gz ]] || wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
    https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

  [[ -f Mus_musculus.GRCm39.ncrna.fa.gz ]] || wget -O Mus_musculus.GRCm39.ncrna.fa.gz \
    https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz

  # Combine cDNA + ncRNA into one transcriptome and build index
  zcat Mus_musculus.GRCm39.cdna.all.fa.gz Mus_musculus.GRCm39.ncrna.fa.gz > all_tx.fa
  kallisto index -i "$IDX" all_tx.fa
  rm -f all_tx.fa

  echo "[INFO] Index built."
fi


# Quantify each paired sample

echo "[INFO] Quantifying all samples in: $TRIM_DIR"
shopt -s nullglob

for R1 in "${TRIM_DIR}"/*_R1_trimmed.fastq.gz; do
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

  echo "[INFO] Finished ${sample]"
  echo
done

echo "[INFO] All samples processed."





