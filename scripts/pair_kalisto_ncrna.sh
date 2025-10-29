#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# -------------- CONFIG -------------
PROJ_DIR=/mnt/vol1/Mouse_model_RNA_Seq
TRIM_DIR="$PROJ_DIR/trimmed_fastq"
INDEX_DIR="$PROJ_DIR/index"
RESULTS_DIR="$PROJ_DIR/kallisto_results_plus_ncRNA"
LOG_DIR="$RESULTS_DIR/logs"
THREADS=8
THREADS_FALLBACK=2

# Combined (cDNA + ncRNA) transcriptome index (Ensembl r115, GRCm39)
IDX="$INDEX_DIR/mouse_transcriptome_plus_ncRNA.k31.idx"

PAIRS_FILE="$RESULTS_DIR/filenames.txt"
MISSING_REPORT="$RESULTS_DIR/missing_pairs.txt"
SAMPLES_MANIFEST="$RESULTS_DIR/samples.txt"
FAILED_LIST="$RESULTS_DIR/failed_samples.txt"
ALIGN_SUM="$RESULTS_DIR/alignment_summary.tsv"
#------------------------------------------------------

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"

# Library
LIB_STRAND="rf"       
BIAS=""         
BOOT="-b 30 --seed 42"
SKIP_IF_DONE=1          

# basic checks
command -v kallisto >/dev/null || { echo "ERROR: kallisto not found"; exit 1; }
command -v jq >/dev/null || { echo "ERROR: jq not found (needed for summary)"; exit 1; }
[[ -d "$TRIM_DIR" ]] || { echo "ERROR: TRIM_DIR not found: $TRIM_DIR"; exit 1; }



# build pairs list (+ report missing mates)
printf "" > "$PAIRS_FILE"
printf "" > "$MISSING_REPORT"
printf "" > "$SAMPLES_MANIFEST"
printf "" > "$FAILED_LIST"


echo "[INFO] Scanning for pairs in: $TRIM_DIR"
shopt -s nullglob
for R1 in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
  R2="${R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
  if [[ -f "$R2" ]]; then
    printf "%s\t%s\n" "$R1" "$R2" >> "$PAIRS_FILE"
  else
    echo "$(basename "$R1")  -> missing $(basename "$R2")" >> "$MISSING_REPORT"
  fi
done
if [[ -s "$MISSING_REPORT" ]]; then
  echo "[ERROR] Missing mates found. See $MISSING_REPORT"
  exit 2
fi
echo "[INFO] Pairs: $(wc -l < "$PAIRS_FILE")"


# ---- build mouse cDNA+ncRNA index (Ensembl r115, GRCm39) ----
if [[ ! -s "$IDX" ]]; then
  echo "[INFO] Index not found. Building at: $IDX"
  cd "$INDEX_DIR"
  [[ -f Mus_musculus.GRCm39.cdna.all.fa.gz ]] || wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
    https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
  [[ -f Mus_musculus.GRCm39.ncrna.fa.gz ]] || wget -O Mus_musculus.GRCm39.ncrna.fa.gz \
    https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz

  zcat Mus_musculus.GRCm39.cdna.all.fa.gz Mus_musculus.GRCm39.ncrna.fa.gz > all_tx.fa
  kallisto index -i "$IDX" all_tx.fa
  rm -f all_tx.fa
  echo "[INFO] Index built."
fi

echo "[INFO] Quantifying with $THREADS threads…"

# quant loop (fallback = drop bootstraps)
while IFS=$'\t' read -r FQ1 FQ2; do
  base=$(basename "$FQ1"); sample="${base%_R1_trimmed.fastq.gz}"
  out="$RESULTS_DIR/$sample"
  log="$LOG_DIR/${sample}.kallisto.log"
  mkdir -p "$out"

  if [[ "$SKIP_IF_DONE" -eq 1 && -s "$out/abundance.tsv" && -s "$out/run_info.json" ]]; then
    echo "[INFO] $sample: already quantified, skipping."
    echo "$sample" >> "$SAMPLES_MANIFEST"
    continue
  fi

  echo "[INFO] $sample"
  strand=""; [[ -n "${LIB_STRAND:-}" ]] && strand="--${LIB_STRAND}-stranded"

  set +e
  # try 1: with bootstraps (if BOOT set)
  kallisto quant -i "$IDX" -o "$out" -t "$THREADS" ${BIAS:-} ${BOOT:-} $strand "$FQ1" "$FQ2" &> "$log"
  rc=$?
  # fallback: drop bootstraps only
  if (( rc != 0 )) && [[ -n "${BOOT:-}" ]]; then
    echo "[WARN] $sample failed (rc=$rc). Retrying without bootstraps…" | tee -a "$log"
    rm -f "$out/abundance.tsv" "$out/run_info.json"
    kallisto quant -i "$IDX" -o "$out" -t "$THREADS" ${BIAS:-} $strand "$FQ1" "$FQ2" &>> "$log"
    rc=$?
  fi
  set -e

  if (( rc == 0 )) && [[ -s "$out/abundance.tsv" && -s "$out/run_info.json" ]]; then
    echo "[INFO] done $sample"
    echo "$sample" >> "$SAMPLES_MANIFEST"
  else
    echo "[ERROR] $sample failed. See $log"
    echo "$sample" >> "$FAILED_LIST"
  fi
done < "$PAIRS_FILE"

echo "[INFO] Done."
[[ -s "$FAILED_LIST" ]] && echo "[INFO] Failures listed in $FAILED_LIST" || echo "[INFO] All samples completed."