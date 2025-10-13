#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ========== CONFIG ==========
PROJ_DIR=/mnt/vol1/Mouse_model_RNA_Seq
TRIM_DIR="$PROJ_DIR/trimmed_fastq"
INDEX_DIR="$PROJ_DIR/index"
RESULTS_DIR="$PROJ_DIR/kallisto_results"
LOG_DIR="$RESULTS_DIR/logs"
THREADS=16
IDX="$INDEX_DIR/mouse_transcriptome.k31.idx" # index filename
PAIRS_FILE="$RESULTS_DIR/filenames.txt"     # two columns: R1  R2
MISSING_REPORT="$RESULTS_DIR/missing_pairs.txt"
THREECOL="$RESULTS_DIR/abundance_3col.tsv.gz"
SAMPLES_MANIFEST="$RESULTS_DIR/samples.txt"
# ============================

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"

# ---- sanity checks ----
command -v kallisto >/dev/null 2>&1 || { echo "ERROR: kallisto not found in PATH"; exit 1; }
[[ -d "$TRIM_DIR" ]] || { echo "ERROR: TRIM_DIR not found: $TRIM_DIR"; exit 1; }

# ---- preflight: build two-column pairs list & validate ----
: > "$PAIRS_FILE"
: > "$MISSING_REPORT"

echo "[INFO] Scanning for pairs in: $TRIM_DIR"
shopt -s nullglob
R1_FILES=( "$TRIM_DIR"/*_R1_trimmed.fastq.gz )
if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No R1 files found matching *_R1_trimmed.fastq.gz in $TRIM_DIR"
  exit 1
fi

MISSING=0
for R1 in "${R1_FILES[@]}"; do
  R2="${R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
  if [[ -f "$R2" ]]; then
    # write two columns: absolute paths
    printf "%s\t%s\n" "$R1" "$R2" >> "$PAIRS_FILE"
  else
    echo "Missing pair for: $(basename "$R1")  (expected $(basename "$R2"))" >> "$MISSING_REPORT"
    MISSING=1
  fi
done

if [[ "$MISSING" -ne 0 ]]; then
  echo "[ERROR] Pairing check failed. Missing seq detected:"
  cat "$MISSING_REPORT"
  exit 2
fi

PAIR_COUNT=$(wc -l < "$PAIRS_FILE" | tr -d ' ')
echo "[INFO] Pairing check OK. Found $PAIR_COUNT pairs."
cp "$PAIRS_FILE" "$RESULTS_DIR/filenames.checked.txt"

# ---- build mouse cDNA index (Ensembl r114, GRCm39) once ----
if [[ ! -f "$IDX" ]]; then
  echo "[INFO] Index not found. Building at: $IDX"
  cd "$INDEX_DIR"
  if [[ ! -f Mus_musculus.GRCm39.cdna.all.fa.gz ]]; then
    echo "[INFO] Downloading Ensembl r115 mouse cDNA FASTA…"
    wget -O Mus_musculus.GRCm39.cdna.all.fa.gz \
      https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
  fi
  gunzip -c Mus_musculus.GRCm39.cdna.all.fa.gz > Mus_musculus.GRCm39.cdna.all.fa
  kallisto index -i "$IDX" Mus_musculus.GRCm39.cdna.all.fa
  echo "[INFO] Index built."
fi

# ---- quantify all pairs ----
: > "$SAMPLES_MANIFEST"
echo "[INFO] Starting quantification with $THREADS threads…"

# reads two names per line (two columns from the preflight list)
while IFS=$'\t' read -r FQZ1 FQZ2; do
  # infer sample name from R1 base (drop suffix)
  base=$(basename "$FQZ1")
  sample="${base%_R1_trimmed.fastq.gz}"
  OUTDIR="$RESULTS_DIR/$sample"
  mkdir -p "$OUTDIR"
  echo "$sample" >> "$SAMPLES_MANIFEST"

  echo "[INFO] $sample"
  {
    echo "kallisto quant"
    echo "sample: $sample"
    echo "R1: $FQZ1"
    echo "R2: $FQZ2"
    echo "threads: $THREADS"
    date
    kallisto quant \
      -i "$IDX" \
      -o "$OUTDIR" \
      -t "$THREADS" \
      "$FQZ1" "$FQZ2"
    date
    echo "done $sample"
  } &> "$LOG_DIR/${sample}.kallisto.log"

done < "$PAIRS_FILE"

echo "[INFO] Quantification complete."
























