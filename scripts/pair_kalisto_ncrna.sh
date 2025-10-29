#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ========== CONFIG ==========
PROJ_DIR=/mnt/vol1/Mouse_model_RNA_Seq
TRIM_DIR="$PROJ_DIR/trimmed_fastq"
INDEX_DIR="$PROJ_DIR/index"
RESULTS_DIR="$PROJ_DIR/kallisto_results_plus_ncRNA"
LOG_DIR="$RESULTS_DIR/logs"
THREADS=8
THREADS_FALLBACK=2

# Combined (cDNA + ncRNA) transcriptome index (Ensembl r115, GRCm39)
IDX="$INDEX_DIR/mouse_transcriptome_plus_ncRNA.k31.idx"

# Library / options (set empty "" to disable)
LIB_STRAND="rf"          # "rf", "fr", or "" if unknown
BIAS="--bias"            # "" to disable
BOOT="-b 50 --seed 42"   # "" to disable (drop bootstraps if you get crashes)

# Control behaviour
SKIP_IF_DONE=1           # 1 = skip sample if abundance.tsv already exists
RETRY_ON_FAIL=1          # 1 = try progressively safer re-runs on failure

PAIRS_FILE="$RESULTS_DIR/filenames.txt"
MISSING_REPORT="$RESULTS_DIR/missing_pairs.txt"
SAMPLES_MANIFEST="$RESULTS_DIR/samples.txt"
FAILED_LIST="$RESULTS_DIR/failed_samples.txt"
ALIGN_SUM="$RESULTS_DIR/alignment_summary.tsv"
# ============================

mkdir -p "$INDEX_DIR" "$RESULTS_DIR" "$LOG_DIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }; }
need kallisto
need jq
need zcat

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

# ---- build mouse cDNA+ncRNA index (Ensembl r115, GRCm39) once ----
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

# ---- init manifests ----
: > "$SAMPLES_MANIFEST"
: > "$FAILED_LIST"
echo -e "sample\tn_processed\tn_pseudoaligned\trate_percent\tn_targets\tindex_version\tkallisto_version" > "$ALIGN_SUM"

echo "[INFO] Starting quantification with $THREADS threads…"

# Helper: run one kallisto command with given args safely
run_kallisto() {
  local threads="$1"; shift
  local outdir="$1"; shift
  local log="$1"; shift
  set +e
  kallisto quant -t "$threads" "$@" &> "$log"
  local rc=$?
  set -e
  echo "$rc"
}

# Helper: parse run_info.json and update alignment summary
update_summary() {
  local sample="$1"; local outdir="$2"
  local info="$outdir/run_info.json"
  [[ -s "$info" ]] || return 1
  local npr npa pp ver idxv ntg rate
  npr=$(jq -r '.n_processed // .num_processed // empty' "$info")
  npa=$(jq -r '.n_pseudoaligned // .num_pseudoaligned // empty' "$info")
  pp=$(jq -r '.p_pseudoaligned // empty' "$info")
  ver=$(jq -r '.kallisto_version // empty' "$info")
  idxv=$(jq -r '.index_version // empty' "$info")
  ntg=$(jq -r '.n_targets // empty' "$info")

  if [[ -n "${pp:-}" ]]; then
    rate=$(awk -v p="$pp" 'BEGIN{printf "%.1f", p*100}')
  elif [[ -n "${npr:-}" && -n "${npa:-}" && "$npr" -gt 0 ]]; then
    rate=$(awk -v a="$npa" -v b="$npr" 'BEGIN{printf "%.1f", (a/b)*100}')
  else
    rate=""
  fi
  echo -e "${sample}\t${npr:-}\t${npa:-}\t${rate:-}\t${ntg:-}\t${idxv:-}\t${ver:-}" >> "$ALIGN_SUM"
}

# reads two names per line (two columns from the preflight list)
while IFS=$'\t' read -r FQZ1 FQZ2; do
  base=$(basename "$FQZ1")
  sample="${base%_R1_trimmed.fastq.gz}"
  OUTDIR="$RESULTS_DIR/$sample"
  log="$LOG_DIR/${sample}.kallisto.log"

  # skip if done
  if (( SKIP_IF_DONE )) && [[ -s "$OUTDIR/abundance.tsv" && -s "$OUTDIR/run_info.json" ]]; then
    echo "[INFO] $sample: already quantified, skipping."
    update_summary "$sample" "$OUTDIR" || true
    echo "$sample" >> "$SAMPLES_MANIFEST"
    continue
  fi

  mkdir -p "$OUTDIR"
  echo "[INFO] $sample"

  # Optional flags
  strand_flag=()
  [[ -n "${LIB_STRAND:-}" ]] && strand_flag=( "--${LIB_STRAND}-stranded" )

  # Base args (index, output, reads)
  args=( -i "$IDX" -o "$OUTDIR" "${strand_flag[@]}" $BIAS $BOOT "$FQZ1" "$FQZ2" )

  # Attempt 1: as configured
  rc=$(run_kallisto "$THREADS" "$OUTDIR" "$log" "${args[@]}")
  if (( rc != 0 )) && (( RETRY_ON_FAIL )); then
    echo "[WARN] $sample: kallisto failed (exit $rc). Retrying without bootstraps…" | tee -a "$log"
    args_noboot=( -i "$IDX" -o "$OUTDIR" "${strand_flag[@]}" $BIAS "$FQZ1" "$FQZ2" )
    rc=$(run_kallisto "$THREADS" "$OUTDIR" "$log" "${args_noboot[@]}")
  fi
  if (( rc != 0 )) && (( RETRY_ON_FAIL )); then
    echo "[WARN] $sample: still failing. Retrying without bootstraps and without bias…" | tee -a "$log"
    args_nobbias=( -i "$IDX" -o "$OUTDIR" "${strand_flag[@]}" "$FQZ1" "$FQZ2" )
    rc=$(run_kallisto "$THREADS" "$OUTDIR" "$log" "${args_nobbias[@]}")
  fi
  if (( rc != 0 )) && (( RETRY_ON_FAIL )); then
    echo "[WARN] $sample: still failing. Retrying with low threads ($THREADS_FALLBACK) and no boots/bias…" | tee -a "$log"
    rc=$(run_kallisto "$THREADS_FALLBACK" "$OUTDIR" "$log" "${args_nobbias[@]}")
  fi

  if (( rc == 0 )) && [[ -s "$OUTDIR/abundance.tsv" && -s "$OUTDIR/run_info.json" ]]; then
    update_summary "$sample" "$OUTDIR" || true
    echo "$sample" >> "$SAMPLES_MANIFEST"
    echo "[INFO] done $sample"
  else
    echo "$sample (exit=$rc)" >> "$FAILED_LIST"
    echo "[ERROR] kallisto failed for $sample (exit $rc). See $log"
    tail -n 40 "$log" || true
    # continue to next sample
  fi

done < "$PAIRS_FILE"

echo "[INFO] Quantification complete."
echo "[INFO] Success: $(wc -l < "$SAMPLES_MANIFEST" | tr -d ' ') | Failed: $(wc -l < "$FAILED_LIST" | tr -d ' ')"
