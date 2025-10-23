#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

#Check whether jq for processing JSON data is installed before start

if command -v jq >/dev/null 2>&1; then
  USE_JQ=1
  echo "[INFO] jq found: using jq for JSON parsing"
else
  USE_JQ=0
  echo "[INFO] jq not found: falling back to grep/cut/awk"
fi

RESULTS_DIR="/mnt/vol1/Mouse_model_RNA_Seq/kallisto_results"

# Outputs
ALIGN_SUM="${RESULTS_DIR}/alignment_summary.tsv"
SAMPLES_TXT="${RESULTS_DIR}/samples.txt"
SAMPLES_MANIFEST="${RESULTS_DIR}/samples_manifest.tsv"

echo -e "sample\tprocessed\tpseudoaligned\trate_percent\tkallisto_version" > "$ALIGN_SUM"
> "$SAMPLES_TXT"
echo -e "sample\tpath" > "$SAMPLES_MANIFEST"

for d in "$RESULTS_DIR"/*; do
  [ -d "$d" ] || continue
  info="$d/run_info.json"; [ -f "$info" ] || continue
  [ -f "$d/abundance.tsv" ] || continue

  sample=$(basename "$d")

  if [ "$USE_JQ" -eq 1 ]; then
    npr=$(jq -r '.n_processed' "$info")
    npa=$(jq -r '.n_pseudoaligned' "$info")
    rate=$(jq -r 'if .n_processed>0 then (.n_pseudoaligned/.n_processed*100) else 0 end' "$info")
    ver=$(jq -r '.kallisto_version // empty' "$info")
  else
    npr=$(grep -o '"n_processed":[0-9]*'      "$info" | cut -d: -f2)
    npa=$(grep -o '"n_pseudoaligned":[0-9]*'  "$info" | cut -d: -f2)
    # compute % only with awk if npr>0
    rate=$(awk -v a="${npa:-0}" -v b="${npr:-0}" 'BEGIN{printf "%.2f", (b>0? a/b*100 : 0)}')
    ver=$(grep -o '"kallisto_version":"[^"]*' "$info" | cut -d: -f2 | tr -d '"')
  fi

  echo -e "${sample}\t${npr}\t${npa}\t${rate}\t${ver}" >> "$ALIGN_SUM"
  echo -e "${sample}\t${d}/abundance.tsv" >> "$SAMPLES_MANIFEST"
  echo "$sample" >> "$SAMPLES_TXT"
done

echo "[INFO] wrote:"
echo "  $ALIGN_SUM"
echo "  $SAMPLES_MANIFEST"
echo "  $SAMPLES_TXT"




