
#!/usr/bin/env bash
set -euo pipefail
# >>> activate project env on /mnt/vol1 >>>
CONDA_DIR="/mnt/vol1/Mouse_model_RNA_seq/software/mambaforge"
source "$CONDA_DIR/etc/profile.d/conda.sh"
conda activate "/mnt/vol1/Mouse_model_RNA_seq/software/envs/rnaseq-mouse"
# <<< env activated <<<

# Usage: ./qc_one_sample.sh BB329 [threads]
SAMPLE="${1:?Usage: $0 <SAMPLE_ID> [THREADS]}"
THREADS="${2:-8}"

IN="/mnt/vol1/Mouse_model_RNA_seq/Mouse/tmp/${SAMPLE}"
OUT_BASE="/mnt/vol1/Mouse_model_RNA_seq"
FASTQC_OUT="${OUT_BASE}/qc/fastqc/${SAMPLE}"
MULTIQC_OUT="${OUT_BASE}/qc/multiqc"

# --- sanity checks ---
command -v fastqc >/dev/null 2>&1 || { echo "ERROR: fastqc not found in PATH"; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc not found in PATH"; exit 1; }
[ -d "$IN" ] || { echo "ERROR: input folder not found: $IN"; exit 1; }

echo "INFO: Sample: $SAMPLE"
echo "INFO: Input : $IN"
echo "INFO: Out   : FASTQC -> $FASTQC_OUT ; MULTIQC -> $MULTIQC_OUT"

# list lanes/reads we have
echo "INFO: Lane/read coverage:"
ls -1 "${IN}/${SAMPLE}"*_L0??_R?.fastq.gz | wc -l
ls -1 "${IN}/${SAMPLE}"*_L0??_R?.fastq.gz \
 | sed -n 's/.*_L\([0-9]\{3\}\)_R\([12]\).*/L\1 R\2/p' \
 | sort | uniq -c

# make output dirs
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

# run FastQC on ALL lanes (R1+R2)
echo "INFO: running FastQC → $FASTQC_OUT"
fastqc -t "$THREADS" -o "$FASTQC_OUT" "${IN}/${SAMPLE}"*_L0??_R?.fastq.gz

# aggregate with MultiQC (one HTML per sample)
DATESTAMP="$(date +%Y%m%d_%H%M%S)"
REPORT_NAME="multiqc_raw_${SAMPLE}_${DATESTAMP}"
echo "INFO: running MultiQC → ${MULTIQC_OUT}/${REPORT_NAME}.html"
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT" -n "$REPORT_NAME"

echo "DONE. Open in RStudio Files:"
echo "  ${MULTIQC_OUT}/${REPORT_NAME}.html"
