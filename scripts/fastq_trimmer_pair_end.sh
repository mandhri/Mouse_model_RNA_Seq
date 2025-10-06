#!/usr/bin/env bash
set -euo pipefail

# >>> activate project env on /mnt/vol1 >>>
CONDA_DIR="/mnt/vol1/Mouse_model_RNA_seq/miniforge3"
source /mnt/vol1/Mouse_model_RNA_Seq/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/vol1/Mouse_model_RNA_Seq/software/envs/rnaseq-mouse

# <<< env activated <<<


# Usage: ./trim_pe_skewer.sh <SAMPLE_ID> [THREADS]
SAMPLE="${1:?Usage: $0 <SAMPLE_ID> [THREADS]}"
THREADS="${2:-16}"

BASE="/mnt/vol1/Mouse_model_RNA_seq/"
IN="${BASE}/tmp/${SAMPLE}"                     # where the raw .fastq.gz live
TRIM="${BASE}/trim/${SAMPLE}"                  # trimmed outputs
QCR="${BASE}/qc/fastqc/trimmed/${SAMPLE}"      # FastQC for trimmed reads
MQC="${BASE}/qc/multiqc"                       # MultiQC folder

#Thresholds 
QCUT=20        # trim bases below Phred 20
LCUT=31        # discard reads shorter than 31 nt
FCUT=5         # (kept from your script; skewer will ignore if not applicable to your version)

mkdir -p "$TRIM" "$QCR" "$MQC"

echo "== Trimming ${SAMPLE} per-lane with Skewer (PE) =="

for L in 001 002 003 004 005 006 007; do
  R1=("${IN}/${SAMPLE}"*_L${L}_R1.fastq.gz)
  R2=("${IN}/${SAMPLE}"*_L${L}_R2.fastq.gz)

  # Skip this lane if it doesn't exist
  [[ -e "${R1[0]}" && -e "${R2[0]}" ]] || { echo "  L${L}: no files, skipping"; continue; }

  OUTPFX="${TRIM}/${SAMPLE}_L${L}"

  echo "  L${L}: ${R1[0]##*/} + ${R2[0]##*/}"

  # Skewer paired-end mode: auto-detect Illumina adapters, trim by quality/QCUT, min length LCUT
  # Note: some skewer builds don't use -f the way your old script did; safe to omit if unknown.
  skewer \
    -m pe \
    -q ${QCUT} \
    -l ${LCUT} \
    -t ${THREADS} \
    -o "${OUTPFX}" \
    "${R1[0]}" "${R2[0]}"

  # Skewer outputs:
  #   ${OUTPFX}-trimmed-pair1.fastq
  #   ${OUTPFX}-trimmed-pair2.fastq
  # Gzip them (FastQC reads .gz fine; saves space/IO)
  gzip -f "${OUTPFX}-trimmed-pair1.fastq"
  gzip -f "${OUTPFX}-trimmed-pair2.fastq"
done

echo "== FastQC on trimmed reads =="
mkdir -p "$QCR"
fastqc -t "$THREADS" -o "$QCR" "${TRIM}/${SAMPLE}"*_L???_trimmed-pair?.fastq.gz

echo "== MultiQC (trimmed only) =="
STAMP=$(date +%Y%m%d_%H%M%S)
multiqc "$QCR" -o "$MQC" -n "multiqc_trimmed_${SAMPLE}_${STAMP}"

echo "DONE."
echo "Open MultiQC: ${MQC}/multiqc_trimmed_${SAMPLE}_${STAMP}.html"
echo "Trimmed FASTQs under: ${TRIM}"
