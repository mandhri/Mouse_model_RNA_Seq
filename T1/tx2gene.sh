#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---- CONFIG  ----
PROJ_DIR="/mnt/vol1/Mouse_model_RNA_Seq/T1"
INDEX_DIR="$PROJ_DIR/index"
GTF_GZ="$INDEX_DIR/Mus_musculus.GRCm39.115.gtf.gz"
GTF="$INDEX_DIR/Mus_musculus.GRCm39.115.gtf"
GTF_URL="https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz"
OUT_CSV="$PROJ_DIR/tx2gene_r115.csv"
RELEASE="115"

mkdir -p "$INDEX_DIR"

# --- Download GTF if missing ---
if [[ ! -f "$GTF" ]]; then
  echo "[INFO] Downloading Ensembl r115 GTF…"
  wget -O "$GTF_GZ" "$GTF_URL"
  gunzip -f "$GTF_GZ"
fi
echo "[INFO] Using GTF: $GTF"

# --- Build TXNAME,GENEID mapping from 'transcript' features (POSIX awk) ---
tmp="$INDEX_DIR/.tx2gene.tmp.csv"
awk 'BEGIN{FS=OFS="\t"}
  $3=="transcript"{
    tid=$9; gid=$9;
    sub(/.*transcript_id "/,"",tid); sub(/".*/,"",tid);
    sub(/.*gene_id "/,"",gid);       sub(/".*/,"",gid);
    if(tid!="" && gid!="") print tid","gid;
  }' "$GTF" | LC_ALL=C sort -u > "$tmp"

{ echo "TXNAME,GENEID"; cat "$tmp"; } > "$OUT_CSV"
rm -f "$tmp"
echo "[INFO] Done: $OUT_CSV"
