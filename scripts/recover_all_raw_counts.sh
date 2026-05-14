#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="/opt/homebrew/Caskroom/miniconda/base/envs/scrna-uc/bin/python"
REF_H5AD="$ROOT_DIR/data/processed/GSE116222_raw.h5ad"
MANIFEST="$ROOT_DIR/data/recovered_raw_counts/sra_manifest.tsv"
OUT_DIR="$ROOT_DIR/data/recovered_raw_counts"
MERGED_H5AD="$OUT_DIR/GSE116222_raw_counts.h5ad"

mkdir -p "$OUT_DIR/triplets" "$OUT_DIR/counts" "$OUT_DIR/h5ad" "$OUT_DIR/logs" "$OUT_DIR/tmp"

emit_triplets_with_retries() {
  local sample="$1"
  local srr="$2"
  local bam_url="$3"
  local triplets="$4"
  local triplets_done="$5"
  local counts="$6"
  local counts_done="$7"
  local sample_h5ad="$8"
  local h5ad_done="$9"
  local log_file="${10}"
  local attempt

  for attempt in 1 2 3; do
    rm -f "$triplets" "$triplets_done" "$counts" "$counts_done" "$sample_h5ad" "$h5ad_done" "$log_file"
    echo "[recover] emit-triplets $sample ($srr) attempt=$attempt"
    if "$PYTHON_BIN" "$ROOT_DIR/scripts/recover_raw_counts_from_bam.py" emit-triplets \
      --reference-h5ad "$REF_H5AD" \
      --sample "$sample" \
      --bam-url "$bam_url" \
      --progress-every 5000000 \
      > "$triplets" 2> "$log_file"; then
      touch "$triplets_done"
      return 0
    fi
  done

  echo "[recover] emit-triplets failed after 3 attempts for $sample" >&2
  return 1
}

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample srr bam_url; do
  triplets="$OUT_DIR/triplets/${sample}.triplets.tsv"
  counts="$OUT_DIR/counts/${sample}.counts.tsv"
  sample_h5ad="$OUT_DIR/h5ad/${sample}_raw_counts.h5ad"
  log_file="$OUT_DIR/logs/${sample}.emit.log"
  triplets_done="$OUT_DIR/triplets/${sample}.done"
  counts_done="$OUT_DIR/counts/${sample}.done"
  h5ad_done="$OUT_DIR/h5ad/${sample}.done"

  if [[ -s "$sample_h5ad" && -f "$h5ad_done" ]]; then
    echo "[recover] reuse h5ad $sample"
    continue
  fi

  if [[ ! -s "$triplets" || ! -f "$triplets_done" ]]; then
    emit_triplets_with_retries "$sample" "$srr" "$bam_url" "$triplets" "$triplets_done" "$counts" "$counts_done" "$sample_h5ad" "$h5ad_done" "$log_file"
  else
    echo "[recover] reuse triplets $sample"
  fi

  if [[ ! -s "$counts" || ! -f "$counts_done" ]]; then
    rm -f "$counts" "$counts_done" "$sample_h5ad" "$h5ad_done"
    echo "[recover] collapse-umis $sample"
    LC_ALL=C sort -u -T "$OUT_DIR/tmp" "$triplets" | cut -f1,2 | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$1}' > "$counts"
    touch "$counts_done"
  else
    echo "[recover] reuse counts $sample"
  fi

  if [[ ! -s "$sample_h5ad" || ! -f "$h5ad_done" ]]; then
    rm -f "$sample_h5ad" "$h5ad_done"
    echo "[recover] write-h5ad $sample"
    "$PYTHON_BIN" "$ROOT_DIR/scripts/recover_raw_counts_from_bam.py" counts-to-h5ad \
      --reference-h5ad "$REF_H5AD" \
      --sample "$sample" \
      --counts-tsv "$counts" \
      --output-h5ad "$sample_h5ad"
    touch "$h5ad_done"
  else
    echo "[recover] reuse h5ad $sample"
  fi
done

sample_h5ads=("$OUT_DIR"/h5ad/*_raw_counts.h5ad)

echo "[recover] merge ${#sample_h5ads[@]} sample h5ad files"
"$PYTHON_BIN" "$ROOT_DIR/scripts/recover_raw_counts_from_bam.py" merge-h5ad \
  --input-h5ads "${sample_h5ads[@]}" \
  --output-h5ad "$MERGED_H5AD"

echo "[recover] done -> $MERGED_H5AD"
