#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT_DIR/data/recovered_raw_counts"
LOG_FILE="$OUT_DIR/monitor.log"
MERGED_H5AD="$OUT_DIR/GSE116222_raw_counts.h5ad"
STATE_FILE="$OUT_DIR/monitor.state"
INTERVAL_SECONDS="${1:-900}"

mkdir -p "$OUT_DIR"

status_line() {
  local now merged_size full_log_tail latest_triplet latest_triplet_name latest_triplet_size latest_triplet_mtime previous_size delta
  now="$(date '+%Y-%m-%d %H:%M:%S %Z')"
  full_log_tail="$(tail -n 1 "$OUT_DIR/full_recovery.log" 2>/dev/null || true)"
  latest_triplet="$(ls -1t "$OUT_DIR"/triplets/*.triplets.tsv 2>/dev/null | head -n 1 || true)"
  latest_triplet_name=""
  latest_triplet_size="0"
  latest_triplet_mtime=""
  delta="na"

  if [[ -n "$latest_triplet" ]]; then
    latest_triplet_name="$(basename "$latest_triplet")"
    latest_triplet_size="$(stat -f '%z' "$latest_triplet")"
    latest_triplet_mtime="$(stat -f '%Sm' -t '%Y-%m-%d %H:%M:%S %Z' "$latest_triplet")"
  fi

  previous_size="0"
  if [[ -f "$STATE_FILE" ]]; then
    previous_size="$(cat "$STATE_FILE" 2>/dev/null || echo 0)"
  fi

  if [[ -n "$latest_triplet" ]]; then
    delta="$(( latest_triplet_size - previous_size ))"
    printf '%s' "$latest_triplet_size" > "$STATE_FILE"
  fi

  if [[ -f "$MERGED_H5AD" ]]; then
    merged_size="$(ls -lh "$MERGED_H5AD" | awk '{print $5}')"
    printf '[%s] merged=yes merged_size=%s latest_triplet=%s latest_triplet_bytes=%s delta_bytes=%s latest_triplet_mtime="%s" tail="%s"\n' \
      "$now" "$merged_size" "$latest_triplet_name" "$latest_triplet_size" "$delta" "$latest_triplet_mtime" "$full_log_tail"
  else
    printf '[%s] merged=no latest_triplet=%s latest_triplet_bytes=%s delta_bytes=%s latest_triplet_mtime="%s" tail="%s"\n' \
      "$now" "$latest_triplet_name" "$latest_triplet_size" "$delta" "$latest_triplet_mtime" "$full_log_tail"
  fi
}

echo "# raw-count recovery monitor started $(date '+%Y-%m-%d %H:%M:%S %Z')" >> "$LOG_FILE"

while true; do
  status_line >> "$LOG_FILE"

  if [[ -f "$MERGED_H5AD" ]]; then
    echo "# raw-count recovery complete $(date '+%Y-%m-%d %H:%M:%S %Z')" >> "$LOG_FILE"
    exit 0
  fi

  sleep "$INTERVAL_SECONDS"
done
