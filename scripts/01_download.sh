#!/usr/bin/env bash
# 01_download.sh — resumable FASTQ download from ENA, with md5 verification
# Usage: bash scripts/01_download.sh SRR_Acc_List.txt data/raw
set -uo pipefail

ACC_LIST="${1:-SRR_Acc_List.txt}"
OUT="${2:-data/raw}"
mkdir -p "$OUT"
LOG="$OUT/download.log"

# md5 works differently on macOS and Linux
md5_of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

while read -r ACC; do
  ACC=$(echo "$ACC" | tr -d '[:space:]')
  [ -z "$ACC" ] && continue

  # Ask ENA for this run's file addresses and md5 checksums
  INFO=$(curl -s --retry 5 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \
    | tail -n 1)
  URLS=$(echo "$INFO" | cut -f2)
  MD5S=$(echo "$INFO" | cut -f3)

  if [ -z "$URLS" ]; then
    echo "$(date '+%F %T') $ACC  ENA_LOOKUP_FAILED" | tee -a "$LOG"; continue
  fi

  IFS=';' read -ra U <<< "$URLS"
  IFS=';' read -ra M <<< "$MD5S"

  for i in "${!U[@]}"; do
    URL="https://${U[$i]}"
    FILE="$OUT/$(basename "${U[$i]}")"
    WANT="${M[$i]}"

    # Already downloaded and verified → skip
    if [ -f "$FILE" ] && [ "$(md5_of "$FILE")" = "$WANT" ]; then
      echo "$(date '+%F %T') $(basename "$FILE")  OK (skipped)" | tee -a "$LOG"; continue
    fi

    # -C - : resume a partial file; --retry : retry when the connection drops
    curl -L -C - --retry 10 --retry-delay 15 -o "$FILE" "$URL"

    if [ "$(md5_of "$FILE")" = "$WANT" ]; then
      echo "$(date '+%F %T') $(basename "$FILE")  OK" | tee -a "$LOG"
    else
      echo "$(date '+%F %T') $(basename "$FILE")  MD5_MISMATCH (rerun the script)" | tee -a "$LOG"
    fi
  done
done < "$ACC_LIST"

echo "--- done. problems:"
grep -E "FAILED|MISMATCH" "$LOG" || echo "none"
