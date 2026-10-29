#!/usr/bin/env bash
# 01_download.sh — download paired FASTQ files from ENA with aria2 (parallel connections),
# survive network outages, resume partial files, and verify every file with md5.
#
# Usage (from the project root, inside the rnaseq env):
#   caffeinate -i bash scripts/01_download.sh SRR_Acc_List.txt data/raw
# Safe to rerun: verified files are skipped, partial files are resumed.
set -uo pipefail

ACC_LIST="${1:-SRR_Acc_List.txt}"
OUT="${2:-data/raw}"
JOBS="${JOBS:-3}"            # files downloaded at the same time
CONN="${CONN:-8}"            # connections per file
MAX_ROUNDS="${MAX_ROUNDS:-100}"  # how many times to restart aria2 after it gives up (network outages do not count)
WAIT="${WAIT:-60}"           # seconds to wait before restarting (e.g. while the network is down)
mkdir -p "$OUT"
LOG="$OUT/download.log"
LIST="$OUT/.aria2_input.txt"
MANIFEST="$OUT/.manifest.tsv"   # url<TAB>file<TAB>expected md5

command -v aria2c >/dev/null || { echo "❌ aria2c not found — run: conda activate rnaseq"; exit 1; }
md5_of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
is_done() { # $1 = file name, $2 = expected md5
  [ -f "$OUT/$1" ] && [ ! -f "$OUT/$1.aria2" ] && [ "$(md5_of "$OUT/$1")" = "$2" ]
}

# --- 1. Ask ENA once for every run's URLs + md5 (retry until all are known)
n_acc=$(grep -c '[^[:space:]]' "$ACC_LIST")
for round in $(seq 1 "$MAX_ROUNDS"); do
  : > "$MANIFEST.tmp"; failed=0
  while read -r ACC; do
    ACC=$(echo "$ACC" | tr -d '[:space:]'); [ -z "$ACC" ] && continue
    INFO=$(curl -s --retry 5 --retry-all-errors \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n 1)
    URLS=$(echo "$INFO" | cut -f2); MD5S=$(echo "$INFO" | cut -f3)
    if [ -z "$URLS" ] || [ "$URLS" = "fastq_ftp" ]; then failed=$((failed+1)); continue; fi
    IFS=';' read -ra U <<< "$URLS"; IFS=';' read -ra M <<< "$MD5S"
    for i in "${!U[@]}"; do printf 'https://%s\t%s\t%s\n' "${U[$i]}" "$(basename "${U[$i]}")" "${M[$i]}" >> "$MANIFEST.tmp"; done
  done < "$ACC_LIST"
  if [ "$failed" -eq 0 ]; then mv "$MANIFEST.tmp" "$MANIFEST"; break; fi
  log "ENA lookup failed for $failed of $n_acc runs (network?) — retrying in ${WAIT}s"
  sleep "$WAIT"
done
[ -s "$MANIFEST" ] || { log "❌ could not get file list from ENA — rerun later"; exit 1; }
log "manifest: $(wc -l < "$MANIFEST" | tr -d ' ') files"

# --- 2. Download rounds: aria2 gives up on some errors (e.g. DNS failure when the network drops),
#        so restart it until every file is verified
wait_for_network() { # block (without using up rounds) until the ENA file server answers
  local waited=0
  until curl -s -o /dev/null -m 20 -I "https://ftp.sra.ebi.ac.uk/vol1/fastq/"; do
    [ "$waited" -eq 0 ] && log "network/ENA unreachable — waiting until it is back (checking every ${WAIT}s)"
    waited=1; sleep "$WAIT"
  done
  [ "$waited" -eq 1 ] && log "network is back"
  return 0
}

for round in $(seq 1 "$MAX_ROUNDS"); do
  wait_for_network
  : > "$LIST"; todo=0
  while IFS=$'\t' read -r url name want; do
    is_done "$name" "$want" && continue
    printf '%s\n  dir=%s\n  out=%s\n' "$url" "$OUT" "$name" >> "$LIST"; todo=$((todo+1))
  done < "$MANIFEST"
  [ "$todo" -eq 0 ] && break
  log "round $round: $todo files to download"
  # -c : resume partial files            --max-tries=0 : retry forever on dropped connections
  # --async-dns=false : use the system DNS resolver (aria2's own resolver fails easily on macOS)
  aria2c -i "$LIST" -j "$JOBS" -x "$CONN" -s "$CONN" -c \
    --max-tries=0 --retry-wait=15 --timeout=60 --async-dns=false \
    --file-allocation=none --auto-file-renaming=false \
    --summary-interval=60 --console-log-level=warn \
    || { log "aria2 stopped with errors — waiting ${WAIT}s, then resuming"; sleep "$WAIT"; }
done

# --- 3. Verify every expected file with md5 (independent of aria2)
ok=0; bad=0
while IFS=$'\t' read -r url name want; do
  f="$OUT/$name"
  if [ ! -f "$f" ];                    then log "$name  MISSING"; bad=$((bad+1))
  elif [ -f "$f.aria2" ];              then log "$name  INCOMPLETE (rerun the script)"; bad=$((bad+1))
  elif [ "$(md5_of "$f")" = "$want" ]; then log "$name  OK"; ok=$((ok+1))
  else                                      log "$name  MD5_MISMATCH (delete it and rerun)"; bad=$((bad+1))
  fi
done < "$MANIFEST"
log "SUMMARY: $ok OK, $bad problems (expected $(wc -l < "$MANIFEST" | tr -d ' ') files)"
