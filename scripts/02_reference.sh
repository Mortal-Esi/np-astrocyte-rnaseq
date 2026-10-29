#!/usr/bin/env bash
# 02_reference.sh — download the rat reference from Ensembl, pinned to ONE release
# so that every reference file comes from the same assembly and annotation.
#
# Usage (from the project root):
#   bash scripts/02_reference.sh              → show transcriptome files (no download)
#   bash scripts/02_reference.sh go           → download transcriptome (cDNA + ncRNA + GTF)
#   bash scripts/02_reference.sh genome       → show genome file (no download)
#   bash scripts/02_reference.sh genome go    → download genome (needed only for a decoy-aware Salmon index)
set -uo pipefail

OUT="data/ref"
mkdir -p "$OUT"
SPECIES="rattus_norvegicus"
VERFILE="$OUT/REFERENCE_VERSION.txt"

MODE="transcriptome"; GO="no"
for a in "$@"; do
  case "$a" in
    genome) MODE="genome" ;;
    go)     GO="yes" ;;
    *) echo "❌ unknown option: $a"; exit 1 ;;
  esac
done

# --- 1. Ensembl release: reuse the one already recorded, so all files stay on the same release
if [ -f "$VERFILE" ] && grep -q '^ensembl_release:' "$VERFILE"; then
  REL=$(grep '^ensembl_release:' "$VERFILE" | awk '{print $2}')
  echo "Using release recorded in $VERFILE: $REL"
else
  REL=$(curl -s --retry 5 "https://rest.ensembl.org/info/data?content-type=application/json" | grep -oE '[0-9]+' | head -1)
  [ -z "$REL" ] && { echo "❌ Could not reach Ensembl REST API. Try with VPN."; exit 1; }
  echo "Current Ensembl release: $REL"
fi
BASE="https://ftp.ensembl.org/pub/release-${REL}"

# --- helpers
pick() { # $1 = server folder, $2 = regex for the file name → prints the first matching file name
  curl -s --retry 5 "$BASE/$1/" | grep -oE 'href="[^"]+"' | cut -d'"' -f2 | grep -E "$2" | head -1
}
set_kv() { # $1 = key, $2 = value → add or replace one line in the version file
  touch "$VERFILE"
  grep -v "^$1:" "$VERFILE" > "$VERFILE.tmp" || true
  echo "$1: $2" >> "$VERFILE.tmp"
  mv "$VERFILE.tmp" "$VERFILE"
}
fetch() { # $1 = path on server → download into data/ref, resumable, verified with gzip -t
  local f="$OUT/$(basename "$1")"
  if [ -f "$f" ] && gzip -t "$f" 2>/dev/null; then echo "✅ $(basename "$f") already complete — skipped"; return; fi
  curl -L -C - --retry 10 --retry-delay 15 -o "$f" "$BASE/$1"
  if gzip -t "$f" 2>/dev/null; then echo "✅ $(basename "$f")  $(du -h "$f" | cut -f1)"
  else echo "❌ $(basename "$f") is incomplete — rerun the same command"; fi
}

# --- 2. Find files on the server
if [ "$MODE" = "transcriptome" ]; then
  CDNA=$(pick  "fasta/$SPECIES/cdna"  '\.cdna\.all\.fa\.gz$')
  NCRNA=$(pick "fasta/$SPECIES/ncrna" '\.ncrna\.fa\.gz$')
  GTF=$(pick   "gtf/$SPECIES"         "\.${REL}\.gtf\.gz$")
  echo "cDNA  (coding + pseudogene transcripts) : ${CDNA:-NOT FOUND}"
  echo "ncRNA (non-coding transcripts)          : ${NCRNA:-NOT FOUND}"
  echo "GTF   (transcript → gene map)           : ${GTF:-NOT FOUND}"
  if [ -z "$CDNA" ] || [ -z "$NCRNA" ] || [ -z "$GTF" ]; then echo "❌ A file was not found — send this output to Claude."; exit 1; fi
else
  # primary_assembly = chromosomes + unplaced scaffolds, without alternative haplotypes (preferred);
  # toplevel is the fallback when a species has no primary_assembly file.
  GENOME=$(pick "fasta/$SPECIES/dna" '\.dna\.primary_assembly\.fa\.gz$')
  [ -z "$GENOME" ] && GENOME=$(pick "fasta/$SPECIES/dna" '\.dna\.toplevel\.fa\.gz$')
  echo "Genome (DNA of all chromosomes) : ${GENOME:-NOT FOUND}"
  [ -z "$GENOME" ] && { echo "❌ Genome file not found — send this output to Claude."; exit 1; }
  echo "(about 0.8 GB compressed)"
fi

if [ "$GO" != "yes" ]; then
  echo; echo "Looks right? Run again with 'go' added to download."; exit 0
fi

# --- 3. Download and record versions
set_kv ensembl_release "$REL"
if [ "$MODE" = "transcriptome" ]; then
  fetch "fasta/$SPECIES/cdna/$CDNA";   set_kv cdna "$CDNA"
  fetch "fasta/$SPECIES/ncrna/$NCRNA"; set_kv ncrna "$NCRNA"
  fetch "gtf/$SPECIES/$GTF";           set_kv gtf "$GTF"
  set_kv transcriptome_downloaded "$(date '+%F')"
else
  fetch "fasta/$SPECIES/dna/$GENOME";  set_kv genome "$GENOME"
  set_kv genome_downloaded "$(date '+%F')"
fi
echo; cat "$VERFILE"
