#!/usr/bin/env bash
# 02_reference.sh — download the rat reference (transcript sequences + annotation) from Ensembl.
# Pins ONE Ensembl release so FASTA and GTF always match.
# Usage:  bash scripts/02_reference.sh          → show what would be downloaded
#         bash scripts/02_reference.sh go       → download
set -uo pipefail

OUT="data/ref"
mkdir -p "$OUT"
SPECIES="rattus_norvegicus"

# 1. Current Ensembl release number (e.g. 115)
REL=$(curl -s --retry 5 "https://rest.ensembl.org/info/data?content-type=application/json" | grep -oE '[0-9]+' | head -1)
[ -z "$REL" ] && { echo "❌ Could not reach Ensembl REST API. Try with VPN."; exit 1; }
BASE="https://ftp.ensembl.org/pub/release-${REL}"

# 2. Find exact file names by listing the server folders
pick() { # $1 = folder, $2 = regex for the file name
  curl -s --retry 5 "$BASE/$1/" | grep -oE 'href="[^"]+"' | cut -d'"' -f2 | grep -E "$2" | head -1
}
CDNA=$(pick "fasta/$SPECIES/cdna"  '\.cdna\.all\.fa\.gz$')
NCRNA=$(pick "fasta/$SPECIES/ncrna" '\.ncrna\.fa\.gz$')
GTF=$(pick "gtf/$SPECIES"          "\.${REL}\.gtf\.gz$")

echo "Ensembl release : $REL"
echo "cDNA  (protein-coding transcripts) : ${CDNA:-NOT FOUND}"
echo "ncRNA (non-coding transcripts)     : ${NCRNA:-NOT FOUND}"
echo "GTF   (transcript → gene map)      : ${GTF:-NOT FOUND}"
if [ -z "$CDNA" ] || [ -z "$NCRNA" ] || [ -z "$GTF" ]; then echo "❌ A file was not found — send this output to Claude."; exit 1; fi

if [ "${1:-}" != "go" ]; then echo; echo "Looks right? Run again with:  bash scripts/02_reference.sh go"; exit 0; fi

# 3. Download (resumable) and test each file is a complete gzip
for pair in "fasta/$SPECIES/cdna/$CDNA" "fasta/$SPECIES/ncrna/$NCRNA" "gtf/$SPECIES/$GTF"; do
  f="$OUT/$(basename "$pair")"
  curl -L -C - --retry 10 --retry-delay 15 -o "$f" "$BASE/$pair"
  if gzip -t "$f" 2>/dev/null; then echo "✅ $(basename "$f")  $(du -h "$f" | cut -f1)"
  else echo "❌ $(basename "$f") is incomplete — rerun the script"; fi
done

# 4. Record exactly what was used (goes into the README later)
{ echo "ensembl_release: $REL"; echo "cdna: $CDNA"; echo "ncrna: $NCRNA"; echo "gtf: $GTF"; echo "downloaded: $(date '+%F')"; } > "$OUT/REFERENCE_VERSION.txt"
echo; cat "$OUT/REFERENCE_VERSION.txt"
