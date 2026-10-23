#!/usr/bin/env bash
# 03_salmon_index.sh — build the Salmon index from the rat transcriptome (cDNA + ncRNA).
# No decoys (genome decoys need >16 GB RAM; this Mac has 8 GB).
# Usage (inside the rnaseq env, from the project root):  bash scripts/03_salmon_index.sh
set -euo pipefail

REF="data/ref"
IDX="$REF/salmon_index"
CDNA="$REF/Rattus_norvegicus.GRCr8.cdna.all.fa.gz"
NCRNA="$REF/Rattus_norvegicus.GRCr8.ncrna.fa.gz"

command -v salmon >/dev/null || { echo "❌ salmon not found — run: conda activate rnaseq"; exit 1; }
for f in "$CDNA" "$NCRNA"; do [ -f "$f" ] || { echo "❌ missing $f"; exit 1; }; done
[ -d "$IDX" ] && { echo "❌ $IDX already exists. Delete it first if you want to rebuild: rm -rf $IDX"; exit 1; }

echo "salmon $(salmon --version | awk '{print $2}') — building index at $IDX"
START=$(date +%s)

# -t  : both FASTA files are indexed together (no need to merge them first)
# -k 31 : k-mer length (default; standard for reads of ~75 bp)
# -p 4  : 4 threads (M1 has 8 cores; leave some for the system)
# --ramLimit 4 : cap index-build memory at 4 GiB (default 8 would fill this Mac's RAM)
salmon index \
  -t "$CDNA" "$NCRNA" \
  -i "$IDX" \
  -k 31 \
  -p 4 \
  --ramLimit 4 \
  2>&1 | tee "$REF/salmon_index.log"

END=$(date +%s)
echo
echo "✅ done in $(( (END-START)/60 )) min — index size: $(du -sh "$IDX" | cut -f1)"
echo "salmon_version: $(salmon --version | awk '{print $2}')" >> "$REF/REFERENCE_VERSION.txt"
echo "salmon_index: k=31, no decoys, cdna+ncrna" >> "$REF/REFERENCE_VERSION.txt"
