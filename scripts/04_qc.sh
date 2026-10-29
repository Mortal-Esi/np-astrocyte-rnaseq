#!/usr/bin/env bash
# 04_qc.sh — quality control of the raw FASTQ files (no reads are modified).
#   FastQC : standard per-file report
#   fastp  : per-sample paired-end report (adapters, duplication, insert size) — report only, no output reads
#   MultiQC: one combined report for all samples and both tools
# Usage (from the project root, inside the rnaseq env):
#   caffeinate -i bash scripts/04_qc.sh
# Safe to rerun: samples that already have reports are skipped.
set -uo pipefail

RAW="data/raw"
QC="results/qc"
THREADS="${THREADS:-4}"
mkdir -p "$QC/fastqc" "$QC/fastp" "$QC/multiqc"

for t in fastqc fastp multiqc; do command -v $t >/dev/null || { echo "❌ $t not found — run: conda activate rnaseq"; exit 1; }; done

runs=$(cut -d, -f1 samples.csv | tail -n +2)
n=$(echo "$runs" | wc -l | tr -d ' ')
i=0
for r in $runs; do
  i=$((i+1))
  r1="$RAW/${r}_1.fastq.gz"; r2="$RAW/${r}_2.fastq.gz"
  [ -f "$r1" ] && [ -f "$r2" ] || { echo "❌ missing FASTQ for $r"; exit 1; }

  # FastQC (one report per file)
  if [ ! -f "$QC/fastqc/${r}_2_fastqc.zip" ]; then
    echo "[$i/$n] FastQC  $r"
    fastqc -q -t "$THREADS" --outdir "$QC/fastqc" "$r1" "$r2"
  fi

  # fastp in report-only mode: no -o/-O, so no trimmed reads are written
  if [ ! -f "$QC/fastp/${r}.json" ]; then
    echo "[$i/$n] fastp   $r"
    fastp -i "$r1" -I "$r2" \
      --detect_adapter_for_pe \
      --overrepresentation_analysis \
      -w "$THREADS" \
      -j "$QC/fastp/${r}.json" -h "$QC/fastp/${r}.html" \
      2> "$QC/fastp/${r}.log"
  fi
done

echo "MultiQC → $QC/multiqc"
multiqc -f -q -o "$QC/multiqc" "$QC/fastqc" "$QC/fastp"
echo "✅ done — open: $QC/multiqc/multiqc_report.html"
