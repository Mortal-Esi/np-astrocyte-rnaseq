#!/usr/bin/env bash
# run_all.sh — runs the whole analysis in order, from raw reads to figures.
#
# Usage (from the project root, inside the conda env "rnaseq"):
#   bash scripts/run_all.sh            # everything (download + quantification take many hours)
#   bash scripts/run_all.sh 06 16      # only steps 06 … 16 (the R analysis; needs results/salmon/ from step 05)
#
# Steps:
#   00 folder check            01 download FASTQ (ENA, md5-checked)   02 Ensembl 116 reference (GRCr8)
#   03 Salmon index            04 QC (FastQC, fastp, MultiQC)          05 Salmon quantification
#   06 tximport (gene counts)  07 / 07b exploration + replicate structure
#   08 DESeq2 (~ replicate + group)        09 comparison with the published tables
#   10 robustness (leave-one-replicate-out, filter)                    11 disagreements with the authors
#   12 reproducing the authors' model      13 GO / KEGG enrichment    14 reactive-astrocyte signatures, KEGG check
#   15 high-confidence DEGs               16 final figures
# Steps 12–14 download KEGG data (internet needed). Every step writes a log next to its results.

set -euo pipefail
cd "$(dirname "$0")/.."
from="${1:-00}"; to="${2:-16}"

steps=(00 01 02 03 04 05 06 07 07b 08 09 10 11 12 13 14 15 16)
cmd_of() {                               # plain case statement: macOS ships bash 3.2 (no associative arrays)
  case "$1" in
    00) echo "bash scripts/00_setup_check.sh" ;;      01) echo "bash scripts/01_download.sh" ;;
    02) echo "bash scripts/02_reference.sh go" ;;     03) echo "bash scripts/03_salmon_index.sh" ;;
    04) echo "bash scripts/04_qc.sh" ;;               05) echo "bash scripts/05_salmon_quant.sh" ;;
    06) echo "Rscript scripts/06_tximport.R" ;;       07) echo "Rscript scripts/07_explore.R" ;;
    07b) echo "Rscript scripts/07b_replicates.R" ;;   08) echo "Rscript scripts/08_deseq2.R" ;;
    09) echo "Rscript scripts/09_validate.R" ;;       10) echo "Rscript scripts/10_robustness.R" ;;
    11) echo "Rscript scripts/11_discordance.R" ;;    12) echo "Rscript scripts/12_authors_model.R" ;;
    13) echo "Rscript scripts/13_enrichment.R" ;;     14) echo "Rscript scripts/14_followup.R" ;;
    15) echo "Rscript scripts/15_confidence.R" ;;     16) echo "Rscript scripts/16_figures.R" ;;
  esac
}
num() { echo "${1//b/.5}"; }            # 07b sorts between 07 and 08

mkdir -p results/logs
for s in "${steps[@]}"; do
  awk -v x="$(num "$s")" -v a="$(num "$from")" -v b="$(num "$to")" 'BEGIN{exit !(x>=a && x<=b)}' || continue
  c=$(cmd_of "$s")
  echo "=== step $s: $c  ($(date '+%Y-%m-%d %H:%M'))"
  start=$(date +%s)
  $c 2>&1 | tee "results/logs/step_${s}.log"
  echo "=== step $s done in $(( $(date +%s) - start )) s"
done

# record exact software versions next to the results
Rscript -e 'suppressPackageStartupMessages(for (p in c("tximport", "DESeq2", "apeglm", "clusterProfiler", "org.Rn.eg.db", "enrichplot", "ggplot2")) library(p, character.only = TRUE)); writeLines(capture.output(sessionInfo()), "results/logs/R_sessionInfo.txt")'
conda env export --no-builds > results/logs/conda_environment.yml 2>/dev/null || echo "(conda env export skipped)"
echo "all requested steps finished"
