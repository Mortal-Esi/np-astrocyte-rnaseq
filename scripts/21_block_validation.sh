#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=results/block_validation
mkdir -p "$out"
python3 scripts/biosample_map.py
Rscript --vanilla scripts/compare_counts.R . "$out"
Rscript --vanilla scripts/block_control.R . "$out"
python3 scripts/summarize_permutations.py
