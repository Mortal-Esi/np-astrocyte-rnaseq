#!/usr/bin/env bash
# 05_salmon_quant.sh — quantify every sample with Salmon (selective alignment, paired-end).
# Output per sample: results/salmon/<run>/quant.sf (transcripts) and quant.genes.sf (genes, cross-check only).
#
# Usage (from the project root, inside the rnaseq env):
#   ONLY=SRR28006120 bash scripts/05_salmon_quant.sh     → one sample (test)
#   caffeinate -i bash scripts/05_salmon_quant.sh        → all samples
# Safe to rerun: samples that already have quant.sf are skipped.
set -uo pipefail

RAW="data/raw"
IDX="data/ref/salmon_index"
GTF="data/ref/Rattus_norvegicus.GRCr8.116.gtf.gz"
OUT="results/salmon"
THREADS="${THREADS:-4}"
ONLY="${ONLY:-}"
mkdir -p "$OUT"

command -v salmon >/dev/null || { echo "❌ salmon not found — run: conda activate rnaseq"; exit 1; }
[ -d "$IDX" ] || { echo "❌ index not found: $IDX"; exit 1; }
[ -f "$GTF" ] || { echo "❌ GTF not found: $GTF"; exit 1; }

runs=$(cut -d, -f1 samples.csv | tail -n +2)
[ -n "$ONLY" ] && runs="$ONLY"
n=$(echo "$runs" | wc -w | tr -d ' '); i=0

for r in $runs; do
  i=$((i+1))
  if [ -f "$OUT/$r/quant.sf" ]; then echo "[$i/$n] $r already done — skipped"; continue; fi
  echo "[$i/$n] $r  started $(date '+%H:%M:%S')"
  # -l A          : detect library type automatically (reported in lib_format_counts.json)
  # --seqBias     : correct random-hexamer priming bias (the FastQC "per base sequence content" failure)
  # --gcBias      : correct fragment GC bias
  # -g + --ignoreTxVersion : also write gene-level estimates (versions differ between FASTA and GTF)
  salmon quant \
    -i "$IDX" \
    -l A \
    -1 "$RAW/${r}_1.fastq.gz" \
    -2 "$RAW/${r}_2.fastq.gz" \
    -p "$THREADS" \
    --seqBias \
    --gcBias \
    -g "$GTF" \
    --ignoreTxVersion \
    -o "$OUT/$r" \
    > "$OUT/${r}.log" 2>&1
  status=$?
  if [ $status -ne 0 ] || [ ! -f "$OUT/$r/quant.sf" ]; then
    echo "   ❌ failed (exit $status) — last lines of $OUT/${r}.log:"; tail -15 "$OUT/${r}.log"; exit 1
  fi
  echo "   ✅ finished $(date '+%H:%M:%S')"
done

# Summary table: reads processed, mapping rate, detected library type, mean fragment length
python3 - "$OUT" <<'PY'
import json, os, sys, glob
out = sys.argv[1]
def find(d, keys):
    """search a nested dict for the first key from `keys`"""
    if isinstance(d, dict):
        for k in keys:
            if k in d: return d[k]
        for v in d.values():
            x = find(v, keys)
            if x is not None: return x
    return None
rows = []
for meta in sorted(glob.glob(f"{out}/*/aux_info/meta_info.json")):
    s = meta.split("/")[-3]
    m = json.load(open(meta))
    lib = {}
    lf = f"{out}/{s}/lib_format_counts.json"
    if os.path.exists(lf): lib = json.load(open(lf))
    fld = find(m, ["frag_length_mean", "fragment_length_mean"])
    rows.append([s, find(m, ["num_processed"]), find(m, ["num_mapped"]), find(m, ["percent_mapped"]),
                 find(lib, ["expected_format"]) or find(m, ["library_types"]), fld])
print("\nsample\tprocessed\tmapped\tpercent_mapped\tlib_type\tfrag_len_mean")
for r in rows: print("\t".join("" if x is None else str(x) for x in r))
with open(f"{out}/salmon_summary.tsv", "w") as f:
    f.write("sample\tprocessed\tmapped\tpercent_mapped\tlib_type\tfrag_len_mean\n")
    for r in rows: f.write("\t".join("" if x is None else str(x) for x in r) + "\n")
print(f"\nsaved: {out}/salmon_summary.tsv")
PY
