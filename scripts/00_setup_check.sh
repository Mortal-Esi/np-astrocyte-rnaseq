#!/usr/bin/env bash
# 00_setup_check.sh — build and verify the np-astrocyte-rnaseq project structure.
# Safe to run several times: it never deletes anything and never overwrites a correct file.
# Usage: bash 00_setup_check.sh

PROJ="$HOME/np-astrocyte-rnaseq"
DL="$HOME/Downloads"
ok=0; bad=0
pass() { echo "  ✅ $1"; ok=$((ok+1)); }
fail() { echo "  ❌ $1"; bad=$((bad+1)); }
info() { echo "  ➜ $1"; }

# Find a file in Downloads, also matching browser renames like "name (1).csv"
find_dl() { ls -t "$DL"/$1 2>/dev/null | head -1; }

echo "== 1. Folders"
for d in data/raw data/ref data/authors results scripts env; do
  mkdir -p "$PROJ/$d"
done
pass "folders ready under $PROJ"
cd "$PROJ" || exit 1

echo "== 2. Fix a misnamed samples.csv (SRA metadata saved under the wrong name)"
if [ -f samples.csv ] && head -1 samples.csv | grep -q '^Run,Assay Type'; then
  if [ ! -f SraRunTable.csv ]; then mv samples.csv SraRunTable.csv; info "renamed old samples.csv → SraRunTable.csv"
  else mv samples.csv samples_WRONG_backup.csv; info "moved wrong samples.csv → samples_WRONG_backup.csv"; fi
fi

echo "== 3. Bring files in from Downloads (only if missing here)"
bring() { # $1 = target path, $2 = pattern in Downloads
  if [ -f "$1" ]; then return; fi
  src=$(find_dl "$2")
  if [ -n "$src" ]; then mv "$src" "$1"; info "moved $(basename "$src") → $1"; fi
}
bring samples.csv            "samples*.csv"
bring SraRunTable.csv        "SraRunTable*.csv"
bring SRR_Acc_List.txt       "SRR_Acc_List*.txt"
bring scripts/01_download.sh "01_download*.sh"
for f in "$DL"/*GSE256038_DifferentialExpression_*.tsv.gz; do
  [ -f "$f" ] && mv "$f" data/authors/ && info "moved $(basename "$f") → data/authors/"
done

echo "== 4. .gitignore"
if [ ! -f .gitignore ]; then
  printf 'data/raw/\ndata/ref/\n*.fastq.gz\n*.fq.gz\n.DS_Store\n.Rhistory\n.RData\n' > .gitignore
  info "created .gitignore"
fi
grep -q '^data/raw/' .gitignore && pass ".gitignore excludes data/raw" || fail ".gitignore missing data/raw/"

echo "== 5. Check contents"
# samples.csv
if [ -f samples.csv ]; then
  hdr=$(head -1 samples.csv | tr -d '\r')
  n=$(grep -c . samples.csv)
  [ "$hdr" = "run,gsm,batch,size_nm,dose_ugml,group" ] && pass "samples.csv header" || fail "samples.csv header is wrong: $hdr"
  [ "$n" -eq 33 ] && pass "samples.csv has 32 samples" || fail "samples.csv has $((n-1)) samples (expected 32)"
  for g in 50A 50B 50C 50D 500A 500B 500C 500D; do
    c=$(cut -d, -f6 samples.csv | tr -d '\r' | grep -cx "$g")
    [ "$c" -eq 4 ] || fail "group $g has $c samples (expected 4)"
  done
else fail "samples.csv not found (download it from the chat)"; fi

# SraRunTable.csv
if [ -f SraRunTable.csv ]; then
  head -1 SraRunTable.csv | grep -q '^Run,' && pass "SraRunTable.csv present" || fail "SraRunTable.csv looks wrong"
else fail "SraRunTable.csv not found"; fi

# Accession list
if [ -f SRR_Acc_List.txt ]; then
  n=$(grep -c '^SRR' SRR_Acc_List.txt)
  [ "$n" -eq 32 ] && pass "SRR_Acc_List.txt has 32 runs" || fail "SRR_Acc_List.txt has $n runs (expected 32)"
  # every run in the list must be in samples.csv
  if [ -f samples.csv ]; then
    miss=$(tr -d '\r' < SRR_Acc_List.txt | grep '^SRR' | while read -r r; do grep -q "^$r," samples.csv || echo "$r"; done)
    [ -z "$miss" ] && pass "accession list matches samples.csv" || fail "runs missing from samples.csv: $miss"
  fi
else fail "SRR_Acc_List.txt not found"; fi

# Download script
if [ -f scripts/01_download.sh ]; then
  bash -n scripts/01_download.sh && pass "scripts/01_download.sh syntax OK" || fail "scripts/01_download.sh has a syntax error"
else fail "scripts/01_download.sh not found"; fi

# Authors' DE tables
n=$(ls data/authors/*_50.tsv.gz data/authors/*_500*.tsv.gz 2>/dev/null | wc -l | tr -d ' ')
[ "$n" -eq 2 ] && pass "both author DE tables in data/authors" || fail "found $n of 2 author DE tables in data/authors"

# Disk space
free_gb=$(df -g "$PROJ" 2>/dev/null | awk 'NR==2{print $4}')
[ -z "$free_gb" ] && free_gb=$(df -BG "$PROJ" | awk 'NR==2{gsub("G","",$4);print $4}')
[ "$free_gb" -ge 60 ] && pass "free disk: ${free_gb} GB" || fail "free disk only ${free_gb} GB (need ~60)"

echo
echo "== Result: $ok passed, $bad problems"
echo "== Tree:"
find . -not -path './.git*' -not -path './data/raw/*' | sort | sed 's|[^/]*/|  |g'
