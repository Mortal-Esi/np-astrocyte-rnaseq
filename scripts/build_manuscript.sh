#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
(cd manuscript && pandoc -f markdown+autolink_bare_uris Manuscript.md -o Manuscript.docx)
python3 scripts/format_manuscript.py
soffice --headless --convert-to pdf --outdir manuscript manuscript/Manuscript.docx manuscript/Cover_letter.docx
python3 - <<'PY'
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
sources = ["fig3_volcano", "fig5_signatures", "fig7_dose_response"]
with ZipFile("manuscript/Figures.zip", "w") as archive:
    for i, name in enumerate(sources, 1):
        archive.write(Path("results/figures") / (name + ".pdf"), f"Figure{i}.pdf")

supplement = [
    "data/marker_sources/SOURCES.md",
    "data/marker_sources/sets_source.tsv",
    "data/marker_sources/liddelow_audit.tsv",
    "results/marker_sensitivity/external_gsea.tsv",
    "results/amplitude/external_matched.tsv",
    "results/amplitude/set_excess_over_amplitude.tsv",
    "results/amplitude/complement_genes.tsv",
]
with ZipFile("manuscript/supplement/Additional_file_2.zip", "w", ZIP_DEFLATED) as archive:
    for source in map(Path, supplement):
        archive.write(source, source.name)
PY
