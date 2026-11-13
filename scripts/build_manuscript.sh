#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
(cd manuscript && pandoc -f markdown+autolink_bare_uris Manuscript.md -o Manuscript.docx)
python3 scripts/format_manuscript.py
soffice --headless --convert-to pdf --outdir manuscript manuscript/Manuscript.docx
python3 - <<'PY'
from pathlib import Path
from zipfile import ZipFile
sources = ["fig3_volcano", "fig5_signatures", "fig7_dose_response"]
with ZipFile("manuscript/Figures.zip", "w") as archive:
    for i, name in enumerate(sources, 1):
        archive.write(Path("results/figures") / (name + ".pdf"), f"Figure{i}.pdf")
PY
