"""Rebuild the study-specific GEO / BioSample / SRA crosswalk from archived records."""
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

root = Path(__file__).resolve().parents[1]
records = {b.get("accession"): b for b in ET.parse(root / "data/biosamples.xml").getroot()}
with (root / "SraRunTable.csv").open() as f:
    biosamples = {row["Run"]: row["BioSample"] for row in csv.DictReader(f)}
rows = []
with (root / "samples.csv").open() as f:
    for sample in csv.DictReader(f):
        run = sample["run"]
        accession = biosamples[run]
        record = records[accession]
        title = record.findtext("Description/Title")
        _, group, suffix = title.split("-")
        gsm = record.findtext("Ids/Id[@db='GEO']")
        if (gsm, group, "rep" + str(int(suffix))) != (sample["gsm"], sample["group"], sample["replicate"]):
            raise ValueError(f"Inconsistent identifiers for {run}")
        rows.append([gsm, run, group, sample["replicate"], accession, title, suffix,
                     "https://www.ncbi.nlm.nih.gov/biosample/" + accession])
if len({row[4] for row in rows}) != len(records) or len(rows) != 32:
    raise ValueError("Expected 32 unique matched BioSamples")
with (root / "results/block_validation/biosample_map.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["GSM","SRR","group","inferred_block","BioSample","submitted_title","suffix","source"])
    writer.writerows(rows)
