"""Check key manuscript counts from the repository's published gene and pathway tables."""
import csv
import gzip
from pathlib import Path
root = Path(__file__).resolve().parents[1] / "results"

def read(name):
    path = root / name
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def significant(value):
    return value != "NA" and float(value) < 0.05

genes = read("interaction/interaction_genes.tsv.gz")
dose = read("dose_response/gene_tests.tsv.gz")
assert [g["gene_id"] for g in genes] == [g["id"] for g in dose]
assert sum(significant(g["omnibus_FDR"]) for g in dose) == 852
assert [sum(significant(g[f"dose{d}_FDR"]) for g in dose) for d in (1,10,100)] == [1,0,990]
assert sum(g["direction_stable"] == "TRUE" for g in genes) == 926
strict = [g for g in genes if g["strict_stable"] == "TRUE"]
assert len(strict) == 144
assert sum(significant(g["unblocked_padj"]) for g in genes) == 89
assert sum(significant(g["unblocked_padj"]) and g["unblocked_same_sign"] == "TRUE" for g in strict) == 50
sets = read("dose_sensitivity/pathway_stability.tsv")
for name, expected in [("ependymal_cilia", [0,5,6]), ("complement_activation", [0,0,6])]:
    for d, count in zip(("1","10","100"),expected):
        fits = [r for r in sets if r["ID"] == name and r["dose"] == d]
        assert len({r["fit"] for r in fits}) == 6
        assert sum(significant(r["FDR_six_tests"]) for r in fits) == count
print("Verified gene-level counts and all six dose-by-pathway summaries.")
