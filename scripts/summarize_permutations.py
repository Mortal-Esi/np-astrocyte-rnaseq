"""Recreate the summary from per-draw counts and recorded sample assignments."""
import csv
from pathlib import Path
from statistics import median

out = Path(__file__).resolve().parents[1] / "results/block_validation"

def read(name):
    with (out / name).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

draws = read("permuted_counts.tsv")
baseline = {row["model"]: row for row in read("control_baselines.tsv")}
groups = {}
for row in read("permuted_labels.tsv"):
    groups.setdefault(row["group"], []).append(row)
if len(draws) != 39 or [int(row["draw"]) for row in draws] != list(range(1, 40)):
    raise ValueError("Expected 39 completed draws in order")
for draw in draws:
    for samples in groups.values():
        if sorted(s["X" + draw["draw"]] for s in samples) != ["rep1", "rep2", "rep3", "rep4"]:
            raise ValueError(f"Unbalanced labels in draw {draw['draw']}")
rows = []
for size in ("50", "500"):
    counts = [int(row["n" + size]) for row in draws]
    observed = int(baseline["inferred"][size])
    exceed = sum(n >= observed for n in counts)
    rows.append([size,baseline["unblocked"][size],observed,len(counts),min(counts),
                 median(counts),max(counts),exceed,(exceed+1)/(len(counts)+1)])
with (out / "permutation_summary.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["size_nm","unblocked","observed","draws","min","median","max",
                     "at_least_observed","one_sided_plus_one"])
    writer.writerows(rows)
