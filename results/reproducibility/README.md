# Processed-data reproducibility audit, 16 November 2026

An empty working directory was populated from the review ZIP (SHA-256 `cd35e7a09c06e47deddf95cb0f968f24f2e264b112d842b0893f936273bc49f9`, revision based on commit c8b648c). The archived results directory was moved outside the analysis tree. Only the tximport object (SHA-256 `68a804c1cdb26aa791cfc95d3032f2b34e7ea012dfe9cb19edbc0fcdff1cb933`) and the gene annotation were copied back. Metadata, original published tables and frozen marker lists came from the bundle. The analysis used the local macOS R installation, R 4.6.1 and DESeq2 1.52.0. This tests computational reproducibility on the same machine starting at tximport. Raw-read processing and cross-platform portability were outside the audit.

The full blocked fit, group-only fit and all four block-omission fits were rebuilt. The joint three-degree-of-freedom test is in step 18. Step 21 checks the block assignments, model counts and 39 within-group permutations. Step 22 now reads fresh plain step-08 tables first and also accepts compressed archived tables.

| Quantity | Recalculated |
|---|---:|
| High-dose interaction discoveries | 990 |
| Direction consistent in all omission fits and individual blocks | 926 |
| Also significant in all four omission fits | 144 |
| Group-only discoveries | 89 |
| Strict discoveries also supported by the group-only model | 50 |
| Three-degree-of-freedom omnibus discoveries | 852 |

All 40 primary and supporting tables listed in [the manifest](table_comparison.tsv) matched byte for byte after decompression where applicable. Every row of manuscript Tables 1 and 2 matched at the displayed precision. The marker means in Additional file 1 matched after accounting for its explicit first-column header. The PNGs for manuscript Figures 1–3 had identical RGB pixels to the archived sources. The comparison covers the specified tables and figures, not every historical analysis in the repository.

At high dose the complement set exceeded its matched background by 0.457115 log2 units (BH-adjusted P = 0.000149993). The corresponding excesses were 0.005825 for ciliary markers (0.479776) and 0.105898 for homeostatic markers (0.376931). None of 39 permuted labelings reached the observed DEG count in either series (one-sided plus-one P = 0.025).

## Repeat the audit from the review bundle

Use a new extraction. Move the results directory to a separate reference directory, then copy back only results/tximport/txi.rds and results/tximport/gene_annotation.tsv. Run the following from the project root with the dependencies recorded in the session logs. Step 07 supplies the marker means. Steps 09, 12, 14 and 15 supply inputs for the first two manuscript figures. Live KEGG analyses in step 14 require network access.

```bash
Rscript --vanilla scripts/07_explore.R
Rscript --vanilla scripts/08_deseq2.R
Rscript --vanilla scripts/09_validate.R
Rscript --vanilla scripts/12_authors_model.R
Rscript --vanilla scripts/14_followup.R
Rscript --vanilla scripts/15_confidence.R
Rscript --vanilla scripts/16_figures.R
bash scripts/run_all.sh 17 23
python3 scripts/verify_results.py
```

The manifest contains SHA-256 digests of uncompressed table contents. SHA-256 for a compressed table can be obtained with `gzip -dc FILE.tsv.gz | shasum -a 256`. Use the fresh plain file when step 08 generated a TSV directly.

## Live pathway annotations

The nine GO tables from step 13 and all three residual pathway tables from step 14 matched their archived files exactly. The online KEGG section of step 13 was rerun separately on the fresh model output. All nine KEGG tables differed from the historical versions. For example, the mapped within-series ORA universe changed from 6,065 to 6,066 genes, and 25 down-enriched pathways passed the cutoff in the 500 nm series compared with 22 originally. This run used the live service and cannot establish which upstream changes account for every difference. The central complement result uses the installed rat GO annotation and the step-22 matched test. See the [pathway comparison](pathway_comparison.tsv) and [dated fresh KEGG results](kegg_2026-11-16/). To regenerate the exploratory pathway tables as a whole, also run step 13 with network access.
