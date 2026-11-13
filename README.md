# Reanalysis of polystyrene nano- and microplastic exposure in rat astrocytes (RNA-seq)

**Mohammad Javad Eslami** · independent reanalysis of public data (GEO [GSE256038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256038) / BioProject PRJNA1077930)

Rat neural-stem-cell-derived astrocytes were exposed for seven days to 50 nm or 500 nm polystyrene at equal nominal mass concentrations (0, 1, 10 and 100 µg/mL). I reprocessed all 32 RNA-seq libraries. Both series show a broad response at 100 µg/mL. Four samples per group have suffixes 01–04 in the submitted BioSample titles. A model using those suffixes as blocks captures cross-group expression structure. The experimental origin of that structure remains unresolved. Direct comparisons of the two responses at equal mass concentration identify candidate differences in homeostatic and ciliary programmes.

## Key findings

1. **Dose response.** The high-dose samples separate along PC1. Across the four comparisons at 1 and 10 µg/mL, the inferred-block model calls 0–39 genes at adjusted P < 0.05. At 100 µg/mL, both series show thousands of changes.
2. **Sample structure.** Each suffix occurs once per group and matches the within-group GSM rank for all 32 libraries. The suffix explains 78% and 91% of the variation in PC2 and PC3, respectively. No gene reaches adjusted P < 0.05 between the two control groups in the blocked comparison. Tests across the four inferred blocks among the eight controls identify 5,482 genes at this threshold. The study reports four samples for RNA extraction and sequencing and describes three biological replicates for its in-vitro assays. It does not identify a shared culture, preparation day or sequencing batch for a given suffix across groups.
3. **Model sensitivity.** Among genes matched to the published tables, blocking increases the number meeting both cutoffs by 1.80-fold for 50 nm and 1.48-fold for 500 nm. The table gives counts for 50 nm / 500 nm at 100 µg/mL against each series’ own control. Adjusted P values come from named DESeq2 coefficients at alpha = 0.05 before rounding. The effect-size cutoff uses apeglm estimates with absolute log2 fold change ≥ log2(1.5).

   | Gene set and model | Adjusted P < 0.05 | With effect-size cutoff |
   |---|---:|---:|
   | Matched genes, group only | 1,299 / 2,659 | 449 / 1,256 |
   | Matched genes, inferred block + group | 3,687 / 5,043 | 806 / 1,858 |
   | All filtered genes, group only | 1,553 / 3,028 | 570 / 1,495 |
   | All filtered genes, inferred block + group | 4,322 / 5,736 | 1,028 / 2,215 |

   As a negative control, I independently shuffled the four labels within each group and refitted DESeq2 39 times (seed 20261113). Using the same contrast extraction for all fits, observed counts were 4,320 / 5,736, versus shuffled medians of 1,506 / 2,967. Neither series reached its observed count in any shuffle. This tests whether the observed labels capture expression structure. It cannot identify a shared culture, experimental day or sequencing batch. The control uses DESeq2 contrast extraction throughout, whereas the table uses named coefficients. Counts near the threshold depend slightly on this choice.

4. **Shared response.** Both series show induction of inflammatory and chemokine transcripts, including *Ccl2*, *Ccl7* and *Cxcl10*, alongside *Gpnmb*, *Lgals3*, *Hspb1* and *Lcn2*. Synapse-related and glutamate-receptor gene sets shift downwards.
5. **Reactive-astrocyte signatures.** Pan-reactive and A2 gene sets show positive enrichment in both series. The A1 set does not reach significance in either series. These are transcriptomic signatures in bulk cultures.
6. **Differences at equal mass concentration.** The direct contrast compares (50 nm high dose − its control) with (500 nm high dose − its control) on the log2 scale. Step 08 reports 990 genes at adjusted P < 0.05. The full-data interaction supports a greater downward shift of homeostatic and ciliary gene sets in the 500 nm series. Complement enrichment favours the 50 nm series. Of the 990 full-model hits, 926 have the same direction in all four leave-one-block-out fits and all four individual blocks, and 144 also retain adjusted P < 0.05 in all four refits. The group-only model calls 89 interaction genes; all 89 agree with the blocked model. In a targeted three-set GSEA, cilia and complement retain direction and significance across all four refits. The homeostatic set retains direction but reaches adjusted P = 0.063 when block 2 is omitted. These within-study sensitivity checks use overlapping samples.
7. **Statistically stable expression changes with exploratory biological interpretation.** *Gpnmb* increases about tenfold in both series. In the 500 nm series, *Uchl1* and *Dnajc6* increase and *Pink1* decreases. These genes were selected for follow-up based on prior disease associations. Their expression changes in this culture cannot establish a Parkinson disease mechanism.
8. ***Npy* induction.** Estimated increases are about 30-fold at 50 nm and 120-fold at 500 nm. All four inferred blocks agree on direction. The published tables omit this gene, and independent confirmation is pending.

## Exploratory follow-up

| Candidate | Why it looked interesting | What the checks showed |
|---|---|---|
| *Snca* (α-synuclein) ↑ | padj 0.023 (50 nm) and 0.044 (500 nm) | The largest change occurs in one inferred block (rep3: +2.2 log2; others +0.3–0.4). Not significant when rep3 is left out (Fig. 6). |
| KEGG "Parkinson disease" pathway (NES +2.7) | Top-ranked pathway at 50 nm | 51% of its genes are shared mitochondrial / ribosome / proteasome genes; AD, HD, prion and ALS pathways rise with the same genes. After removing them, the residual signal is ubiquitin-system and stress genes (*Ubc*, *Uchl1*, *Park7*), not α-synuclein biology. Not significant at 500 nm. |
| "A1 neurotoxic astrocytes" | *C3* (an A1 marker) ↑ 15× at 50 nm | The A1 gene set does not reach significance. *C3* alone is highly variable between replicates (1.2–6.3 log2); its size × dose interaction is not significant. |
| Proliferation ↑ with 500 nm | Cell-cycle GO terms in the interaction contrast | The proliferation marker set does not reach significance in either series. |
| Sample swap (500C rep3 looks like dose 100) | PCA outlier | Decomposing each sample into the 50 nm and 500 nm response patterns (learned from the other 3 replicates) places it closer to 500 nm, and its 50D partner closer to 50 nm — the score favours the submitted label. The method passed 7 of 8 positive controls, so sample identity remains uncertain. |

## Methods

| Step | Tool (version) | Notes |
|---|---|---|
| Download | aria2 1.37, ENA | 64 FASTQ files, md5-verified |
| QC | FastQC 0.13.0, fastp 1.4.0, MultiQC 1.35 | Reads already trimmed by submitters (~2×73 bp); quality scores simplified (SRA Lite); no samples removed |
| Reference | Ensembl release 116, GRCr8 | cDNA + ncRNA transcripts, GTF of the same release |
| Quantification | Salmon 2.8.0, k = 31 | `--seqBias --gcBias`, library type detected as ISR (stranded), mapping 83.9–88.1% |
| Gene level | tximport 1.40.0 | 42,687 genes; versions stripped from transcript IDs |
| Filtering | — | ≥10 reads in ≥4 samples → 16,726 genes |
| Differential expression | DESeq2 1.52.0, apeglm 1.34.0 | `~ replicate + group`; DEG = padj < 0.05 and \|FC\| ≥ 1.5 (shrunken), the authors' criterion |
| Robustness | DESeq2 | four leave-one-block-out fits for within-series comparisons, within-block fold changes, minimal-filter refit |
| Enrichment | clusterProfiler 4.20.3, org.Rn.eg.db 3.23.0 | GO-BP ORA (universe = tested genes) and GSEA (Wald statistic); KEGG ORA/GSEA |
| Environment | R 4.6.1, Bioconductor 3.23, conda (conda-forge + bioconda) | exact versions in `results/logs/` |

**High-confidence DEGs** (`results/confidence/`): a DEG in the main fit **and** significant with the same sign in all four leave-one-replicate-out fits **and** changed in the same direction inside every inferred block. 50 nm: 610 of 1,028 (59%); 500 nm: 1,741 of 2,215 (79%). These percentages describe within-study stability for each high-dose versus control comparison.

## Agreement with the published analysis

- Genes matched to the published tables (by Ensembl ID, then by unique gene name): 12,502 of 17,825 (70.1%, 50 nm) and 12,510 of 18,037 (69.4%, 500 nm). The rest are mostly gene models that changed between Rnor 6.0 (Ensembl 90, used by the authors) and GRCr8 (~2,500), or genes our low-count filter removed (~2,800, but only 7 / 59 of the authors' DEGs).
- Fold-change correlation on shared genes: Spearman 0.93 / 0.94.
- Among matched genes, 474/497 (95.4%) and 1,122/1,156 (97.1%) published DEGs were recovered. All overlapping calls agreed in direction.
- Some disagreements involve low-count genes that the unpaired, shrunken model pulls to fold change ≈ 1 (reproduced in `scripts/12_authors_model.R`).

## Limitations

- *In vitro*, rat, neural-stem-cell-derived astrocyte cultures; bulk RNA-seq. Changes in cell composition (e.g. fewer ciliated cells) cannot be separated from changes inside cells.
- The submitted suffixes support an inferred block. PCA and within-group permutations assess expression structure. Identifying its biological or technical origin requires experimental records.
- No decoy-aware Salmon index (a GRCr8 decoy index needs ~20 GB RAM; not feasible on the 8 GB machine used). Quality scores in the public reads are simplified.
- Marker expression varies across inferred blocks and is compatible with differences in baseline composition. Cell proportions were not measured.
- The A1/A2 classification is a simplification (Escartin *et al.* 2021); signatures were defined in mouse, and 4 A1 genes have no same-named rat gene.
- Many analyses were exploratory; findings beyond the high-confidence set (and *Npy* in particular) should be treated as hypotheses.
- Methods §2.6 reports RNA extraction after seven days. Figure 2B shows astrocyte MTT viability at 100 µg/mL at approximately 70% of control for both sizes (visual estimates, three independent experiments). The authors describe the 500 nm result as a decreasing trend. This separate assay cannot quantify the contribution of viability or composition to any sequenced library.
- Equal nominal mass does not equalize particle number or total geometric surface area. For spheres of the same density at their nominal diameters, 50 nm gives 1,000 times the number and ten times the total area of 500 nm per unit mass. Delivered dose, aggregation and effective contact area were not measured by this calculation.

## Repository layout and how to run

```
scripts/        steps 00–21 via run_all.sh, from raw reads through dose analyses, figures and block validation
results/        logs, tables and figures (large intermediate files are not tracked)
data/           raw reads and reference (not tracked; recreated by steps 01–03)
samples.csv     sample sheet: run, GSM, particle series, dose, group, replicate
```

```bash
conda activate rnaseq                 # conda-forge + bioconda env (results/logs/conda_environment.yml)
bash scripts/run_all.sh 01 21         # from a project clone with dependencies installed
bash scripts/run_all.sh 06 16         # R analysis only, from Salmon output
bash scripts/run_all.sh 17 21         # requires txi.rds, dds.rds and upstream annotation
```

For a processed-data restart, step 06 needs the GTF and per-run Salmon outputs. Steps 07 onward use local RDS objects created by the pipeline. A fresh clone includes compressed result tables, while several downstream scripts expect plain TSVs. Create these inputs before restarting downstream. Step 00 currently assumes a fixed home-directory project path and can move files from Downloads. The Conda export covers command-line tools only in part; install aria2 and the R/Bioconductor packages listed in the session log separately.

The [sample map](results/block_validation/biosample_map.tsv), [model counts](results/block_validation/model_counts.tsv), and [within-group permutation results](results/block_validation/permutation_summary.tsv) document the inferred-block analysis. To regenerate model counts or permutation draws, pass the project root and an existing output directory to `scripts/compare_counts.R` or `scripts/block_control.R`.

## References

- Marcellus KA, Bugiel S, Nunnikhoven A, Curran I, Gill SS. Polystyrene nano- and microplastic particles induce an inflammatory gene expression profile in rat neural stem cell-derived astrocytes in vitro. *Nanomaterials* 2024;14(5):429. doi:10.3390/nano14050429
- Liddelow SA *et al.* Neurotoxic reactive astrocytes are induced by activated microglia. *Nature* 2017;541:481–487.
- Escartin C *et al.* Reactive astrocyte nomenclature, definitions, and future directions. *Nat Neurosci* 2021;24:312–325.
- Patro R *et al.* Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods* 2017;14:417–419.
- Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* 2015;4:1521.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 2014;15:550.
- Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 2019;35:2084–2092.
- Wu T *et al.* clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2021;2:100141.

## Interaction sensitivity outputs

Under `results/interaction/`, `interaction_genes.tsv.gz` contains the full and group-only estimates, four omission fits, individual-block differences in log2 normalized counts (pseudocount 0.5), and pointwise 95% Wald intervals. Positive values mean that the signed response is greater in the 50 nm series at equal nominal mass concentration. The baseMean cutoffs 10, 20 and 50 are reporting filters applied after genome-wide adjustment. Targeted GSEA ranks finite Wald statistics and adjusts across three selected sets within each fit. This is a follow-up to the original discovery analysis. [Selected gene estimates](results/interaction/selected_interactions.pdf)

## Dose-response analysis and manuscript

At 100 µg/mL, the ciliary and complement signatures reached the adjusted significance threshold in the inferred-block model, the group-only model and all four block-omission fits. At 10 µg/mL the ciliary signature retained its direction in all six fits and reached the threshold in five. Omitting block 3 raised its adjusted P value to 0.161. [Sensitivity results](results/dose_sensitivity/pathway_stability.tsv) · [Manuscript Figure 3](results/figures/fig7_dose_response.pdf)

The complete article is available as [editable Word](manuscript/Manuscript.docx), [PDF](manuscript/Manuscript.pdf) and [Markdown source](manuscript/Manuscript.md). Figure 1 uses the original volcano plots, Figure 2 the original signature plots, and Figure 3 the across-dose sensitivity plot. All three source figures are in [results/figures](results/figures). The [journal figure archive](manuscript/Figures.zip) uses manuscript numbering.

Step 21 rebuilds the BioSample map from the archived NCBI XML, compares model counts, runs the within-group permutation control and writes its summary under [results/block_validation](results/block_validation). The [verification script](scripts/verify_results.py) checks the manuscript counts directly against the tracked tables. Run it from the root with `python3 scripts/verify_results.py`. To rebuild Word, PDF and the figure archive from the Markdown and canonical figures, use `bash scripts/build_manuscript.sh` (Pandoc, python-docx and LibreOffice required).
