## Results

### High-dose responses were widespread, while the number of genes distinguishing the two particle series depended on the model

We reanalysed 32 RNA-seq libraries from rat neural stem cell-derived astrocyte cultures exposed to polystyrene for seven days (GSE256038). Each particle series comprised four samples at each of four nominal mass concentrations (0, 1, 10 and 100 µg/mL). The four numbered BioSample suffixes occurred once in each group and matched the within-group GEO accession ranks in all 32 libraries. We used these suffixes as inferred blocks. The records examined did not establish whether matching suffixes denoted a common culture preparation, experimental day or processing batch. [Source study](https://doi.org/10.3390/nano14050429)

At 100 µg/mL, 1,028 genes in the 50 nm series and 2,215 in the 500 nm series met both the adjusted P < 0.05 and 1.5-fold change criteria relative to their respective controls (apeglm estimates). Among published differential-expression calls that could be matched to genes in our analysis, 474 of 497 (95.4%) and 1,122 of 1,156 (97.1%) were recovered, respectively. All recovered calls agreed in direction.

To test whether the responses differed at equal nominal mass concentration, we subtracted the log2 dose-versus-control response in the 500 nm series from the corresponding response in the 50 nm series. The full inferred-block model identified 990 genes at 100 µg/mL (Benjamini–Hochberg adjusted P < 0.05). Of these, 926 retained the direction of the full-model estimate in all four leave-one-block-out fits and in all four individual block contrasts. A subset of 144 also retained significance in all four omission fits. The group-only model identified 89 genes, all of which were significant and concordant in the full model. These fits use the same study and quantify within-study sensitivity.

### Ciliary and complement signatures differed between particle series at 100 µg/mL

We followed up two gene sets identified in the high-dose analysis by ranking genes on the direct response contrast at each dose. A positive normalized enrichment score (NES) indicates enrichment toward a greater signed response in the 50 nm series. At 100 µg/mL, ciliary/ependymal transcripts showed a greater relative decrease in the 500 nm series (full-model NES 1.98, adjusted P = 0.00030). Complement activation genes were enriched toward the 50 nm response (NES 1.99, adjusted P = 0.0016). Both enrichments retained direction and significance in the group-only fit and in all four omission fits (Figure 1B). In individual blocks, the descriptive mean contrast for the seven ciliary genes ranged from −0.03 to 1.61 log2 units. For the 18 complement genes it ranged from 0.40 to 0.79 (Figure 1A). Omitting one member gene at a time preserved positive enrichment in all 25 high-dose tests (largest BH-adjusted P = 0.0027).

At 10 µg/mL, the ciliary set showed positive enrichment in all six fits (Table 1). Adjusted P was 0.0056 in the full model and 0.035 in the group-only model. Omitting block 3 attenuated the NES from 1.82 to 1.46 and increased adjusted P to 0.161. Five genes contributed to the full-model leading edge (*Rsph4a*, *Dynlrb2*, *Dnah5*, *Foxj1* and *Dnah9*). The complement set did not reach the adjusted significance threshold at this dose in any fit. At 1 µg/mL neither set reached the threshold. These dose-specific tests describe the observed evidence at the sampled concentrations. They do not estimate a biological onset threshold or directly test a difference between the dose thresholds of the two signatures.

**Table 1. Sensitivity of enrichment to the inferred-block model.** Counts give fits with BH-adjusted P < 0.05 out of six overlapping fits. Adjustment covers the two selected gene sets at three doses within each fit.

| Nominal concentration (µg/mL) | Ciliary/ependymal | Complement activation |
|---:|---:|---:|
| 1 | 0/6 | 0/6 |
| 10 | 5/6 | 0/6 |
| 100 | 6/6 | 6/6 |

**Figure 1. Dose-specific differences in transcriptional response between the 50 nm and 500 nm series.** A. Each point gives the mean gene-level difference of dose-versus-control responses in one inferred block, calculated from log2(normalized count + 0.5). Positive values indicate a larger signed change at 50 nm. Lines connect concentrations sharing a sample suffix. B. Colour encodes the GSEA normalized enrichment score from genes ranked by the DESeq2 Wald statistic. Dots mark adjusted P < 0.05. The full model includes group and inferred block, the group-only model uses all 32 libraries, and each omission fit uses 24. The targeted analyses use signatures selected in this dataset. Bulk RNA-seq cannot distinguish altered expression within a cell type from a change in the relative abundance of cells expressing these markers.
