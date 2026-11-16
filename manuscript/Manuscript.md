---
title: "Complement transcript responses to two polystyrene particle sizes in rat astrocyte cultures after accounting for response amplitude"
---

Mohammad Javad Eslami

Independent researcher, Tehran, Iran

Correspondence: Mohammad Javad Eslami, mohamadjavadeslami10@gmail.com

## Abstract

**Background** Repeated exposure to polystyrene particles alters inflammatory gene expression in rat astrocyte cultures. Whether the magnitude or the pattern of the response differs between particle sizes requires a direct comparison at the same mass concentration, and requires separating a difference in overall response amplitude from a difference in which genes respond.

**Results** We compared seven-day responses to 50 nm and 500 nm polystyrene in 32 libraries at equal nominal mass concentrations. At 100 µg/mL the contrast between the two responses identified 990 genes in a model with inferred sample blocks. Of these, 144 retained significance in all four block-omission fits and agreed in direction in every block; 50 of the 144 were also significant without blocking. To assess whether the three selected gene sets differed more than expected from the larger overall response in the 500 nm series, we compared them with genes matched on their signed 500 nm response. Complement transcripts exceeded this background by 0.457 log₂ units (one-sided empirical P = 0.000050; adjusted P = 0.00015). The estimated excess was 0.006 for ciliary/ependymal and 0.106 for homeostatic markers (adjusted P = 0.480 and 0.377). Ciliary enrichment at 10 µg/mL weakened when block 3 was omitted. At 100 µg/mL, higher nestin and lower aquaporin-4 expression in the 500 nm cultures were compatible with differences in cell state or composition.

**Conclusions** Of three exploratory sets, complement showed the clearest excess response difference after matching on the estimated 500 nm effect. The greater ciliary and homeostatic decline in the 500 nm series was compatible with the genome-wide response pattern. Bulk RNA measurements do not establish complement activity or identify the cellular source of the marker shifts.

**Keywords** Polystyrene, astrocyte cultures, RNA sequencing, nanoplastics, microplastics, cilia, complement, transcriptomics

## Background

Marcellus and colleagues found that repeated exposure to polystyrene particles altered inflammatory and metabolic gene expression in rat neural stem cell-derived astrocyte cultures [1]. At 100 µg/mL, they reported 531 differentially expressed genes for the 50 nm series and 1,274 for the 500 nm series. These counts describe each exposure relative to its own control. They leave a separate question open: for which genes or gene sets does the change from control differ between the two series? The public sequencing libraries and published gene tables allow that contrast to be estimated directly [2].

Estimating it raises a question of interpretation that is easy to overlook. If one preparation simply produces a larger response than the other at the same nominal mass, then every gene that responds in both series will show a non-zero between-series difference, and any coherently regulated gene set will appear to differ. Distinguishing a difference in amplitude from a difference in which genes respond therefore requires the overall relationship between the two responses to be described before set-level results are interpreted.

A distinction between shared and divergent responses is biologically relevant in this culture system. Inflammatory and stress responses have also been observed in primary human astrocytes exposed to polystyrene under different conditions [3]. The interpretation of ciliary transcripts is less straightforward. In the postnatal mouse brain, *Foxj1* contributes to ependymal differentiation and marks a subset of astrocytes [4]. In a bulk culture, a decrease in *Foxj1* and associated transcripts can therefore reflect a change in marker expression, a change in the representation of marker-positive cells, or both.

We asked whether equal nominal mass concentrations of the two particle sizes produced different transcriptional responses after seven days. A recurring suffix in the sample records suggested a possible block shared across exposure groups. We examined how the choice of model affected gene-level discoveries, described the overall relationship between the two responses, and examined all three sets at high dose, then assessed ciliary and complement signals at the two lower doses. The gene sets were selected after examining the high-dose data, so the dose-specific analyses are exploratory.

## Methods

### Source material and sample structure

GSE256038 comprises paired-end libraries from rat fetal neural stem cell-derived astrocyte cultures exposed to polystyrene for seven days [1,2]. Each series (50 nm and 500 nm nominal diameter) comprised a control and exposures at 1, 10 and 100 µg/mL, with four libraries per group. All 32 libraries were retained. The sample size was fixed by the source experiment. The original study reported four samples for RNA extraction and three independent experiments for the separate viability assay [1].

Each GEO sample was linked to its SRA run and BioSample accession. Within each of the eight groups, BioSample titles carried the suffixes 01–04, matching the order of the GEO accessions. These suffixes defined four inferred blocks. The available metadata did not establish whether a suffix denoted a shared culture, preparation day or technical batch. Principal component analysis used the 500 most variable genes after variance stabilization. Separate regressions of component scores on suffix assessed their association.

### Quantification and comparison with the published analysis

The two FASTQ files per library were verified against ENA MD5 checksums. FastQC and fastp provided per-file and per-library quality reports, which were aggregated with MultiQC [5,6]. Salmon 2.8.0 quantified the paired reads with automatic library-type detection and sequence- and GC-bias correction [7]. The index contained Ensembl release 116 rat cDNA and noncoding RNA sequences (GRCr8, k = 31), without genomic decoys. Transcript identifiers were matched to the same release's gene annotation after removal of version suffixes. Gene-level counts and sample-specific transcript-length offsets were imported with tximport 1.40.0 [8]. Genes with at least ten counts in at least four libraries entered the DESeq2 1.52.0 analysis [9].

The blocked negative-binomial model included an eight-level exposure group and a four-level block as additive fixed effects. A group-only model used the same libraries and filtered gene universe. Within each series, high-dose discoveries required a Benjamini–Hochberg (BH) adjusted Wald P value below 0.05 and an absolute apeglm-shrunken log₂ fold change of at least log₂(1.5) [10,11]. Published genes were matched by Ensembl identifier, with unique gene symbols as a fallback when identifiers were absent from the updated annotation. Ambiguous matches and symbol matches conflicting with a direct identifier match were excluded. Recovery was calculated among published discoveries represented in the comparison. For the within-series stability display in Fig. 1, highlighted genes also retained significance and direction in all four block-omission fits and agreed in direction within every observed block.

We shuffled the four suffix labels independently within each group and refitted DESeq2 for 39 draws (seed 20261113). Each fit retained the same group composition and number of block parameters. Observed and shuffled fits used identical contrast extraction and a BH threshold of 0.05 alone. This comparison assessed whether the observed suffix arrangement captured expression structure in these samples.

### Relationship between the two within-series responses

At 100 µg/mL we summarized each series by the median absolute shrunken log₂ fold change and regressed 50 nm on 500 nm estimates through the origin. For an exploratory conditional test we extracted unshrunken coefficients from the fitted DESeq2 model and ordered genes by their signed estimated 500 nm effect. For each represented set member we drew a gene uniformly from its rank neighbourhood (±50 ranks, clipped at the boundaries, including the member). Draws were independent across members, allowing repeats. In 20,000 draws per set (seed 20261115) we compared the observed unweighted mean interaction with the matched-set means. One-sided P values used (1 + number of draws at least as large as observed)/20,001, followed by BH adjustment across three sets. Matching conditions on an estimated component of the contrast, does not propagate its uncertainty and does not preserve between-gene correlation. This sensitivity analysis uses model-derived effects; it is not a new fitted interaction model.

### Differences between exposure responses

At each nonzero dose, the gene-level contrast was (50 nm exposure − 50 nm control) − (500 nm exposure − 500 nm control) on the log₂ scale. Positive estimates indicate that the signed transcriptional response was greater in the 50 nm series. Because this contrast is a difference of signed responses, a positive estimate can arise from induction in the 50 nm series, from a decline in the 500 nm series, or from both; within-series estimates are therefore reported alongside it. Direct contrasts used unshrunken coefficients, DESeq2 independent filtering and BH adjustment within each contrast. Pointwise 95% Wald intervals describe uncertainty in individual coefficients. A three-degree-of-freedom likelihood-ratio test compared a block + series × categorical dose model with a reduced model containing block, series and dose main effects.

Each block-omission fit retained 24 libraries and the original gene universe. Normalization and dispersions were re-estimated. Descriptive within-block contrasts were calculated from log₂(normalized count + 0.5). We recorded sign agreement separately from significance and compared the blocked-model discoveries with those from the group-only model. Sensitivity to low abundance was assessed by restricting reported results to genes with a full-model mean normalized count of at least 10, 20 or 50 after multiple-testing adjustment.

### Gene sets and sensitivity analyses

Exploratory analyses examined Gene Ontology biological processes, KEGG pathways and a panel of eight astrocyte and cell-type marker sets. The latter included pan-reactive and historical A1 and A2 marker lists derived from mouse experiments [12]. Within each high-dose contrast, the panel analysis ranked mapped genes with evaluable adjusted P values by a Wald statistic reconstructed from the exported effect estimate and standard error.

Three sets were then selected for targeted follow-up after inspection of the high-dose results: a ciliary/ependymal marker list, a homeostatic astrocyte marker list and the complement activation term (GO:0006956). The author constructed the homeostatic and ciliary lists after examining the high-dose data, including the Aqp4 change and GO cilium-movement/axoneme signal. Some ciliary members (including Foxj1 and Dnah9) contributed to the original signal. The proliferation, oligodendrocyte and neuronal panels were also hand-selected. For PAN/A1/A2, the hypothesis preceded inspection but the author reconstructed the gene lists from memory afterwards. The resulting membership and selection timing are documented in sets_source.tsv. Complement members were obtained from the rat annotation package org.Rn.eg.db, including genes annotated to descendant terms. The high-dose analysis across model fits tested all three sets, with BH adjustment covering the three sets within each fit. Overlapping high-dose tests used the same rank-based NES and raw P from the three-set analysis before separate family adjustment. The subsequent across-dose analysis carried forward the ciliary and complement sets only, with BH adjustment covering those two sets at all three nonzero doses within each fit.

For the targeted analysis, preranked gene-set enrichment analysis (GSEA) used clusterProfiler and full-precision, finite DESeq2 Wald statistics [13,14]. Where multiple gene identifiers shared a symbol, the entry with the largest absolute statistic was used, with identifier order breaking ties. We tested sets of 4–500 represented genes. Seven ciliary, ten homeostatic and 18 complement members were represented. A positive normalized enrichment score (NES) denotes enrichment toward the 50 nm response. We examined the leading edge (members contributing to the enrichment peak), and repeated high-dose GSEA after removing each represented gene from its set. The 25 resulting P values were adjusted together; because these tests share a single ranking and differ by one member, the adjustment is a stability display rather than a control of error across independent tests. Block-level scores in Fig. 3 are unweighted means of member-gene log₂ response differences. All omission analyses used overlapping samples or sets.

### Sensitivity to externally defined marker lists

Before rerunning enrichment, we timestamped the marker-selection rules and locked memberships in version control. We transcribed PAN/A1/A2 from the labelled panels of Liddelow et al. (Fig. 1a and Extended Data Fig. 12a) and checked each label against rat and mouse primer tables in their supplement [12]. The published panels use Ugt1a, whereas the original analysis used Ugt1a1. The panels label Hspb1, while the primer tables print Hsbp1; the panel label was used and primer specificity was not inferred. We retained the submitted punctuation normalization for the two H2 labels. The marker audit records unmatched genes. For a separate sensitivity analysis, all mouse-supported astrocyte and ependymal markers were taken from the 27 March 2020 PanglaoDB release [19]. We used case-insensitive symbol matching to the rat ranking without imputing missing orthologues. These are broader cell-identity proxies for the original homeostatic and motile-cilia sets. Full blocked-model GSEA at high dose used the same ranking and tie-break rules and BH across five external lists per contrast. A further exploratory test applied signed-effect matching to the two external identity sets, adjusting over two tests. The same samples were used throughout.

### Computational assistance

OpenAI Codex assisted with code execution, verification and manuscript preparation. The supplied external review used Claude to audit archived tables and propose analyses and textual revisions. The revised computations were run locally using the original tximport data and the archived results were checked separately. Author review and approval are required before submission.

## Results

### Reprocessing recovered high-dose responses in both series

All 32 libraries had an inferred inward, stranded reverse orientation; mapping rates ranged from 83.9% to 88.1%. Transcript aggregation yielded estimates for 42,687 genes, of which 16,726 passed the count filter. The sample suffix was associated with the second and third principal components (R² = 0.78 and 0.91, respectively). At 100 µg/mL, the observed block labels yielded 4,320 and 5,736 discoveries at adjusted P < 0.05 for the 50 nm and 500 nm comparisons with their controls. Across the 39 shuffles, the corresponding counts ranged from 1,382 to 1,662 and from 2,839 to 3,186. This association motivated reporting the blocked and group-only models together.

Adding the 1.5-fold effect-size criterion gave 1,028 high-dose discoveries in the 50 nm series and 2,215 in the 500 nm series (Fig. 1). Of these, 610 and 1,741, respectively, met the within-series stability criteria. Among the 531 and 1,274 discoveries reported in the original study, 497 and 1,156 could be compared after gene matching and filtering. The reanalysis recovered 474 of 497 (95.4%) and 1,122 of 1,156 (97.1%). Every recovered gene agreed in direction. Within these two high-dose contrasts, Spearman correlations between published and reanalysed log₂ fold changes across all matched genes were 0.93 and 0.94, respectively.

### The two series differed in overall response amplitude

The median absolute shrunken log₂ fold change was 0.113 for 50 nm and 0.162 for 500 nm. The through-origin slope of 50 nm on 500 nm estimates was 0.483 across 16,726 genes, 0.484 among discoveries in the 500 nm series, and 0.638 among discoveries in both series. For unshrunken estimates the all-gene slope was 0.453; the corresponding medians were 0.207 and 0.223. These asymmetric slopes describe an association between two estimated quantities and do not imply a uniform twofold response for individual genes.

At high dose the mean unshrunken between-series contrast for 18 complement genes was +0.589 log₂ units, compared with +0.132 for sets matched on their 500 nm responses (excess +0.457; BH-adjusted P = 0.00015). The corresponding observed and matched means were +0.745 and +0.739 for seven ciliary markers (excess +0.006; adjusted P = 0.480) and +0.470 and +0.364 for ten homeostatic markers (excess +0.106; adjusted P = 0.377). The positive raw contrasts for the latter two sets were therefore also common in their matched backgrounds.

### A subset of high-dose differences persisted across model specifications

The joint interaction test identified 852 genes whose response differed between the two series across doses. In dose-specific tests, one gene met the genome-wide adjusted threshold at 1 µg/mL, none at 10 µg/mL and 990 at 100 µg/mL, among the 13,807 genes with evaluable adjusted P values in the high-dose contrast. Of the 990 high-dose discoveries, 707 had positive and 283 negative contrast estimates. All 990 retained their estimated direction in each leave-one-block-out fit. For 926, the observed contrast in every individual block also agreed in direction. Of these, 144 retained adjusted P < 0.05 in each omission fit. Restricting the report to mean normalized counts of at least 20 retained 987 of the 990 discoveries and all 144 of the more stringent set.

Removing the block term reduced the high-dose list to 89 genes, all significant and concordant in the full model. Fifty of these also belonged to the 144-gene stringent set. For *Foxj1*, the difference in response was 0.81 log₂ units (pointwise 95% interval 0.51–1.12); for *Gja1*, it was 1.15 (0.74–1.57). Both met the stringent block-omission criteria, but they behaved differently without blocking: *Gja1* remained significant in the group-only model (adjusted P = 0.022) whereas *Foxj1* did not (adjusted P = 0.76). *C3* had a larger point estimate (2.65 log₂ units) but greater uncertainty and a genome-wide adjusted P value of 0.090.

### Ciliary, homeostatic and complement signatures at the highest concentration

Both series showed positive enrichment of pan-reactive and A2 marker sets after high-dose exposure (Fig. 2). Ciliary/ependymal and homeostatic astrocyte markers shifted downward more strongly in the 500 nm series. The exploratory direct contrast also identified complement activation toward the 50 nm response (GO enrichment NES 1.88, adjusted P = 0.0084).

In the targeted three-set analysis at 100 µg/mL, all three sets were enriched toward the 50 nm response in the full model (Table 1). Ciliary and complement sets met the adjusted threshold in all six fits; the homeostatic set met it in five, reaching adjusted P = 0.063 when block 2 was omitted. After individual genes were removed from the ciliary and complement sets, the largest adjusted P value across all 25 high-dose omissions was 0.0027.

The within-series estimates distinguish these three results. For the ciliary set, five of seven represented members declined significantly in the 500 nm series (*Dnah9* −2.24, *Dnah5* −1.45, *Rsph4a* −1.38, *Dynlrb2* −1.25 and *Foxj1* −0.79 log₂ units), whereas only two declined significantly in the 50 nm series and by smaller amounts (*Dnah5* −0.80, *Dnah9* −0.71). The homeostatic set behaved similarly: six of ten members declined significantly in the 500 nm series (including *Atp1a2* −1.56, *Aqp4* −1.43, *Glul* −1.25 and *Gja1* −1.24) and none did so in the 50 nm series. The complement set differed in kind. Six of 18 members rose significantly in the 50 nm series, including *C3* (+3.94), *A2m* (+0.86) and *Cfb* (+0.66), while nine declined significantly in the 500 nm series, including *Clu* (−1.17) and *C1r* (−1.01). The complement contrast therefore combined induction in one series with repression in the other, whereas the ciliary and homeostatic contrasts were driven principally by the larger declines in the 500 nm series.

### External marker lists

Of the 63 mouse-supported astrocyte and 59 ependymal markers in PanglaoDB, 58 and 44 entered the rat gene ranking. Both sets were depleted in the 500 nm high-dose series and enriched toward the 50 nm series on the direct contrast (NES 2.51 and 2.26; adjusted P = 1.0 × 10⁻⁹ and 9.6 × 10⁻⁷). In a secondary matched-background comparison, the broader astrocyte set retained a positive excess (+0.179 log₂, adjusted P = 0.021), while the ependymal set gave +0.144 (adjusted P = 0.067). Thus the ependymal direction generalized to an external list; evidence beyond matched background depended on the identity programme tested. The source-checked Liddelow panel retained positive PAN and A2 enrichment in each within-series contrast and no significant A1 enrichment. None of the three differed significantly on the direct interaction. Gene coverage was 11/13, 7/12 and 9/12 for PAN, A1 and A2.

### Dose dependence and its sensitivity to one block

The average high-dose complement contrast was positive in each individual block (0.40–0.79 log₂ units). Ciliary scores showed more variation (−0.03 to 1.61), with the largest contrast in block 3 (Fig. 3A). At 10 µg/mL, ciliary enrichment was positive in all six fits. Omitting block 3 reduced the NES from 1.82 to 1.46 and raised adjusted P from 0.0056 to 0.161. The full-model leading edge comprised *Rsph4a*, *Dynlrb2*, *Dnah5*, *Foxj1* and *Dnah9*. Complement enrichment at 10 µg/mL, and both sets at 1 µg/mL, did not reach significance in any fit (Table 2).

### Cell-type marker expression differed between the two high-dose groups

Group means of cell-type markers provide descriptive evidence relevant to culture state. Relative to its own control, the 500 nm high-dose group showed a 3.3-fold increase in the progenitor marker *Nes* (4,072 to 13,240 normalized counts) and a decrease in the mature astrocyte marker *Aqp4* (35,518 to 13,003). The corresponding changes in the 50 nm high-dose group were smaller (*Nes* 3,335 to 5,848; *Aqp4* 35,650 to 28,951), and *Gfap* rose more in the 50 nm series (+1.04 log₂ units) than in the 500 nm series (+0.18). These are descriptive group means from bulk libraries and were not tested formally (Additional file 1).

**Table 1. Targeted enrichment of the direct 50 nm minus 500 nm response contrast at 100 µg/mL.**

| Gene set | Represented members | NES | Adjusted P | Fits meeting criterion |
|---|---:|---:|---:|---:|
| Ciliary/ependymal | 7 | 1.98 | 0.00015 | 6/6 |
| Complement | 18 | 1.97 | 0.00058 | 6/6 |
| Homeostatic astrocyte | 10 | 1.92 | 0.00058 | 5/6 |

NES and adjusted P are from the full blocked model. The final column counts fits reaching adjusted P < 0.05 out of six overlapping fits. BH adjustment covered three targeted sets per fit.

**Table 2. Enrichment of the direct 50 nm minus 500 nm response contrast across doses.**

| Dose (µg/mL) | Gene set | NES | Adjusted P | Fits meeting criterion |
|---:|---|---:|---:|---:|
| 1 | Ciliary/ependymal | 0.51 | 0.974 | 0/6 |
| 1 | Complement | 1.09 | 0.377 | 0/6 |
| 10 | Ciliary/ependymal | 1.82 | 0.0056 | 5/6 |
| 10 | Complement | 1.33 | 0.202 | 0/6 |
| 100 | Ciliary/ependymal | 1.98 | 0.0003 | 6/6 |
| 100 | Complement | 1.97 | 0.0016 | 6/6 |

NES and adjusted P are from the full blocked model. The final column counts fits reaching adjusted P < 0.05 out of six overlapping fits. BH adjustment covered six targeted tests per fit. This analysis carried forward two of the three sets in Table 1.

## Discussion

This reanalysis separates two features of the original experiment: a broad response shared by both particle series and differences in that response at the same nominal mass concentration. The shared response is supported by close agreement with the published high-dose fold changes despite different quantification pipelines [1].

The high-dose response depended on the metric used to describe its amplitude. The 500 nm estimates were more dispersed after shrinkage, but the slope also depended on which series was used as predictor and which genes were included. Matching directly on the signed 500 nm coefficient avoided assuming that one global slope described every gene. This changed the interpretation of the ciliary and homeostatic findings: their larger declines remained observable, but sets with comparable 500 nm effects yielded similar contrasts.

Complement retained an excess over the matched background. Its members were heterogeneous: for example, C3 increased strongly in the 50 nm series and had a smaller positive estimate in the 500 nm series, while Clu declined strongly in the latter. These patterns support divergence at the level of the set without implying that individual genes generally reversed direction. The one-sided test used sets selected from this dataset; the small empirical P value does not establish independent replication or complement protein activity. Coregulation among set members remains a limitation of gene-based resampling [17].

The cell-type marker data bear on the ciliary and homeostatic signals in particular. The 500 nm high-dose cultures showed a marked rise in *Nes* and fall in *Aqp4* relative to control, changes that were much smaller in the 50 nm high-dose cultures. A shift in the differentiation state or composition of the culture would reduce transcripts associated with mature astrocyte and ependymal identity together, which is the pattern observed for *Foxj1*, *Aqp4*, *Gja1*, *Glul* and *Atp1a2*. *Foxj1* is relevant to both readings because of its role in ependymal differentiation and its expression in a subset of astrocytes [4]. Bulk sequencing cannot locate a response within a particular cell population, and cell proportions were not measured; these marker comparisons are descriptive and cannot quantify the contribution of composition. They do, however, indicate that a compositional explanation is not merely a theoretical alternative in this dataset. The pan-reactive and A2 scores in Fig. 2 likewise describe similarities to published expression signatures; assigning a protective or harmful phenotype would require functional measurements [12,16].

Within the concentrations sampled, the intermediate-dose evidence was stronger for ciliary genes, and its dependence on block 3 matters. The other three blocks together retained a positive enrichment score, yet the statistical support weakened. Two further considerations limit how much weight this observation can carry. First, the three doses were analysed as separate contrasts; the study therefore tests whether an interaction is detectable at a sampled concentration, and does not establish a monotonic dose trend, a difference in onset between pathways or a biological threshold. Second, block 3 at 10 µg/mL contains the one library whose identity our own analysis leaves unresolved: sample 500C-03 has a first principal component score of 22.1, within the range of the 500 nm high-dose samples (18.6 to 35.5) and far from its three group mates. A decomposition of each sample into the two series' response patterns favoured the submitted label, but that method passed only seven of eight positive controls, so the label cannot be confirmed. The intermediate-dose result should be read as an observation that depends on a single block, and the submitted label was retained; removing the entire block does not isolate the influence of that one library.

The ciliary list was selected partly to explain a GO signal already observed in these same samples; overlapping genes make part of that comparison circular. External identity signatures provide an independent definition, while reusing the same cultures. The two original marker lists were specified in code, and all three sets were chosen after inspection of the high-dose results. Their adjusted P values quantify evidence within the follow-up test family; they do not correct for exploratory selection, nor for the choice of which genes to include in a hand-assembled list.

Exposure conditions further constrain interpretation. For ideal spheres of equal density, reducing the diameter tenfold at fixed mass increases particle number a thousandfold and total geometric surface area tenfold. Sedimentation and diffusion can also affect the dose reaching cells [18]. The experiment therefore compares two exposure preparations at a specified mass concentration, and the amplitude difference reported above may reflect delivered dose rather than any property of particle size as such. The separate seven-day viability assay showed an MTT signal of a reduction relative to controls at 100 µg/mL in the original Fig. 2B; the authors described the 500 nm result as a decreasing trend [1]. These measurements cannot be assigned to individual sequencing libraries. In primary human astrocytes, shorter exposure to different polystyrene preparations also elicited inflammatory and oxidative stress responses [3]. That study provides biological context for the shared inflammatory response but did not analyse ciliary or complement transcripts.

## Conclusions

Among three exploratory gene sets, complement showed a high-dose response difference exceeding that observed for genes matched on their estimated 500 nm response. Ciliary and homeostatic markers declined more in the 500 nm series, but their mean contrasts were consistent with their matched backgrounds. The intermediate-dose ciliary signal was sensitive to omission of the block containing 500C-03. The marker pattern warrants attention to culture state and composition in interpreting these bulk data.

## Abbreviations

BH, Benjamini–Hochberg; GEO, Gene Expression Omnibus; GSEA, gene set enrichment analysis; NES, normalized enrichment score; SRA, Sequence Read Archive; MTT, methylthiazolyldiphenyl-tetrazolium bromide.

## Declarations

**Ethics approval and consent to participate.** This secondary analysis used public rat sequencing data. No new animals, specimens or human participants were involved. The experimental procedures for the source study are reported by Marcellus et al. [1].

**Consent for publication.** Not applicable.

**Availability of data and materials.** Raw reads and original differential-expression tables are available under GSE256038 [2]. The preceding analysis (v1.2.0) is archived in Zenodo [15]. Supplementary revision files contain the updated analysis code, results, marker table and processed R objects used for the present revision.

**Competing interests.** The author declares no competing interests.

**Funding.** The author funded the study with personal resources.

**Authors' contributions.** MJE designed the reanalysis, analysed and interpreted the data, and prepared the manuscript.

**Acknowledgements.** This work used the sequencing libraries and gene-level results deposited by the original investigators [1,2].

## References

1. Marcellus KA, Bugiel S, Nunnikhoven A, Curran I, Gill SS. Polystyrene Nano- and Microplastic Particles Induce an Inflammatory Gene Expression Profile in Rat Neural Stem Cell-Derived Astrocytes In Vitro. Nanomaterials. 2024;14:429. https://doi.org/10.3390/nano14050429

2. Marcellus KA, Bugiel S, Nunnikhoven A, Curran IH, Gill S. Polystyrene nano- and microplastic particles induce an inflammatory gene expression profile in neural stem cell-derived astrocytes. Gene Expression Omnibus. 2024. GSE256038. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256038

3. Kiran S, Thilanka U, Xue Y, Sang QXA. Polystyrene Micro- and Nanoplastic Exposure Triggers an Activation and Stress Response in Human Astrocytes. International Journal of Molecular Sciences. 2025;26:11273. https://doi.org/10.3390/ijms262311273

4. Jacquet BV, Salinas-Mondragon R, Liang H, Therit B, Buie JD, Dykstra M, et al. FoxJ1-dependent gene expression is required for differentiation of radial glia into ependymal cells and a subset of astrocytes in the postnatal brain. Development. 2009;136:4021-4031. https://doi.org/10.1242/dev.041129

5. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34:i884-i890. https://doi.org/10.1093/bioinformatics/bty560

6. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32:3047-3048. https://doi.org/10.1093/bioinformatics/btw354

7. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods. 2017;14:417-419. https://doi.org/10.1038/nmeth.4197

8. Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research. 2016;4:1521. https://doi.org/10.12688/f1000research.7563.2

9. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8

10. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society Series B: Statistical Methodology. 1995;57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

11. Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 2019;35:2084-2092. https://doi.org/10.1093/bioinformatics/bty895

12. Liddelow SA, Guttenplan KA, Clarke LE, Bennett FC, Bohlen CJ, Schirmer L, et al. Neurotoxic reactive astrocytes are induced by activated microglia. Nature. 2017;541:481-487. https://doi.org/10.1038/nature21029

13. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences. 2005;102:15545-15550. https://doi.org/10.1073/pnas.0506580102

14. Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021;2:100141. https://doi.org/10.1016/j.xinn.2021.100141

15. Eslami MJ. RNA-seq reanalysis of ciliary and complement responses to polystyrene exposure in rat astrocyte cultures. Version 1.2.0. Zenodo. 2026. https://doi.org/10.5281/zenodo.24197910

16. Escartin C, Galea E, Lakatos A, O'Callaghan JP, Petzold GC, Serrano-Pozo A, et al. Reactive astrocyte nomenclature, definitions, and future directions. Nature Neuroscience. 2021;24:312-325. https://doi.org/10.1038/s41593-020-00783-4

17. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Research. 2012;40:e133. https://doi.org/10.1093/nar/gks461

18. Hinderliter PM, Minard KR, Orr G, Chrisler WB, Thrall BD, Pounds JG, et al. ISDD: A computational model of particle sedimentation, diffusion and target cell dosimetry for in vitro toxicity studies. Particle and Fibre Toxicology. 2010;7:36. https://doi.org/10.1186/1743-8977-7-36

19. Franzén O, Gan LM, Björkegren JLM. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. Database. 2019;2019:baz046. https://doi.org/10.1093/database/baz046

## Figure legends

**Figure 1. Within-series transcriptional responses at 100 µg/mL.** Each point represents a gene. The horizontal axis gives the apeglm-shrunken log₂ fold change relative to the series-specific control; the vertical axis gives −log₁₀ adjusted P (capped at 40). Grey points are genes not meeting the discovery criteria, orange points are discoveries that did not pass all stability checks, and red points are discoveries that also retained significance and direction in every block-omission fit and agreed in direction within every observed block. Dashed vertical lines mark the 1.5-fold effect-size criterion.

![](../results/figures/fig3_volcano.png){width=6.1in}

**Figure 2. Exploratory high-dose enrichment of astrocyte and cell-type marker sets.** Columns show the two within-series responses and their direct contrast. Fill gives the normalized enrichment score, point size reflects −log₁₀ adjusted P, and a black outline marks adjusted P < 0.05 in the exploratory panel. Positive values in the contrast column favour the 50 nm response. A1 and A2 refer to historical marker lists [12].

![](../results/figures/fig5_signatures.png){width=5.7in}

**Figure 3. Ciliary and complement responses across doses and model fits.** A. Each point is the mean difference between series in log₂ response for genes in a set, calculated within one inferred block. Lines join doses with the same suffix. B. Fill gives the normalized enrichment score. Dots mark BH-adjusted P < 0.05 among six tests per fit. Positive values favour the 50 nm response. The full and group-only models use all 32 libraries. Each omission fit uses 24.

![](../results/figures/fig7_dose_response.png){width=5.6in}

## Additional files

**Additional file 1.** Group mean normalized counts for astrocyte, oligodendrocyte-lineage, neuronal, progenitor and proliferation markers across the eight exposure groups (source: `results/explore/marker_group_means.tsv`).


**Additional file 2.** Frozen panel definitions and member-level source audit (data/marker_sources). External-marker enrichment and matched-background results are archived in results/marker_sensitivity and results/amplitude.
