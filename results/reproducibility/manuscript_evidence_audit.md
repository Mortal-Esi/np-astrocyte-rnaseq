# Manuscript evidence audit, 20 November 2026

The 21 references in the audited manuscript were checked by identifier and title. Reference 22 was added after publication of the v1.3.1 archive. The table retains the numbering of that audit. In the submission manuscript, FastQC is reference 6 and the subsequent numbers increase by one. The checks below distinguish the existence of a work from its relevance to the claim. Software publications describe the methods; the versions actually run are documented in the project environment and session logs.

| Reference | Record checked | Use in the manuscript |
|---|---|---|
| 1 | [Publisher](https://doi.org/10.3390/nano14050429) | Source experiment, exposure duration, published counts, overlap and metabolic assay |
| 2 | [NCBI GEO series metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256038) | Submitted sample series and original differential-expression tables |
| 3 | [Crossref](https://doi.org/10.4049/jimmunol.139.7.2361) | Mouse astrocyte C3 and factor B production |
| 4 | [Crossref](https://doi.org/10.3390/ijms262311273) | Human astrocyte response to polystyrene under different experimental conditions |
| 5 | [Crossref](https://doi.org/10.1242/dev.041129) | Postnatal mouse Foxj1 and cell identity |
| 6 | [Crossref](https://doi.org/10.1093/bioinformatics/bty560) | fastp quality reports |
| 7 | [Crossref](https://doi.org/10.1093/bioinformatics/btw354) | MultiQC aggregation |
| 8 | [Crossref](https://doi.org/10.1038/nmeth.4197) | Salmon quantification and bias correction |
| 9 | [Crossref](https://doi.org/10.12688/f1000research.7563.2) | tximport gene-level estimates and offsets |
| 10 | [Crossref](https://doi.org/10.1186/s13059-014-0550-8) | DESeq2 model and inference |
| 11 | [Crossref](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) | BH multiple-testing adjustment |
| 12 | [Crossref](https://doi.org/10.1093/bioinformatics/bty895) | apeglm shrinkage |
| 13 | [Crossref](https://doi.org/10.1038/nature21029) | Historical mouse reactive signatures, with a separate member-level source audit |
| 14 | [Crossref](https://doi.org/10.1073/pnas.0506580102) | Preranked enrichment framework |
| 15 | [Crossref](https://doi.org/10.1016/j.xinn.2021.100141) | clusterProfiler software |
| 16 | [Crossref](https://doi.org/10.1093/database/baz046) | PanglaoDB resource; 2020 marker file identified separately |
| 17 | [Crossref](https://doi.org/10.1038/s41593-020-00783-4) | Functional interpretation of reactive astrocyte signatures |
| 18 | [Crossref](https://doi.org/10.1093/nar/gks461) | Why between-gene correlation matters in competitive tests |
| 19 | [Crossref](https://doi.org/10.1186/1743-8977-7-36) | Particle transport and delivered dose in culture |
| 20 | [DataCite and archived Zenodo metadata](https://doi.org/10.5281/zenodo.24197910) | Earlier project deposit |
| 21 | [DataCite and archived Zenodo metadata](https://doi.org/10.5281/zenodo.24361377) | Processed data, code and archived analysis |
| 22 | [Published Zenodo record](https://doi.org/10.5281/zenodo.24390580) | Current archive. The deposited ZIP matches the local v1.3.1 bundle in size and MD5 checksum. |

## Corrections and checks

- The four high-dose 500 nm PC1 scores in results/explore/pca_all_scores.tsv are 18.59, 17.58, 35.52 and 22.84. Sample 500C-03 scores 22.11. The published draft's lower bound of 18.6 has been corrected to 17.6.
- The diagnostic in script 11 estimates similarity to the two particle-series profiles learned from the other three blocks. Seven of eight high-dose libraries had a larger coefficient for their recorded series. It cannot verify the concentration or authenticate the identity of a library.
- Scripts 17–19 test enrichment across six overlapping fits. Script 22 matches genes using the full blocked model only. The manuscript now explicitly distinguishes these procedures and limits the matched result to the set mean. The larger external astrocyte set does not establish a distinct biological programme.
- In script 09 the published list is defined from deposited adjusted P values and signed fold changes. The reanalysis uses apeglm-shrunken fold changes. Recovery equals 474/497 and 1,122/1,156 among comparable published discoveries. FastQC and the dated PanglaoDB file now have direct source links.
- Script verify_results.py passed on the stored gene and enrichment tables. These reproduce 990, 926, 144, 89, 50 and 852, the six dose-by-set stability counts and the three matched excesses. The 18 per-gene estimates yield the unshrunken log₂ means +0.319 and −0.269 and the direct difference +0.589. The stored matched background is +0.132 (excess +0.457, BH-adjusted P = 0.00015).
- All three sets were chosen after inspection. Matching uses an estimated component of the contrast and ignores inter-gene dependence. Its empirical P value is conditional on this resampling procedure. The study measures RNA abundance in bulk cultures at equal nominal mass and supplies no direct test of complement activity, cell-type proportions or an isolated diameter effect.

This review checked source records, code, archived outputs and selected calculations. No new FASTQ processing or model fitting was performed. The earlier restart from tximport is documented separately in the processed-data audit. Historical list-selection timing relies in part on the author's account.
