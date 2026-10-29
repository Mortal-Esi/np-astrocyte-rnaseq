# 12_authors_model.R — can we reproduce the AUTHORS' numbers by changing our model?
# Hypothesis (from 11_discordance.R): the authors (1) did not include the replicate (paired experiment) term and
# (2) reported SHRUNKEN fold changes. Without the replicate term, between-experiment differences become "noise";
# shrinkage then pulls noisy / low-count genes to FC ~ 1 — which would make the vertical band in fc_scatter.
# Test: refit with design ~ group (no replicate), with and without shrinkage, and see which version matches the authors best.
# Run from the project root:  source("scripts/12_authors_model.R", echo = FALSE)       (~5 min)
# Needs: results/deseq2/dds.rds and *_vs_*.tsv (08), data/authors/*.tsv.gz, results/tximport/gene_annotation.tsv
# Outputs: results/authors_model/ (authors_model_log.txt, model_comparison.tsv)

suppressPackageStartupMessages(library(DESeq2))
out <- "results/authors_model"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
sink(file.path(out, "authors_model_log.txt"), split = TRUE)

dds <- readRDS("results/deseq2/dds.rds")
ann_all <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
to_log2 <- function(fc) sign(fc) * log2(abs(fc))
read_authors <- function(f) { a <- read.delim(gzfile(f), stringsAsFactors = FALSE); a$EnsemblGene <- sub("\\.[0-9]+$", "", a$EnsemblGene); a }
auth <- list(`50`  = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_50.tsv.gz"),
             `500` = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_500.tsv.gz"))
match_ids <- function(a) {               # same rule as 09_validate.R / 11_discordance.R
  direct <- a$EnsemblGene %in% ann_all$gene
  uniq_ours <- ann_all$gene_name[!duplicated(ann_all$gene_name) & !duplicated(ann_all$gene_name, fromLast = TRUE)]
  uniq_auth <- a$Gene_name[!duplicated(a$Gene_name) & !duplicated(a$Gene_name, fromLast = TRUE)]
  by_name <- !direct & a$Gene_name %in% uniq_ours & a$Gene_name %in% uniq_auth & !grepl("^ENSRNOG|^AABR|^AC[0-9]", a$Gene_name)
  id <- ifelse(direct, a$EnsemblGene, NA)
  id[by_name] <- ann_all$gene[match(a$Gene_name[by_name], ann_all$gene_name)]
  clash <- by_name & id %in% id[direct]; id[clash] <- NA
  id
}

# ---------------------------------------------------------------- refit without the replicate term
cat("== Refit with design ~ group (no replicate term)\n")
dds_nr <- dds
design(dds_nr) <- ~ group
dds_nr <- DESeq(dds_nr, quiet = TRUE)
dds_nr500 <- dds_nr
dds_nr500$group <- relevel(dds_nr500$group, ref = "500A")
dds_nr500 <- nbinomWaldTest(dds_nr500, quiet = TRUE)

variants <- function(d, coef_name) {
  r <- results(d, name = coef_name)
  ape <- tryCatch(lfcShrink(d, coef = coef_name, type = "apeglm", quiet = TRUE)$log2FoldChange,
                  error = function(e) { cat("apeglm failed:", conditionMessage(e), "\n"); rep(NA, nrow(r)) })
  nrm <- tryCatch(lfcShrink(d, coef = coef_name, type = "normal", quiet = TRUE)$log2FoldChange, error = function(e) rep(NA, nrow(r)))
  data.frame(gene_id = rownames(r), noRep_unshrunk = r$log2FoldChange, noRep_apeglm = ape, noRep_normal = nrm,
             noRep_padj = r$padj, stringsAsFactors = FALSE)
}
fits <- list(`50` = variants(dds_nr, "group_50D_vs_50A"), `500` = variants(dds_nr500, "group_500D_vs_500A"))

# ---------------------------------------------------------------- compare every version with the authors
rows <- list()
for (sz in names(auth)) {
  a <- auth[[sz]]
  A <- data.frame(gene_id = match_ids(a), authors_l2 = to_log2(a$FoldChange_D_vs_A), authors_padj = a$padj_D_vs_A)
  A <- A[!is.na(A$gene_id) & !duplicated(A$gene_id), ]
  o <- read.delim(file.path("results/deseq2", sprintf("%sD_vs_%sA.tsv", sz, sz)), stringsAsFactors = FALSE)
  o <- data.frame(gene_id = o$gene_id, baseMean = o$baseMean, paired_unshrunk = o$log2FC, paired_apeglm = o$log2FC_shrunk, paired_padj = o$padj)
  m <- merge(merge(A, o, by = "gene_id"), fits[[sz]], by = "gene_id")
  # the "vertical band": our paired analysis calls a DEG, authors report ~0
  band <- !is.na(m$paired_padj) & m$paired_padj < 0.05 & abs(m$paired_apeglm) >= log2(1.5) & abs(m$authors_l2) < 0.3
  auth_deg <- m$authors_padj < 0.05 & abs(m$authors_l2) >= log2(1.5)
  for (v in c("paired_unshrunk", "paired_apeglm", "noRep_unshrunk", "noRep_apeglm", "noRep_normal")) {
    x <- m[[v]]; ok <- !is.na(x)
    p <- if (startsWith(v, "paired")) m$paired_padj else m$noRep_padj
    deg <- !is.na(p) & p < 0.05 & abs(x) >= log2(1.5)
    rows[[length(rows) + 1]] <- data.frame(size = sz, version = v,
      pearson = round(cor(m$authors_l2[ok], x[ok]), 3),
      spearman = round(cor(m$authors_l2[ok], x[ok], method = "spearman"), 3),
      median_abs_diff = round(median(abs(m$authors_l2[ok] - x[ok])), 3),
      band_genes_median_abs_l2 = round(median(abs(x[band & ok])), 2),
      n_DEG = sum(deg), authors_DEG_recovered_pct = round(100 * sum(deg & auth_deg) / sum(auth_deg), 1),
      jaccard_with_authors = round(sum(deg & auth_deg) / sum(deg | auth_deg), 3))
  }
  cat(sprintf("\n%s nm: %d shared genes, %d 'vertical band' genes\n", sz, nrow(m), sum(band)))
  cat(sprintf("Spearman of -log10 padj with the authors: paired %.3f | no replicate %.3f\n",
              cor(-log10(m$authors_padj), -log10(m$paired_padj), use = "complete.obs", method = "spearman"),
              cor(-log10(m$authors_padj), -log10(m$noRep_padj), use = "complete.obs", method = "spearman")))
}
res <- do.call(rbind, rows)
cat("\n== Which version of our analysis looks most like the authors' table?\n")
cat("(band_genes_median_abs_l2: for the 'vertical band' genes, how big is the fold change in this version — ~0 means this version reproduces the band)\n")
print(res, row.names = FALSE)
write.table(res, file.path(out, "model_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- is Npy hiding under another name?
cat("\n== Authors' genes with a dose-100 fold change >= 16x (either direction) — looking for Npy under another ID/name\n")
for (sz in names(auth)) {
  a <- auth[[sz]]; big <- a[abs(a$FoldChange_D_vs_A) >= 16, c("EnsemblGene", "Gene_name", "FoldChange_D_vs_A", "padj_D_vs_A")]
  big$in_our_annotation <- big$EnsemblGene %in% ann_all$gene
  cat(sprintf("\n%s nm (%d genes):\n", sz, nrow(big)))
  print(big[order(-abs(big$FoldChange_D_vs_A)), ], row.names = FALSE)
}
sink()
cat("done — see", out, "\n")
