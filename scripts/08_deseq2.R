# 08_deseq2.R — differential expression with DESeq2, blocking on replicate
#
# Design: ~ replicate + group
#   replicate = the 4 paired experiment runs shared by all 8 groups (found in 07_explore / 07b)
#   group     = 50A..50D, 500A..500D  (A = control, B/C/D = 1/10/100 ug/ml)
#
# Contrasts:
#   dose vs its own control, per particle size: 50B/50C/50D vs 50A, 500B/500C/500D vs 500A
#   same dose, different size:                  50D vs 500D
#   difference of responses (interaction):      (50D - 50A) - (500D - 500A)
# Plus a sensitivity check: the same 50D/500D contrasts WITHOUT replicate in the model.
#
# Run from the project root:  source("scripts/08_deseq2.R", echo = FALSE)
# Outputs: results/deseq2/<contrast>.tsv, deg_summary.tsv, ma_*.png, volcano_*.png, deseq2_log.txt

suppressPackageStartupMessages({ library(tximport); library(DESeq2) })
out <- "results/deseq2"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out, "deseq2_log.txt"), split = TRUE)

shrink_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
cat("LFC shrinkage method:", shrink_type, "\n")

# ---------------------------------------------------------------- inputs
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
if (!"replicate" %in% names(samples)) stop("samples.csv has no 'replicate' column — run scripts/07b_replicates.R first")
rownames(samples) <- samples$run
samples$group     <- factor(samples$group, levels = c("50A","50B","50C","50D","500A","500B","500C","500D"))
samples$replicate <- factor(samples$replicate)
txi <- readRDS("results/tximport/txi.rds")
stopifnot(identical(colnames(txi$counts), samples$run))
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------- fit
cat("\n== 1. Fit: ~ replicate + group\n")
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ replicate + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
cat("genes tested:", nrow(dds), "\n")
dds <- DESeq(dds, quiet = TRUE)                        # reference level = 50A
dds500 <- dds; dds500$group <- relevel(dds500$group, "500A")
dds500 <- nbinomWaldTest(dds500, quiet = TRUE)         # same dispersions, coefficients vs 500A
cat("coefficients (ref 50A):", paste(resultsNames(dds), collapse = ", "), "\n")

# ---------------------------------------------------------------- helpers
signed_fc <- function(l2) ifelse(l2 >= 0, 2^l2, -2^(-l2))   # authors' format: 2 = doubled, -2 = halved
tidy <- function(res, shr) {
  df <- data.frame(gene_id = rownames(res),
                   gene_name = ann$gene_name[match(rownames(res), ann$gene)],
                   gene_biotype = ann$gene_biotype[match(rownames(res), ann$gene)],
                   baseMean = round(res$baseMean, 2),
                   log2FC = round(res$log2FoldChange, 4),
                   lfcSE = round(res$lfcSE, 4),
                   log2FC_shrunk = round(shr$log2FoldChange, 4),
                   FC_shrunk_signed = round(signed_fc(shr$log2FoldChange), 4),
                   pvalue = signif(res$pvalue, 4),
                   padj = signif(res$padj, 4),
                   stringsAsFactors = FALSE)
  df[order(df$padj), ]
}
summarise <- function(name, df, res) {
  s <- df$padj < 0.05 & !is.na(df$padj)
  a <- s & abs(df$FC_shrunk_signed) >= 1.5
  data.frame(contrast = name, tested = sum(!is.na(df$padj)),
             padj05 = sum(s), up = sum(s & df$log2FC > 0), down = sum(s & df$log2FC < 0),
             padj05_FC1.5 = sum(a),
             cooks_outlier_genes = sum(is.na(res$pvalue) & res$baseMean > 0),
             low_count_filtered = sum(!is.na(res$pvalue) & is.na(res$padj)))
}
plot_ma_volcano <- function(name, df) {
  sig <- !is.na(df$padj) & df$padj < 0.05
  png(file.path(out, paste0("ma_", name, ".png")), width = 1600, height = 1200, res = 200)
  plot(log10(df$baseMean + 1), df$log2FC_shrunk, pch = 16, cex = 0.35,
       col = ifelse(sig, "#de2d26", "grey70"), xlab = "log10(mean normalised count + 1)",
       ylab = "shrunken log2 fold change", main = paste("MA:", name))
  abline(h = 0, col = "grey30"); invisible(dev.off())
  png(file.path(out, paste0("volcano_", name, ".png")), width = 1600, height = 1400, res = 200)
  y <- -log10(pmax(df$padj, 1e-300)); y[is.na(y)] <- 0
  plot(df$log2FC_shrunk, y, pch = 16, cex = 0.4, col = ifelse(sig, "#de2d26", "grey70"),
       xlab = "shrunken log2 fold change", ylab = "-log10(padj)", main = paste("Volcano:", name))
  top <- head(df[sig, ], 15)
  if (nrow(top)) text(top$log2FC_shrunk, -log10(pmax(top$padj, 1e-300)), top$gene_name, cex = 0.5, pos = 3)
  abline(h = -log10(0.05), lty = 2, col = "grey40"); invisible(dev.off())
}
run_coef <- function(d, coef, name) {
  res <- results(d, name = coef, alpha = 0.05)
  shr <- lfcShrink(d, coef = coef, type = shrink_type, quiet = TRUE)
  df <- tidy(res, shr)
  write.table(df, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_ma_volcano(name, df)
  summarise(name, df, res)
}
run_contrast <- function(d, contrast, name, shrink = TRUE) {
  res <- results(d, contrast = contrast, alpha = 0.05)
  # apeglm needs a single coefficient; for contrasts use "normal" shrinkage.
  # The interaction contrast cannot be shrunk this way → unshrunken LFC is reported in both columns.
  shr <- if (shrink) lfcShrink(d, contrast = contrast, res = res, type = "normal", quiet = TRUE) else res
  df <- tidy(res, shr)
  write.table(df, file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_ma_volcano(name, df)
  summarise(name, df, res)
}

# ---------------------------------------------------------------- 2. contrasts
cat("\n== 2. Contrasts\n")
summ <- rbind(
  run_coef(dds,    "group_50B_vs_50A",   "50B_vs_50A"),
  run_coef(dds,    "group_50C_vs_50A",   "50C_vs_50A"),
  run_coef(dds,    "group_50D_vs_50A",   "50D_vs_50A"),
  run_coef(dds500, "group_500B_vs_500A", "500B_vs_500A"),
  run_coef(dds500, "group_500C_vs_500A", "500C_vs_500A"),
  run_coef(dds500, "group_500D_vs_500A", "500D_vs_500A"),
  run_contrast(dds, c("group", "50D", "500D"), "50D_vs_500D"),
  run_contrast(dds, list(c("group_50D_vs_50A", "group_500A_vs_50A"), c("group_500D_vs_50A")),
               "response50_minus_response500", shrink = FALSE)
)
print(summ, row.names = FALSE)
write.table(summ, file.path(out, "deg_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("note: response50_minus_response500 uses UNSHRUNKEN log2FC (no shrinkage possible for this contrast)\n")
cat("\n(authors, same criterion padj<0.05 & |FC|>=1.5, dose 100: 50 nm = 531, 500 nm = 1274)\n")

# ---------------------------------------------------------------- 3. sensitivity: no replicate term
cat("\n== 3. Sensitivity check: same dose-100 contrasts WITHOUT replicate in the model\n")
d0 <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
d0 <- d0[rownames(dds), ]
d0 <- DESeq(d0, quiet = TRUE)
r50  <- results(d0, contrast = c("group", "50D", "50A"), alpha = 0.05)
r500 <- results(d0, contrast = c("group", "500D", "500A"), alpha = 0.05)
cat("padj<0.05 without replicate term — 50D vs 50A:", sum(r50$padj < 0.05, na.rm = TRUE),
    "| 500D vs 500A:", sum(r500$padj < 0.05, na.rm = TRUE), "\n")
cat("padj<0.05 with replicate term    — 50D vs 50A:", summ$padj05[summ$contrast == "50D_vs_50A"],
    "| 500D vs 500A:", summ$padj05[summ$contrast == "500D_vs_500A"], "\n")

# ---------------------------------------------------------------- 4. top genes
cat("\n== 4. Top 15 genes (by padj) for dose 100\n")
for (nm in c("50D_vs_50A", "500D_vs_500A")) {
  df <- read.delim(file.path(out, paste0(nm, ".tsv")))
  cat("\n--", nm, "\n"); print(head(df[, c("gene_name", "baseMean", "log2FC_shrunk", "FC_shrunk_signed", "padj")], 15), row.names = FALSE)
}

saveRDS(dds, file.path(out, "dds.rds"))
sink()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("done — see", out, "\n")
