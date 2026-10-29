# 07b_replicates.R — the replicate structure found in the PCA
#   PCA showed that replicate k of every group clusters with replicate k of the other groups,
#   in BOTH particle sizes. So the 4 replicates are 4 paired experiment runs shared by all groups.
#   This script:
#     1. adds a `replicate` column (01..04, from the GSM order within each group) to samples.csv
#     2. measures how much of each top PC is explained by replicate vs. group (R²)
#     3. re-tests control vs control with replicate as a blocking factor (paired test)
#     4. draws the PCA coloured by replicate
# Run from the project root:  source("scripts/07b_replicates.R", echo = FALSE)

suppressPackageStartupMessages({ library(tximport); library(DESeq2) })
out <- "results/explore"
sink(file.path(out, "replicates_log.txt"), split = TRUE)

# ---------------------------------------------------------------- 1. replicate column
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
samples$replicate <- ave(samples$gsm, samples$group,
                         FUN = function(g) sprintf("rep%d", rank(g)))   # lowest GSM in a group = rep1
cat("== 1. replicate assignment (GSM order within each group)\n")
print(table(samples$group, samples$replicate))
write.csv(samples, "samples.csv", row.names = FALSE, quote = FALSE)
cat("samples.csv updated with a 'replicate' column\n")

rownames(samples) <- samples$run
samples$group     <- factor(samples$group, levels = c("50A","50B","50C","50D","500A","500B","500C","500D"))
samples$batch     <- factor(samples$batch, levels = c("b50", "b500"))
samples$replicate <- factor(samples$replicate)

txi <- readRDS("results/tximport/txi.rds")
stopifnot(identical(colnames(txi$counts), samples$run))
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ replicate + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
v <- apply(mat, 1, var); sel <- order(v, decreasing = TRUE)[1:500]
p <- prcomp(t(mat[sel, ])); ve <- 100 * p$sdev^2 / sum(p$sdev^2)

# ---------------------------------------------------------------- 2. what drives each PC?
cat("\n== 2. Share of each PC explained by replicate vs group (R² of a one-factor linear model)\n")
r2 <- t(sapply(1:5, function(k) c(
  PC = k, pct_variance = round(ve[k], 1),
  R2_replicate = round(summary(lm(p$x[, k] ~ samples$replicate))$r.squared, 2),
  R2_group     = round(summary(lm(p$x[, k] ~ samples$group))$r.squared, 2),
  R2_size_batch = round(summary(lm(p$x[, k] ~ samples$batch))$r.squared, 2))))
print(r2, row.names = FALSE)

# ---------------------------------------------------------------- 3. paired control test
cat("\n== 3. Controls 500A vs 50A, now blocking on replicate (paired)\n")
ctrl <- samples$run[samples$dose_ugml == 0]
dc <- dds[, ctrl]; dc$batch <- droplevels(dc$batch); dc$replicate <- droplevels(dc$replicate)
design(dc) <- ~ replicate + batch
dc <- DESeq(dc, quiet = TRUE)
rc <- results(dc, contrast = c("batch", "b500", "b50"))
cat("genes padj<0.05:", sum(rc$padj < 0.05, na.rm = TRUE),
    "| also |FC|>=1.5:", sum(rc$padj < 0.05 & abs(rc$log2FoldChange) >= log2(1.5), na.rm = TRUE), "\n")
dr <- dds[, ctrl]; dr$replicate <- droplevels(dr$replicate)
design(dr) <- ~ replicate
dr <- DESeq(dr, test = "LRT", reduced = ~ 1, quiet = TRUE)
rr <- results(dr)
cat("for comparison — genes that differ BETWEEN REPLICATES among the 8 controls (LRT, padj<0.05):",
    sum(rr$padj < 0.05, na.rm = TRUE), "\n")

# ---------------------------------------------------------------- 4. PCA by replicate
rep_col <- c(rep1 = "#1b9e77", rep2 = "#7570b3", rep3 = "#d95f02", rep4 = "#e7298a")
dose_pch <- c("0" = 1, "1" = 2, "10" = 0, "100" = 16)
for (ab in list(c(1, 2), c(2, 3))) {
  a <- ab[1]; b <- ab[2]
  png(file.path(out, sprintf("pca_by_replicate_PC%d_PC%d.png", a, b)), width = 1700, height = 1400, res = 200)
  par(mar = c(5, 5, 4, 9), xpd = TRUE)
  plot(p$x[, a], p$x[, b], col = rep_col[as.character(samples$replicate)],
       pch = dose_pch[as.character(samples$dose_ugml)], cex = 1.8, lwd = 2,
       xlab = sprintf("PC%d (%.1f%%)", a, ve[a]), ylab = sprintf("PC%d (%.1f%%)", b, ve[b]),
       main = sprintf("PCA coloured by replicate: PC%d vs PC%d", a, b))
  text(p$x[, a], p$x[, b], labels = samples$group, pos = 3, cex = 0.55, col = "grey30")
  u <- par("usr")
  legend(u[2] + 0.02 * diff(u[1:2]), u[4], title = "replicate", legend = names(rep_col),
         col = rep_col, pch = 15, bty = "n", cex = 0.85)
  legend(u[2] + 0.02 * diff(u[1:2]), u[3] + 0.45 * diff(u[3:4]), title = "dose (ug/ml)",
         legend = names(dose_pch), pch = dose_pch, bty = "n", cex = 0.85)
  invisible(dev.off())
}
cat("\nsaved pca_by_replicate_PC1_PC2.png and pca_by_replicate_PC2_PC3.png\n")
sink()
