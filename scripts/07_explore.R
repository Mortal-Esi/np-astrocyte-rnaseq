# 07_explore.R — exploratory analysis before differential expression
#   1. filter + variance-stabilising transform (VST)
#   2. PCA of all samples and of the 8 controls only (batch question)
#   3. control-vs-control DE (how big is the batch effect, in genes?)
#   4. sample-to-sample distance heatmap
#   5. replicate consistency (outlier check)
#   6. cell-composition markers across groups (are oligodendrocyte-lineage cells changing?)
#
# Run from the project root (RStudio project open):
#   source("scripts/07_explore.R", echo = TRUE)
# Outputs: results/explore/*.png, *.tsv and explore_log.txt

suppressPackageStartupMessages({ library(tximport); library(DESeq2) })

out <- "results/explore"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out, "explore_log.txt"), split = TRUE)

# ---------------------------------------------------------------- inputs
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
txi <- readRDS("results/tximport/txi.rds")
stopifnot("sample order differs between samples.csv and txi" = identical(colnames(txi$counts), samples$run))
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)

rownames(samples) <- samples$run
samples$group <- factor(samples$group, levels = c("50A", "50B", "50C", "50D", "500A", "500B", "500C", "500D"))
samples$batch <- factor(samples$batch, levels = c("b50", "b500"))
samples$dose  <- factor(samples$dose_ugml, levels = c(0, 1, 10, 100))

# ---------------------------------------------------------------- 1. filter + VST
cat("== 1. Filter and transform\n")
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)
cat("genes kept (>=10 reads in >=4 samples):", nrow(dds), "of", length(keep), "\n")
vsd <- vst(dds, blind = TRUE)          # blind = ignore the design: unbiased look at the data
mat <- assay(vsd)

# ---------------------------------------------------------------- helpers
dose_col <- c("0" = "grey55", "1" = "#9ecae1", "10" = "#3182bd", "100" = "#de2d26")
batch_pch <- c(b50 = 16, b500 = 17)
run_pca <- function(m, ntop = 500) {
  v <- apply(m, 1, var)
  sel <- order(v, decreasing = TRUE)[seq_len(min(ntop, length(v)))]
  p <- prcomp(t(m[sel, ]), center = TRUE, scale. = FALSE)
  list(x = p$x, ve = 100 * p$sdev^2 / sum(p$sdev^2))
}
plot_pca <- function(p, meta, file, title, a = 1, b = 2) {
  png(file, width = 1700, height = 1400, res = 200)
  par(mar = c(5, 5, 4, 9), xpd = TRUE)
  plot(p$x[, a], p$x[, b], col = dose_col[as.character(meta$dose_ugml)],
       pch = batch_pch[as.character(meta$batch)], cex = 1.8,
       xlab = sprintf("PC%d (%.1f%%)", a, p$ve[a]), ylab = sprintf("PC%d (%.1f%%)", b, p$ve[b]),
       main = title)
  text(p$x[, a], p$x[, b], labels = meta$group, pos = 3, cex = 0.55, col = "grey30")
  u <- par("usr")
  legend(u[2] + 0.02 * diff(u[1:2]), u[4], title = "dose (ug/ml)", legend = names(dose_col),
         col = dose_col, pch = 15, bty = "n", cex = 0.85)
  legend(u[2] + 0.02 * diff(u[1:2]), u[3] + 0.45 * diff(u[3:4]), title = "batch / size",
         legend = c("b50 (50 nm)", "b500 (500 nm)"), pch = batch_pch, bty = "n", cex = 0.85)
  invisible(dev.off())
}

# ---------------------------------------------------------------- 2. PCA
cat("\n== 2. PCA (top 500 most variable genes)\n")
p_all <- run_pca(mat)
cat("variance explained, PC1..PC6 (%):", round(p_all$ve[1:6], 1), "\n")
plot_pca(p_all, samples, file.path(out, "pca_all_PC1_PC2.png"), "All samples: PC1 vs PC2")
plot_pca(p_all, samples, file.path(out, "pca_all_PC3_PC4.png"), "All samples: PC3 vs PC4", 3, 4)
pc_tab <- data.frame(samples[, c("run", "group", "batch", "dose_ugml")], round(p_all$x[, 1:4], 2))
write.table(pc_tab, file.path(out, "pca_all_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("mean PC1/PC2 per group:\n")
print(round(aggregate(cbind(PC1, PC2) ~ group, data = pc_tab, FUN = mean)[, -1], 1) |>
        `rownames<-`(levels(samples$group)))

ctrl <- samples$run[samples$dose_ugml == 0]
p_ctrl <- run_pca(mat[, ctrl])
cat("\ncontrols only — variance explained PC1..PC4 (%):", round(p_ctrl$ve[1:4], 1), "\n")
plot_pca(p_ctrl, samples[ctrl, ], file.path(out, "pca_controls.png"), "Controls only (50A vs 500A)")
cat("controls PC1 scores:\n"); print(round(setNames(p_ctrl$x[, 1], samples[ctrl, "group"]), 1))

# ---------------------------------------------------------------- 3. control vs control
cat("\n== 3. Batch effect in genes: 500A controls vs 50A controls\n")
dds_c <- dds[, ctrl]
dds_c$batch <- droplevels(dds_c$batch)
design(dds_c) <- ~ batch
dds_c <- DESeq(dds_c, quiet = TRUE)
res_c <- results(dds_c, contrast = c("batch", "b500", "b50"))
n05 <- sum(res_c$padj < 0.05, na.rm = TRUE)
n15 <- sum(res_c$padj < 0.05 & abs(res_c$log2FoldChange) >= log2(1.5), na.rm = TRUE)
cat("control-vs-control genes with padj<0.05:", n05, "| also |FC|>=1.5:", n15, "\n")
cat("(for scale: the authors report 531 and 1274 DEGs for dose 100 vs its own control)\n")
top_c <- head(res_c[order(res_c$padj), ], 15)
top_c <- data.frame(gene_id = rownames(top_c), gene_name = ann$gene_name[match(rownames(top_c), ann$gene)],
                    log2FC = round(top_c$log2FoldChange, 2), padj = signif(top_c$padj, 3))
cat("top control-vs-control genes:\n"); print(top_c, row.names = FALSE)
write.table(as.data.frame(res_c), file.path(out, "controls_b500_vs_b50.tsv"), sep = "\t", quote = FALSE)

# ---------------------------------------------------------------- 4. distance heatmap
cat("\n== 4. Sample distance heatmap\n")
d <- as.matrix(dist(t(mat)))
lab <- paste(samples$group, samples$run, sep = " | ")
dimnames(d) <- list(lab, lab)
png(file.path(out, "sample_distance_heatmap.png"), width = 2000, height = 2000, res = 200)
heatmap(d, symm = TRUE, col = hcl.colors(50, "Blues", rev = FALSE), margins = c(11, 11),
        ColSideColors = dose_col[as.character(samples$dose_ugml)],
        RowSideColors = ifelse(samples$batch == "b50", "black", "orange"),
        main = "Euclidean distance (VST)")
invisible(dev.off())
cat("saved sample_distance_heatmap.png (row colour: black = b50, orange = b500; column colour = dose)\n")

# ---------------------------------------------------------------- 5. replicate consistency
cat("\n== 5. Replicate consistency (outlier check)\n")
cors <- cor(mat)
rep_cor <- sapply(samples$run, function(s) {
  mates <- samples$run[samples$group == samples[s, "group"] & samples$run != s]
  mean(cors[s, mates])
})
rc <- data.frame(run = samples$run, group = samples$group, mean_cor_with_replicates = round(rep_cor, 4))
rc <- rc[order(rc$mean_cor_with_replicates), ]
cat("lowest 5 (a clear outlier would stand well below the rest):\n"); print(head(rc, 5), row.names = FALSE)
cat("range:", round(range(rep_cor), 4), "\n")

# ---------------------------------------------------------------- 6. cell composition markers
cat("\n== 6. Cell-type markers across groups (DESeq2-normalised counts, group means)\n")
markers <- c(Gfap = "astro", Aqp4 = "astro", Slc1a3 = "astro", Aldh1l1 = "astro",
             Pdgfra = "OPC", Cspg4 = "OPC", Olig1 = "OL-lineage", Olig2 = "OL-lineage", Sox10 = "OL-lineage",
             Mbp = "mature OL", Plp1 = "mature OL", Mog = "mature OL",
             Rbfox3 = "neuron", Tubb3 = "neuron", Nes = "progenitor", Sox2 = "progenitor", Mki67 = "proliferation")
norm <- counts(dds, normalized = TRUE)
ids <- ann$gene[match(names(markers), ann$gene_name)]
ok <- !is.na(ids) & ids %in% rownames(norm)
grp_mean <- sapply(levels(samples$group), function(g) rowMeans(norm[ids[ok], samples$group == g, drop = FALSE]))
rownames(grp_mean) <- paste0(names(markers)[ok], " (", markers[ok], ")")
print(round(grp_mean))
lfc <- data.frame(
  `50nm_D_vs_A`  = log2((grp_mean[, "50D"] + 1) / (grp_mean[, "50A"] + 1)),
  `500nm_D_vs_A` = log2((grp_mean[, "500D"] + 1) / (grp_mean[, "500A"] + 1)), check.names = FALSE)
cat("\nlog2(dose 100 / control) — rough look, not a test:\n"); print(round(lfc, 2))
if (any(!ok)) cat("markers not found:", names(markers)[!ok], "\n")
write.table(round(grp_mean, 1), file.path(out, "marker_group_means.tsv"), sep = "\t", quote = FALSE)

sink()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("done — see", out, "\n")
