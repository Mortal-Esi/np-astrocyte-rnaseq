# 10_robustness.R — are the key results robust?
#   1. per-sample expression of key genes, with each replicate's samples connected (paired view)
#   2. per-replicate dose-100 fold changes (does every replicate agree?)
#   3. Snca in detail (expression level, all contrasts)
#   4. leave-one-replicate-out: refit DESeq2 without each replicate in turn
#   5. low-count filter sensitivity: refit with almost no filtering
# Run from the project root:  source("scripts/10_robustness.R", echo = FALSE)
# Needs results/deseq2/dds.rds from 08_deseq2.R. Takes ~10 min (5 extra DESeq2 fits).

suppressPackageStartupMessages({ library(tximport); library(DESeq2) })
out <- "results/robustness"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out, "robustness_log.txt"), split = TRUE)

dds <- readRDS("results/deseq2/dds.rds")
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
key <- c("Snca", "C3", "Gfap", "Aqp4", "Vim", "Nes", "Npy", "Ccl2", "Gpnmb", "Lgals3", "Clu", "Abat")
key_id <- setNames(ann$gene[match(key, ann$gene_name)], key)
key_id <- key_id[!is.na(key_id) & key_id %in% rownames(dds)]
cat("key genes found:", paste(names(key_id), collapse = ", "), "\n")
cd <- as.data.frame(colData(dds))
norm <- counts(dds, normalized = TRUE)
rep_col <- c(rep1 = "#1b9e77", rep2 = "#7570b3", rep3 = "#d95f02", rep4 = "#e7298a")

# ---------------------------------------------------------------- 1. paired per-sample plots
x_of <- c("50A" = 1, "50B" = 2, "50C" = 3, "50D" = 4, "500A" = 6, "500B" = 7, "500C" = 8, "500D" = 9)
png(file.path(out, "key_genes_paired.png"), width = 2400, height = 1800, res = 170)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (g in names(key_id)) {
  y <- log2(norm[key_id[g], ] + 1)
  plot(NA, xlim = c(0.5, 9.5), ylim = range(y), xaxt = "n", xlab = "", ylab = "log2(normalised count + 1)", main = g)
  axis(1, at = x_of, labels = names(x_of), las = 2, cex.axis = 0.75)
  abline(v = 5, col = "grey85")
  for (r in levels(cd$replicate)) for (sz in c("50", "500")) {
    idx <- which(cd$replicate == r & grepl(paste0("^", sz, "[ABCD]$"), cd$group))
    idx <- idx[order(x_of[as.character(cd$group[idx])])]
    lines(x_of[as.character(cd$group[idx])], y[idx], col = rep_col[r], lwd = 1.5)
    points(x_of[as.character(cd$group[idx])], y[idx], col = rep_col[r], pch = 16)
  }
  # legend inside the first panel: a separate plot.new() would start a 2nd page and overwrite the PNG
  if (g == names(key_id)[1]) legend("topleft", legend = names(rep_col), col = rep_col, lwd = 2, pch = 16, cex = 0.7, bty = "n")
}
invisible(dev.off())
cat("saved key_genes_paired.png (one line per replicate; left = 50 nm, right = 500 nm)\n")

# ---------------------------------------------------------------- 2. per-replicate fold changes
cat("\n== 2. Dose 100 vs control, separately in each replicate (log2 of normalised counts)\n")
per_rep <- do.call(rbind, lapply(names(key_id), function(g) {
  v <- sapply(c("50", "500"), function(sz) sapply(levels(cd$replicate), function(r) {
    d <- norm[key_id[g], cd$replicate == r & cd$group == paste0(sz, "D")]
    a <- norm[key_id[g], cd$replicate == r & cd$group == paste0(sz, "A")]
    log2((d + 0.5) / (a + 0.5))
  }))
  data.frame(gene = g, size = c("50", "500"), t(round(v, 2)),
             all_same_sign = apply(v, 2, function(z) all(z > 0) || all(z < 0)), check.names = FALSE)
}))
print(per_rep, row.names = FALSE)
write.table(per_rep, file.path(out, "per_replicate_dose100_log2FC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- 3. Snca in detail
cat("\n== 3. Snca in detail\n")
if ("Snca" %in% names(key_id)) {
  sid <- key_id["Snca"]
  cat("mean normalised count, controls:", round(mean(norm[sid, cd$dose_ugml == 0])), "\n")
  cat("for scale — Gfap:", round(mean(norm[key_id["Gfap"], cd$dose_ugml == 0])),
      "| median over all tested genes:", round(median(rowMeans(norm[, cd$dose_ugml == 0]))), "\n")
  for (nm in c("50B_vs_50A", "50C_vs_50A", "50D_vs_50A", "500B_vs_500A", "500C_vs_500A", "500D_vs_500A",
               "50D_vs_500D", "response50_minus_response500")) {
    f <- file.path("results/deseq2", paste0(nm, ".tsv"))
    if (!file.exists(f)) f <- paste0(f, ".gz")
    t <- read.delim(f); r <- t[t$gene_id == sid, ]
    cat(sprintf("  %-30s log2FC %6.2f  shrunk %6.2f  padj %s\n", nm, r$log2FC, r$log2FC_shrunk, signif(r$padj, 3)))
  }
} else cat("Snca not among tested genes\n")

# ---------------------------------------------------------------- helper: refit and summarise
fit_summary <- function(d, label) {
  d <- DESeq(d, quiet = TRUE)
  r50  <- results(d, contrast = c("group", "50D", "50A"), alpha = 0.05)
  r500 <- results(d, contrast = c("group", "500D", "500A"), alpha = 0.05)
  row <- data.frame(fit = label, genes = nrow(d),
                    DEG_50D = sum(r50$padj < 0.05, na.rm = TRUE), DEG_500D = sum(r500$padj < 0.05, na.rm = TRUE))
  for (g in names(key_id)) {
    row[[paste0(g, "_50D_l2fc")]]   <- round(r50[key_id[g], "log2FoldChange"], 2)
    row[[paste0(g, "_50D_padj")]]   <- signif(r50[key_id[g], "padj"], 2)
    row[[paste0(g, "_500D_l2fc")]]  <- round(r500[key_id[g], "log2FoldChange"], 2)
    row[[paste0(g, "_500D_padj")]]  <- signif(r500[key_id[g], "padj"], 2)
  }
  list(row = row, sig50 = rownames(r50)[which(r50$padj < 0.05)], sig500 = rownames(r500)[which(r500$padj < 0.05)])
}
main <- list(sig50 = rownames(dds)[which(results(dds, contrast = c("group", "50D", "50A"))$padj < 0.05)],
             sig500 = rownames(dds)[which(results(dds, contrast = c("group", "500D", "500A"))$padj < 0.05)])
overlap <- function(a, b) round(100 * length(intersect(a, b)) / length(b), 1)

# ---------------------------------------------------------------- 4. leave one replicate out
cat("\n== 4. Leave-one-replicate-out (24 samples per fit; unshrunken log2FC)\n")
txi <- readRDS("results/tximport/txi.rds")
loo <- list()
for (r in levels(cd$replicate)) {
  keep_s <- cd$run[cd$replicate != r]
  d <- DESeqDataSetFromTximport(lapply(txi[c("abundance", "counts", "length")], function(m) m[rownames(dds), keep_s]) |>
                                  c(list(countsFromAbundance = txi$countsFromAbundance)),
                                colData = droplevels(cd[keep_s, ]), design = ~ replicate + group)
  fs <- fit_summary(d, paste("without", r))
  fs$row$pct_main_DEG50_kept  <- overlap(fs$sig50, main$sig50)
  fs$row$pct_main_DEG500_kept <- overlap(fs$sig500, main$sig500)
  loo[[r]] <- fs$row
  cat(" done:", r, "\n")
}
loo_tab <- do.call(rbind, loo)
print(t(loo_tab[, c("genes", "DEG_50D", "DEG_500D", "pct_main_DEG50_kept", "pct_main_DEG500_kept")]) |> `colnames<-`(loo_tab$fit))
cat("\nkey genes across leave-one-out fits (log2FC / padj):\n")
for (g in names(key_id)) cat(sprintf("  %-7s 50D: %s\n          500D: %s\n", g,
  paste(sprintf("%5.2f/%-8s", loo_tab[[paste0(g, "_50D_l2fc")]], loo_tab[[paste0(g, "_50D_padj")]]), collapse = " "),
  paste(sprintf("%5.2f/%-8s", loo_tab[[paste0(g, "_500D_l2fc")]], loo_tab[[paste0(g, "_500D_padj")]]), collapse = " ")))
write.table(loo_tab, file.path(out, "leave_one_replicate_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- 5. filter sensitivity
cat("\n== 5. Low-count filter sensitivity: keep every gene with >= 10 reads in total (vs >= 10 in >= 4 samples)\n")
d_all <- DESeqDataSetFromTximport(txi, colData = cd[colnames(txi$counts), ], design = ~ replicate + group)
d_all <- d_all[rowSums(counts(d_all)) >= 10, ]
fa <- fit_summary(d_all, "minimal filter")
cat("genes tested:", fa$row$genes, "(main analysis:", nrow(dds), ")\n")
cat("DEG 50D:", fa$row$DEG_50D, "(main", length(main$sig50), ") | DEG 500D:", fa$row$DEG_500D, "(main", length(main$sig500), ")\n")
cat("main-analysis DEGs still significant — 50D:", overlap(fa$sig50, main$sig50), "% | 500D:", overlap(fa$sig500, main$sig500), "%\n")
cat("new DEGs that the main filter had removed — 50D:", length(setdiff(fa$sig50, rownames(dds))),
    "| 500D:", length(setdiff(fa$sig500, rownames(dds))), "\n")
if ("Snca" %in% names(key_id)) cat("Snca with minimal filter — 50D padj:", fa$row$Snca_50D_padj, "| 500D padj:", fa$row$Snca_500D_padj, "\n")
write.table(fa$row, file.path(out, "filter_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink()
cat("done — see", out, "\n")
