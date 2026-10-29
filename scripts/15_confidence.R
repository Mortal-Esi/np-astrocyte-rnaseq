# 15_confidence.R — which results survive the strictest checks?
#   1. Parkinson-related genes (Uchl1, Pink1, Gpnmb, ...) + Npy: per-replicate fold changes and a paired plot
#   2. Leave-one-replicate-out (4 refits on 24 samples), keeping the FULL result tables this time
#   3. A "high-confidence DEG" set for 50D and 500D. A gene qualifies only if it is
#        a DEG in the main analysis (padj < 0.05 and |FC| >= 1.5, shrunken)
#      + significant (padj < 0.05) in ALL 4 leave-one-out fits, with the same direction
#      + changes in the same direction inside ALL 4 replicates (dose 100 vs its own control)
#   4. Verdict per key gene: robust / fragile
# Run from the project root:  source("scripts/15_confidence.R", echo = FALSE)      (~10 min: 4 DESeq2 fits)
# Needs: results/deseq2/dds.rds and *.tsv (08), results/tximport/txi.rds
# Outputs: results/confidence/ (confidence_log.txt, key_genes_verdict.tsv, high_confidence_DEGs_<size>.tsv, pd_genes_paired.png)

suppressPackageStartupMessages({ library(tximport); library(DESeq2) })
out <- "results/confidence"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
sink(file.path(out, "confidence_log.txt"), split = TRUE)

dds <- readRDS("results/deseq2/dds.rds")
cd <- as.data.frame(colData(dds))
norm <- counts(dds, normalized = TRUE)
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
key <- c("Gpnmb", "Uchl1", "Pink1", "Fbxo7", "Dnajc6", "Park7", "Lrrk2", "Snca", "Npy")
key_id <- setNames(ann$gene[match(key, ann$gene_name)], key)
key_id <- key_id[!is.na(key_id) & key_id %in% rownames(dds)]
reps <- levels(cd$replicate)
main_tab <- list(`50` = read.delim("results/deseq2/50D_vs_50A.tsv", stringsAsFactors = FALSE),
                 `500` = read.delim("results/deseq2/500D_vs_500A.tsv", stringsAsFactors = FALSE))

# ---------------------------------------------------------------- 1. per-replicate fold changes (all genes) + paired plot
cat("== 1. Dose 100 vs control inside each replicate\n")
per_rep <- lapply(c(`50` = "50", `500` = "500"), function(sz) sapply(reps, function(r) {
  d <- norm[, cd$replicate == r & cd$group == paste0(sz, "D")]; a <- norm[, cd$replicate == r & cd$group == paste0(sz, "A")]
  log2((d + 0.5) / (a + 0.5))
}))
for (g in names(key_id)) for (sz in c("50", "500"))
  cat(sprintf("  %-7s %3s nm: %s\n", g, sz, paste(sprintf("%6.2f", per_rep[[sz]][key_id[g], ]), collapse = " ")))

x_of <- c("50A" = 1, "50B" = 2, "50C" = 3, "50D" = 4, "500A" = 6, "500B" = 7, "500C" = 8, "500D" = 9)
rep_col <- c(rep1 = "#1b9e77", rep2 = "#7570b3", rep3 = "#d95f02", rep4 = "#e7298a")
png(file.path(out, "pd_genes_paired.png"), width = 2000, height = 1800, res = 170)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (g in names(key_id)) {
  y <- log2(norm[key_id[g], ] + 1)
  plot(NA, xlim = c(0.5, 9.5), ylim = range(y), xaxt = "n", xlab = "", ylab = "log2(normalised count + 1)", main = g)
  axis(1, at = x_of, labels = names(x_of), las = 2, cex.axis = 0.75); abline(v = 5, col = "grey85")
  for (r in reps) for (sz in c("50", "500")) {
    idx <- which(cd$replicate == r & grepl(paste0("^", sz, "[ABCD]$"), cd$group))
    idx <- idx[order(x_of[as.character(cd$group[idx])])]
    lines(x_of[as.character(cd$group[idx])], y[idx], col = rep_col[r], lwd = 1.5)
    points(x_of[as.character(cd$group[idx])], y[idx], col = rep_col[r], pch = 16)
  }
  if (g == names(key_id)[1]) legend("topleft", legend = names(rep_col), col = rep_col, lwd = 2, pch = 16, cex = 0.7, bty = "n")
}
invisible(dev.off())
cat("saved pd_genes_paired.png\n")

# ---------------------------------------------------------------- 2. leave one replicate out (full tables)
cat("\n== 2. Leave-one-replicate-out refits\n")
txi <- readRDS("results/tximport/txi.rds")
loo <- list()
for (r in reps) {
  keep_s <- cd$run[cd$replicate != r]
  d <- DESeqDataSetFromTximport(c(lapply(txi[c("abundance", "counts", "length")], function(m) m[rownames(dds), keep_s]),
                                  list(countsFromAbundance = txi$countsFromAbundance)),
                                colData = droplevels(cd[keep_s, ]), design = ~ replicate + group)
  d <- DESeq(d, quiet = TRUE)
  loo[[r]] <- list(`50`  = as.data.frame(results(d, contrast = c("group", "50D", "50A")))[, c("log2FoldChange", "padj")],
                   `500` = as.data.frame(results(d, contrast = c("group", "500D", "500A")))[, c("log2FoldChange", "padj")])
  cat(" done: without", r, "\n")
}

# ---------------------------------------------------------------- 3. high-confidence DEG sets
cat("\n== 3. High-confidence DEGs (main DEG + significant & same sign in all 4 leave-one-out fits + same sign in all 4 replicates)\n")
verdict_rows <- list()
for (sz in c("50", "500")) {
  m <- main_tab[[sz]]; rownames(m) <- m$gene_id
  main_deg <- m$gene_id[!is.na(m$padj) & m$padj < 0.05 & abs(m$FC_shrunk_signed) >= 1.5]
  s <- sign(m[main_deg, "log2FC"])
  loo_ok <- sapply(reps, function(r) { x <- loo[[r]][[sz]][main_deg, ]; !is.na(x$padj) & x$padj < 0.05 & sign(x$log2FoldChange) == s })
  rep_ok <- sign(per_rep[[sz]][main_deg, ]) == s
  hc <- main_deg[rowSums(loo_ok) == 4 & rowSums(rep_ok) == 4]
  cat(sprintf("\n%s nm: main DEGs %d | pass all 4 leave-one-out %d | pass all 4 replicates %d | BOTH (high-confidence) %d (%.0f%%; %d up, %d down)\n",
              sz, length(main_deg), sum(rowSums(loo_ok) == 4), sum(rowSums(rep_ok) == 4), length(hc),
              100 * length(hc) / length(main_deg), sum(m[hc, "log2FC_shrunk"] > 0), sum(m[hc, "log2FC_shrunk"] < 0)))
  h <- m[hc, c("gene_id", "gene_name", "gene_biotype", "baseMean", "log2FC_shrunk", "FC_shrunk_signed", "padj")]
  h$max_padj_leave_one_out <- signif(apply(sapply(reps, function(r) loo[[r]][[sz]][hc, "padj"]), 1, max), 3)
  h <- h[order(-abs(h$log2FC_shrunk)), ]
  write.table(h, file.path(out, sprintf("high_confidence_DEGs_%snm.tsv", sz)), sep = "\t", quote = FALSE, row.names = FALSE)
  for (dir in c("up", "down")) {
    hh <- if (dir == "up") h[h$log2FC_shrunk > 0, ] else h[h$log2FC_shrunk < 0, ]
    cat(sprintf("  top 15 %s (by |log2FC|): %s\n", dir, paste(head(hh$gene_name, 15), collapse = ", ")))
  }

  for (g in names(key_id)) {
    id <- key_id[g]
    lp <- sapply(reps, function(r) loo[[r]][[sz]][id, "padj"])
    ll <- sapply(reps, function(r) loo[[r]][[sz]][id, "log2FoldChange"])
    pr <- per_rep[[sz]][id, ]
    is_deg <- id %in% main_deg
    robust <- all(!is.na(lp) & lp < 0.05) && length(unique(sign(c(ll, m[id, "log2FC"])))) == 1 && length(unique(sign(c(pr, m[id, "log2FC"])))) == 1
    verdict_rows[[paste(g, sz)]] <- data.frame(gene = g, size = sz, baseMean = round(m[id, "baseMean"]),
      log2FC_shrunk = round(m[id, "log2FC_shrunk"], 2), padj = signif(m[id, "padj"], 2), main_DEG_1.5x = is_deg,
      loo_padj = paste(signif(lp, 2), collapse = " / "), per_replicate_l2fc = paste(round(pr, 2), collapse = " / "),
      verdict = ifelse(is.na(m[id, "padj"]) || m[id, "padj"] >= 0.05, "not significant", ifelse(robust, "ROBUST", "fragile")))
  }
}

# ---------------------------------------------------------------- 4. verdict table
cat("\n== 4. Verdict for Parkinson-related genes and Npy\n")
cat("(ROBUST = padj < 0.05 in the main fit AND in all 4 leave-one-out fits, same direction in every fit and every replicate)\n")
vt <- do.call(rbind, verdict_rows); vt <- vt[order(vt$gene, vt$size), ]
print(vt, row.names = FALSE)
write.table(vt, file.path(out, "key_genes_verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink()
cat("done — see", out, "\n")
