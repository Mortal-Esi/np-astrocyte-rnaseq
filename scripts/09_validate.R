# 09_validate.R — compare our DESeq2 results with the authors' published DE tables (GEO supplementary)
#   - how many of the authors' genes exist in our analysis (different assembly / filtering)
#   - correlation of fold changes for shared genes
#   - overlap of DEG lists with the authors' own criterion (padj < 0.05 & |FC| >= 1.5)
#   - a table of key genes side by side
# Run from the project root:  source("scripts/09_validate.R", echo = FALSE)
# Outputs: results/validation/*.png, *.tsv, validation_log.txt

out <- "results/validation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out, "validation_log.txt"), split = TRUE)

to_log2 <- function(fc) sign(fc) * log2(abs(fc))          # authors' signed FC → log2 FC
read_authors <- function(f) {
  a <- read.delim(gzfile(f), stringsAsFactors = FALSE)
  a$EnsemblGene <- sub("\\.[0-9]+$", "", a$EnsemblGene)
  a
}
auth <- list(`50`  = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_50.tsv.gz"),
             `500` = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_500.tsv.gz"))
ours <- function(name) read.delim(file.path("results/deseq2", paste0(name, ".tsv")), stringsAsFactors = FALSE)

key_genes <- c("Gpnmb", "Lgals3", "Hspb1", "Atf3", "Ccl2", "Ccl7", "Cxcl10", "C3", "Lcn2", "Gfap", "Vim",
               "Nes", "Aqp4", "Lix1", "Abat", "Slc7a11", "Slc38a3", "Psat1", "Gstm1", "Clu", "Prdx6", "Npy", "Snca")
summary_rows <- list(); key_tab <- list()

for (size in c("50", "500")) {
  for (dose in c("B", "C", "D")) {
    name <- sprintf("%s%s_vs_%sA", size, dose, size)
    a <- auth[[size]]
    fc_col <- sprintf("FoldChange_%s_vs_A", dose); p_col <- sprintf("padj_%s_vs_A", dose)
    o <- ours(name)
    m <- merge(data.frame(gene_id = a$EnsemblGene, gene_name_a = a$Gene_name,
                          l2_a = to_log2(a[[fc_col]]), padj_a = a[[p_col]]),
               o[, c("gene_id", "gene_name", "baseMean", "log2FC_shrunk", "FC_shrunk_signed", "padj")],
               by = "gene_id")
    deg_a <- m$padj_a < 0.05 & abs(m$l2_a) >= log2(1.5)
    deg_o <- !is.na(m$padj) & m$padj < 0.05 & abs(m$FC_shrunk_signed) >= 1.5
    n_a_all <- sum(a[[p_col]] < 0.05 & abs(a[[fc_col]]) >= 1.5)
    row <- data.frame(contrast = name,
      authors_genes = nrow(a), our_genes = nrow(o), shared_genes = nrow(m),
      pct_authors_found = round(100 * nrow(m) / nrow(a), 1),
      cor_pearson = round(cor(m$l2_a, m$log2FC_shrunk, use = "complete.obs"), 3),
      cor_spearman = round(cor(m$l2_a, m$log2FC_shrunk, use = "complete.obs", method = "spearman"), 3),
      authors_DEG = n_a_all, authors_DEG_in_shared = sum(deg_a), our_DEG_in_shared = sum(deg_o),
      both = sum(deg_a & deg_o), only_authors = sum(deg_a & !deg_o), only_ours = sum(!deg_a & deg_o),
      pct_authors_DEG_recovered = if (sum(deg_a)) round(100 * sum(deg_a & deg_o) / sum(deg_a), 1) else NA,
      same_direction_in_both = if (sum(deg_a & deg_o)) round(100 * mean(sign(m$l2_a[deg_a & deg_o]) == sign(m$log2FC_shrunk[deg_a & deg_o])), 1) else NA)
    summary_rows[[name]] <- row

    if (dose == "D") {
      png(file.path(out, paste0("fc_scatter_", name, ".png")), width = 1500, height = 1500, res = 200)
      col <- ifelse(deg_a & deg_o, "#7b3294", ifelse(deg_a, "#e66101", ifelse(deg_o, "#5e3c99", "grey75")))
      col[!(deg_a | deg_o)] <- "grey80"
      plot(m$l2_a, m$log2FC_shrunk, pch = 16, cex = 0.35, col = col,
           xlab = "authors: log2 fold change", ylab = "this analysis: shrunken log2 fold change",
           main = sprintf("%s  (Pearson r = %.2f)", name, row$cor_pearson))
      abline(0, 1, lty = 2, col = "grey30"); abline(h = 0, v = 0, col = "grey60")
      legend("topleft", bty = "n", cex = 0.8, pch = 16,
             col = c("#7b3294", "#e66101", "#5e3c99", "grey80"),
             legend = c("DEG in both", "DEG authors only", "DEG ours only", "neither"))
      invisible(dev.off())
      k <- m[match(key_genes, m$gene_name), c("gene_name", "l2_a", "padj_a", "log2FC_shrunk", "padj")]
      k$gene_name <- key_genes
      names(k) <- c("gene", paste0(size, "nm_authors_log2FC"), paste0(size, "nm_authors_padj"),
                    paste0(size, "nm_ours_log2FC"), paste0(size, "nm_ours_padj"))
      key_tab[[size]] <- k
    }
  }
}

summ <- do.call(rbind, summary_rows)
cat("== 1. Agreement with the authors (criterion: padj < 0.05 & |FC| >= 1.5)\n")
tt <- t(summ[, -1]); colnames(tt) <- summ$contrast
print(noquote(format(tt, drop0trailing = TRUE)))
write.table(summ, file.path(out, "validation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n== 2. Key genes at dose 100 (log2 fold change; padj)\n")
kt <- merge(key_tab[["50"]], key_tab[["500"]], by = "gene", sort = FALSE)
kt <- kt[match(key_genes, kt$gene), ]
fmt <- kt; for (j in 2:ncol(fmt)) fmt[[j]] <- if (grepl("padj", names(fmt)[j])) signif(fmt[[j]], 2) else round(fmt[[j]], 2)
print(fmt, row.names = FALSE)
write.table(kt, file.path(out, "key_genes_dose100.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink()
cat("done — see", out, "\n")
