# 16_figures.R — final, presentation-quality figures (PNG for the README + PDF vector versions)
#   fig1_pca.png              PCA: dose drives PC1; replicate (experiment) drives PC2/PC3; particle size does not
#   fig2_paired_design.png    DEGs with the paired model vs the authors-like unpaired model vs the published table
#   fig3_volcano.png          dose 100 vs control, both sizes; high-confidence DEGs highlighted, key genes labelled
#   fig4_size_specificity.png 50 nm response vs 500 nm response for every DEG; size-dependent genes (interaction) coloured
#   fig5_signatures.png       GSEA of reactive-astrocyte and cell-type gene sets
#   fig6_replicate_check.png  key genes: change inside each of the 4 experiments (robust vs fragile)
# Run from the project root:  source("scripts/16_figures.R", echo = FALSE)      (~1 min)
# Needs outputs of 08, 09, 12, 14, 15. Packages: DESeq2, ggplot2 (ggrepel optional — nicer labels)

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
out <- "results/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
label_layer <- function(data, ...) if (has_repel) ggrepel::geom_text_repel(data = data, ..., max.overlaps = 50, min.segment.length = 0, seed = 1) else geom_text(data = data, ..., vjust = -0.6)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h)
  cat("saved", name, "\n")
}
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))
size_col <- c(`50 nm` = "#5e3c99", `500 nm` = "#e66101")
rep_col <- c(rep1 = "#1b9e77", rep2 = "#7570b3", rep3 = "#d95f02", rep4 = "#e7298a")

dds <- readRDS("results/deseq2/dds.rds")
cd <- as.data.frame(colData(dds)); cd$group <- as.character(cd$group); cd$replicate <- as.character(cd$replicate)
tab <- function(cn) read.delim(file.path("results/deseq2", paste0(cn, ".tsv")), stringsAsFactors = FALSE)
r50 <- tab("50D_vs_50A"); r500 <- tab("500D_vs_500A"); rint <- tab("response50_minus_response500")
key <- c("Gpnmb", "Npy", "Ccl2", "Cxcl10", "Hspb1", "Lcn2", "C3", "Gfap", "Aqp4", "Clu", "Uchl1", "Pink1", "Abat", "Nes", "Snca")

# ---------------------------------------------------------------- fig1: PCA
vs <- vst(dds, blind = TRUE)
pc <- prcomp(t(assay(vs)[order(-apply(assay(vs), 1, var))[1:2000], ]))
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
# shape = particle series (50 nm or 500 nm); the controls (dose 0) of each series contain no particles
pd <- data.frame(pc$x[, 1:3], size = paste(sub("[ABCD]$", "", cd$group), "nm series"), dose = factor(cd$dose_ugml), replicate = cd$replicate)
p1a <- ggplot(pd, aes(PC1, PC2, colour = dose, shape = size)) + geom_point(size = 2.6) +
  scale_colour_manual(values = c("0" = "grey60", "1" = "#fdb863", "10" = "#e66101", "100" = "#7f2704"), name = "dose (ug/ml)") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), title = "A. PC1 = dose 100")
p1b <- ggplot(pd, aes(PC2, PC3, colour = replicate, shape = size)) + geom_point(size = 2.6) +
  scale_colour_manual(values = rep_col) +
  labs(x = sprintf("PC2 (%.1f%%)", ve[2]), y = sprintf("PC3 (%.1f%%)", ve[3]), title = "B. PC2 / PC3 = experiment (replicate)")
if (requireNamespace("patchwork", quietly = TRUE)) {
  save_fig(patchwork::wrap_plots(p1a, p1b, nrow = 1), "fig1_pca", 10, 4)
} else { save_fig(p1a, "fig1a_pca_dose", 5.2, 4); save_fig(p1b, "fig1b_pca_replicate", 5.2, 4) }

# ---------------------------------------------------------------- fig2: paired vs unpaired design
mc <- read.delim("results/authors_model/model_comparison.tsv", stringsAsFactors = FALSE)
vs_auth <- read.delim("results/validation/validation_summary.tsv", stringsAsFactors = FALSE)
d2 <- rbind(
  data.frame(size = paste(mc$size, "nm"), model = mc$version, DEG = mc$n_DEG)[mc$version %in% c("paired_apeglm", "noRep_apeglm"), ],
  data.frame(size = c("50 nm", "500 nm"), model = "published table",
             DEG = vs_auth$authors_DEG_in_shared[match(c("50D_vs_50A", "500D_vs_500A"), vs_auth$contrast)]))
d2$model <- factor(d2$model, levels = c("published table", "noRep_apeglm", "paired_apeglm"),
                   labels = c("published\n(authors)", "our data,\nunpaired model", "our data,\npaired model"))
p2 <- ggplot(d2, aes(model, DEG, fill = model)) + geom_col(width = 0.7) + geom_text(aes(label = DEG), vjust = -0.3, size = 3.2) +
  facet_wrap(~ size) + scale_fill_manual(values = c("grey60", "#92c5de", "#2166ac"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "DEGs (padj < 0.05, |FC| >= 1.5)",
       title = "Modelling the paired experiments increases power",
       subtitle = "genes present in both our annotation and the published table; dose 100 vs control")
save_fig(p2, "fig2_paired_design", 7, 4)

# ---------------------------------------------------------------- fig3: volcano plots
hc <- list(`50 nm` = read.delim("results/confidence/high_confidence_DEGs_50nm.tsv")$gene_id,
           `500 nm` = read.delim("results/confidence/high_confidence_DEGs_500nm.tsv")$gene_id)
vol <- rbind(cbind(r50, size = "50 nm"), cbind(r500, size = "500 nm"))
vol <- vol[!is.na(vol$padj), ]
vol$y <- pmin(-log10(vol$padj), 40)
vol$class <- ifelse(mapply(function(g, s) g %in% hc[[s]], vol$gene_id, vol$size), "high-confidence DEG",
             ifelse(vol$padj < 0.05 & abs(vol$FC_shrunk_signed) >= 1.5, "DEG (not robust to all checks)", "not DEG"))
vol$class <- factor(vol$class, levels = c("not DEG", "DEG (not robust to all checks)", "high-confidence DEG"))
lab <- vol[vol$gene_name %in% key, ]
p3 <- ggplot(vol[order(vol$class), ], aes(log2FC_shrunk, y, colour = class)) + geom_point(size = 0.5, alpha = 0.7) +
  scale_colour_manual(values = c("grey80", "#fdb863", "#b2182b"), name = NULL) +
  label_layer(lab, aes(label = gene_name), colour = "black", size = 2.6) +
  geom_vline(xintercept = c(-1, 1) * log2(1.5), linetype = 2, colour = "grey50") +
  facet_wrap(~ size) + labs(x = "log2 fold change (shrunken), dose 100 vs control", y = "-log10 padj (capped at 40)") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(override.aes = list(size = 2.5)))
save_fig(p3, "fig3_volcano", 9, 5)

# ---------------------------------------------------------------- fig4: size specificity
m <- merge(r50[, c("gene_id", "gene_name", "log2FC_shrunk", "padj")], r500[, c("gene_id", "log2FC_shrunk", "padj")], by = "gene_id", suffixes = c("_50", "_500"))
m <- merge(m, data.frame(gene_id = rint$gene_id, padj_int = rint$padj), by = "gene_id", all.x = TRUE)
m <- m[(!is.na(m$padj_50) & m$padj_50 < 0.05) | (!is.na(m$padj_500) & m$padj_500 < 0.05), ]
m$class <- ifelse(!is.na(m$padj_int) & m$padj_int < 0.05,
                  ifelse(abs(m$log2FC_shrunk_50) > abs(m$log2FC_shrunk_500), "stronger in 50 nm", "stronger in 500 nm"), "shared response")
m$class <- factor(m$class, levels = c("shared response", "stronger in 50 nm", "stronger in 500 nm"))
p4 <- ggplot(m[order(m$class), ], aes(log2FC_shrunk_50, log2FC_shrunk_500, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 0.7, alpha = 0.7) + scale_colour_manual(values = c("grey70", "#5e3c99", "#e66101"), name = "size x dose interaction\n(padj < 0.05)") +
  label_layer(m[m$gene_name %in% key, ], aes(label = gene_name), colour = "black", size = 2.6) +
  labs(x = "50 nm: log2FC dose 100 vs control", y = "500 nm: log2FC dose 100 vs control",
       title = "Most of the response is shared; a subset depends on particle size",
       subtitle = sprintf("%d genes significant at dose 100 in either size", nrow(m))) +
  coord_equal() + guides(colour = guide_legend(override.aes = list(size = 2.5)))
save_fig(p4, "fig4_size_specificity", 7, 6)

# ---------------------------------------------------------------- fig5: gene-set signatures
sg <- read.delim("results/followup/signature_gsea.tsv", stringsAsFactors = FALSE, check.names = FALSE)
names(sg) <- sub("^X", "", names(sg))
long <- do.call(rbind, lapply(names(sg)[-1], function(cn) {
  v <- sg[[cn]]; ok <- v != "-"
  data.frame(set = sg$set[ok], contrast = cn, NES = as.numeric(sub(" .*", "", v[ok])), padj = as.numeric(gsub(".*\\(|\\)", "", v[ok])))
}))
long$contrast <- factor(long$contrast, levels = c("50D_vs_50A", "500D_vs_500A", "response50_minus_response500"),
                        labels = c("50 nm\ndose 100", "500 nm\ndose 100", "difference\n(50 minus 500)"))
long$set <- factor(long$set, levels = rev(sg$set))
p5 <- ggplot(long, aes(contrast, set, fill = NES, size = -log10(padj))) + geom_point(shape = 21, colour = "grey30") +
  geom_point(data = long[long$padj < 0.05, ], shape = 21, colour = "black", stroke = 1.1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) + scale_size_continuous(range = c(2, 9)) +
  labs(x = NULL, y = NULL, title = "Reactive-astrocyte and cell-type gene sets (GSEA)",
       subtitle = "black outline = padj < 0.05; signatures: Liddelow et al. 2017")
save_fig(p5, "fig5_signatures", 6.5, 4.5)

# ---------------------------------------------------------------- fig6: replicate check for key genes
norm <- counts(dds, normalized = TRUE)
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
g6 <- c("Gpnmb", "Npy", "Uchl1", "Pink1", "Dnajc6", "C3", "Aqp4", "Snca")
rows <- list()
for (g in g6) {
  id <- ann$gene[match(g, ann$gene_name)]; if (is.na(id) || !id %in% rownames(norm)) next
  for (sz in c("50", "500")) for (r in names(rep_col)) {
    d <- norm[id, cd$replicate == r & cd$group == paste0(sz, "D")]; a <- norm[id, cd$replicate == r & cd$group == paste0(sz, "A")]
    rows[[length(rows) + 1]] <- data.frame(gene = g, size = paste(sz, "nm"), replicate = r, l2fc = log2((d + 0.5) / (a + 0.5)))
  }
}
d6 <- do.call(rbind, rows)
est <- rbind(data.frame(gene = r50$gene_name, size = "50 nm", l2fc = r50$log2FC_shrunk, padj = r50$padj),
             data.frame(gene = r500$gene_name, size = "500 nm", l2fc = r500$log2FC_shrunk, padj = r500$padj))
est <- est[est$gene %in% g6, ]
d6$gene <- factor(d6$gene, levels = rev(g6)); est$gene <- factor(est$gene, levels = rev(g6))
p6 <- ggplot(d6, aes(l2fc, gene)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(colour = replicate), size = 2.2, position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_point(data = est, shape = 23, size = 3, fill = "black") +
  scale_colour_manual(values = rep_col) + facet_wrap(~ size) +
  labs(x = "log2 fold change, dose 100 vs control", y = NULL, title = "Does every experiment agree?",
       subtitle = "dots: inside each replicate; diamond: model estimate (shrunken). Snca is driven by rep3 alone")
save_fig(p6, "fig6_replicate_check", 8, 4.5)
cat("done — see", out, "\n")
