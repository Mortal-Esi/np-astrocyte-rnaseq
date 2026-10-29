# 14_followup.R — three follow-ups to the enrichment results (13_enrichment.R)
#   A. Is the KEGG "Parkinson disease" signal (and AD / HD / prion / ALS) more than oxidative phosphorylation,
#      ribosome and proteasome genes? Test: remove those "shared machinery" genes from each disease pathway and
#      re-run GSEA on the RESIDUAL pathway. Also: where do the classic PD genes (Snca, Lrrk2, Pink1 ...) sit?
#   B. Reactive astrocyte signatures (Liddelow et al. 2017, Nature 541:481): pan-reactive, "A1", "A2"
#      + cell-type / cell-state marker sets (homeostatic astrocyte, ependymal/ciliated, proliferation,
#      oligodendrocyte lineage, neuronal). GSEA per contrast + a per-sample score, one line per replicate.
#   C. Ciliated / ependymal markers gene by gene: do they all drop together (= fewer ciliated cells)
#      or only some (= cells switching off part of a program)?
# Run from the project root:  source("scripts/14_followup.R", echo = FALSE)      (~3–5 min; part A needs KEGG = internet/VPN)
# Needs: results/deseq2/dds.rds and *.tsv (08). Packages: DESeq2, clusterProfiler, org.Rn.eg.db, ggplot2
# Outputs: results/followup/ (followup_log.txt, *.tsv, signature_scores_paired.png, signature_genes_heatmap.png)

suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler); library(org.Rn.eg.db); library(ggplot2)
})
out <- "results/followup"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
sink(file.path(out, "followup_log.txt"), split = TRUE)
set.seed(2024)

contrasts <- c("50D_vs_50A", "500D_vs_500A", "response50_minus_response500")
all_contrasts <- c("50B_vs_50A", "50C_vs_50A", "50D_vs_50A", "500B_vs_500A", "500C_vs_500A", "500D_vs_500A")
read_tab <- function(cn) {
  d <- read.delim(file.path("results/deseq2", paste0(cn, ".tsv")), stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj), ]; d$stat <- d$log2FC / d$lfcSE; d
}
res <- setNames(lapply(contrasts, read_tab), contrasts)
rank_by <- function(d, id) {                      # ranked list (Wald statistic), one value per ID
  d <- d[!is.na(d[[id]]) & is.finite(d$stat), ]
  d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d[[id]]), ]
  sort(setNames(d$stat, d[[id]]), decreasing = TRUE)
}
fmt_gsea <- function(g) {                         # one-line "NES (padj)" per set
  df <- as.data.frame(g)
  setNames(sprintf("%+.2f (%.1e)", df$NES, df$p.adjust), df$ID)
}

# ================================================================ A. neurodegeneration pathways vs shared machinery
cat("== A. KEGG neurodegeneration pathways: real disease signal or shared mitochondrial / ribosome / proteasome genes?\n")
ens_keys <- keys(org.Rn.eg.db, keytype = "ENSEMBL"); sym_keys <- keys(org.Rn.eg.db, keytype = "SYMBOL")
to_entrez <- function(ids, names) {               # same rule as 13_enrichment.R
  e <- rep(NA_character_, length(ids)); ok <- ids %in% ens_keys
  if (any(ok)) e[ok] <- suppressMessages(mapIds(org.Rn.eg.db, keys = ids[ok], keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first"))
  ok2 <- is.na(e) & names %in% sym_keys
  if (any(ok2)) e[ok2] <- suppressMessages(mapIds(org.Rn.eg.db, keys = names[ok2], keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
  e
}
for (cn in contrasts) res[[cn]]$entrez <- to_entrez(res[[cn]]$gene_id, res[[cn]]$gene_name)
sym_of <- function(eg) unname(suppressMessages(mapIds(org.Rn.eg.db, keys = eg, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")))

disease <- c(Parkinson = "rno05012", Alzheimer = "rno05010", Huntington = "rno05016", Prion = "rno05020",
             ALS = "rno05014", Neurodegeneration_multi = "rno05022")
machinery <- c(OXPHOS = "rno00190", Ribosome = "rno03010", Proteasome = "rno03050")
pd_genes <- c("Snca", "Lrrk2", "Prkn", "Pink1", "Park7", "Uchl1", "Gba1", "Gba", "Vps35", "Atp13a2", "Fbxo7", "Synj1", "Dnajc6", "Gpnmb")

kegg_ok <- TRUE
for (cn in contrasts) {
  gl <- rank_by(res[[cn]], "entrez")
  gk <- tryCatch(gseKEGG(gl, organism = "rno", minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE),
                 error = function(e) { cat("KEGG unreachable (", conditionMessage(e), ") — part A skipped; retry with VPN.\n"); NULL })
  if (is.null(gk)) { kegg_ok <- FALSE; break }
  gs <- gk@geneSets
  shared <- unique(unlist(gs[machinery]))
  t2g <- do.call(rbind, c(
    lapply(names(machinery), function(n) data.frame(term = n, gene = gs[[machinery[n]]])),
    lapply(names(disease), function(n) data.frame(term = paste0(n, "_full"), gene = gs[[disease[n]]])),
    lapply(names(disease), function(n) data.frame(term = paste0(n, "_residual"), gene = setdiff(gs[[disease[n]]], shared)))))
  g <- GSEA(gl, TERM2GENE = t2g, minGSSize = 10, maxGSSize = 2000, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
  df <- as.data.frame(g)
  cat(sprintf("\n-- %s\n", cn))
  tab <- data.frame(pathway = names(disease),
    genes_in_list = sapply(names(disease), function(n) sum(gs[[disease[n]]] %in% names(gl))),
    pct_shared_machinery = sapply(names(disease), function(n) round(100 * mean(gs[[disease[n]]][gs[[disease[n]]] %in% names(gl)] %in% shared))),
    full = sapply(names(disease), function(n) { r <- df[df$ID == paste0(n, "_full"), ]; if (nrow(r)) sprintf("%+.2f (%.1e)", r$NES, r$p.adjust) else "-" }),
    residual = sapply(names(disease), function(n) { r <- df[df$ID == paste0(n, "_residual"), ]; if (nrow(r)) sprintf("%+.2f (%.1e)", r$NES, r$p.adjust) else "-" }),
    residual_leading_edge_top = sapply(names(disease), function(n) { r <- df[df$ID == paste0(n, "_residual"), ]
      if (!nrow(r)) return("-"); paste(head(sym_of(strsplit(r$core_enrichment, "/")[[1]]), 8), collapse = ",") }))
  print(tab, row.names = FALSE)
  cat("machinery sets: "); m <- df[df$ID %in% names(machinery), ]; cat(paste0(m$ID, " ", sprintf("%+.2f (%.1e)", m$NES, m$p.adjust), collapse = " | "), "\n")
  write.table(df[, setdiff(names(df), "core_enrichment")], file.path(out, sprintf("neurodegeneration_residual_%s.tsv", cn)), sep = "\t", quote = FALSE, row.names = FALSE)

  d <- res[[cn]]; pd <- d[d$gene_name %in% pd_genes, ]
  le_pd <- if (any(df$ID == "Parkinson_full")) strsplit(df$core_enrichment[df$ID == "Parkinson_full"], "/")[[1]] else character(0)
  pd$rank_pct <- round(100 * (match(pd$entrez, names(gl)) / length(gl)), 1)       # 0 = top (most up), 100 = bottom (most down)
  pd$in_PD_leading_edge <- pd$entrez %in% le_pd
  cat("classic PD genes (rank_pct: 0 = most up, 100 = most down):\n")
  print(data.frame(gene = pd$gene_name, baseMean = round(pd$baseMean), log2FC = round(pd$log2FC, 2), padj = signif(pd$padj, 2),
                   rank_pct = pd$rank_pct, in_PD_leading_edge = pd$in_PD_leading_edge)[order(pd$rank_pct), ], row.names = FALSE)
}

# ================================================================ B. reactive astrocyte signatures + cell-type / state markers
cat("\n\n== B. Gene-set scores: reactive astrocyte signatures (Liddelow 2017) and cell-type / state markers\n")
sets <- list(
  PAN_reactive = c("Lcn2", "Steap4", "S1pr3", "Timp1", "Hspb1", "Cxcl10", "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap"),
  A1           = c("H2-T23", "Serping1", "H2-D1", "Ggta1", "Iigp1", "Gbp2", "Fbln5", "Ugt1a1", "Fkbp5", "Psmb8", "Srgn", "Amigo2"),
  A2           = c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14"),
  homeostatic_astro = c("Aqp4", "Slc1a2", "Slc1a3", "Glul", "Kcnj10", "Aldh1l1", "Gja1", "Sparcl1", "Atp1a2", "Slc4a4"),
  ependymal_cilia   = c("Foxj1", "Dnah5", "Dnah9", "Rsph1", "Rsph4a", "Ccdc153", "Tmem212", "Dynlrb2", "Spag6", "Sntn"),
  proliferation     = c("Mki67", "Top2a", "Ccnb1", "Cdk1", "Birc5", "Aurkb", "Cenpf", "Pcna"),
  oligo_lineage     = c("Pdgfra", "Olig2", "Cspg4", "Sox10", "Plp1", "Mbp", "Mog"),
  neuronal          = c("Rbfox3", "Tubb3", "Syt1", "Snap25", "Stmn2", "Dcx"))
dds <- readRDS("results/deseq2/dds.rds")
cd <- as.data.frame(colData(dds)); cd$group <- as.character(cd$group); cd$replicate <- as.character(cd$replicate)
norm <- counts(dds, normalized = TRUE)
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
id_of <- function(sym) { id <- ann$gene[match(sym, ann$gene_name)]; id[!is.na(id) & id %in% rownames(norm)] }

cat("genes found in our data (expressed = mean normalised count >= 10):\n")
for (s in names(sets)) {
  found <- sets[[s]][sets[[s]] %in% ann$gene_name[ann$gene %in% rownames(norm)]]
  expr <- found[rowMeans(norm[id_of(found), , drop = FALSE]) >= 10]
  missing <- setdiff(sets[[s]], found)
  cat(sprintf("  %-18s %2d/%2d found, %2d expressed%s\n", s, length(found), length(sets[[s]]), length(expr),
              if (length(missing)) paste0(" | not found (no rat gene with this symbol): ", paste(missing, collapse = ",")) else ""))
}

t2g_sym <- do.call(rbind, lapply(names(sets), function(s) data.frame(term = s, gene = sets[[s]])))
gsea_sets <- sapply(contrasts, function(cn) {
  g <- GSEA(rank_by(res[[cn]], "gene_name"), TERM2GENE = t2g_sym, minGSSize = 4, maxGSSize = 500, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
  f <- fmt_gsea(g); f[names(sets)]
})
rownames(gsea_sets) <- names(sets); gsea_sets[is.na(gsea_sets)] <- "-"
cat("\nGSEA per gene set — NES (padj); positive = shifted up; interaction + = response bigger in 50 nm:\n")
print(noquote(gsea_sets))
write.table(data.frame(set = rownames(gsea_sets), gsea_sets), file.path(out, "signature_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# per-sample score: for each expressed gene, log2 expression centred on the mean of ALL 8 control samples and scaled by its SD;
# score = mean over the set's genes. 0 = like an average control.
L <- log2(norm + 1); ctrl <- cd$group %in% c("50A", "500A")
score <- sapply(names(sets), function(s) {
  ids <- id_of(sets[[s]]); ids <- ids[rowMeans(norm[ids, , drop = FALSE]) >= 10]
  if (length(ids) < 3) return(rep(NA, ncol(L)))
  z <- (L[ids, , drop = FALSE] - rowMeans(L[ids, ctrl, drop = FALSE])) / apply(L[ids, , drop = FALSE], 1, sd)
  colMeans(z)
})
rownames(score) <- colnames(L)
write.table(data.frame(sample = rownames(score), group = cd$group, replicate = cd$replicate, round(score, 3)),
            file.path(out, "signature_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

x_of <- c("50A" = 1, "50B" = 2, "50C" = 3, "50D" = 4, "500A" = 6, "500B" = 7, "500C" = 8, "500D" = 9)
rep_col <- c(rep1 = "#1b9e77", rep2 = "#7570b3", rep3 = "#d95f02", rep4 = "#e7298a")
png(file.path(out, "signature_scores_paired.png"), width = 2400, height = 1300, res = 170)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (s in names(sets)) {
  y <- score[, s]; if (all(is.na(y))) { plot.new(); title(paste(s, "(too few genes)")); next }
  plot(NA, xlim = c(0.5, 9.5), ylim = range(y, na.rm = TRUE), xaxt = "n", xlab = "", ylab = "score (mean z vs controls)", main = s)
  axis(1, at = x_of, labels = names(x_of), las = 2, cex.axis = 0.75); abline(v = 5, col = "grey85"); abline(h = 0, col = "grey85", lty = 2)
  for (r in names(rep_col)) for (sz in c("50", "500")) {
    idx <- which(cd$replicate == r & grepl(paste0("^", sz, "[ABCD]$"), cd$group))
    idx <- idx[order(x_of[cd$group[idx]])]
    lines(x_of[cd$group[idx]], y[idx], col = rep_col[r], lwd = 1.5); points(x_of[cd$group[idx]], y[idx], col = rep_col[r], pch = 16)
  }
  if (s == names(sets)[1]) legend("topleft", legend = names(rep_col), col = rep_col, lwd = 2, pch = 16, cex = 0.7, bty = "n")
}
invisible(dev.off())
cat("\nsaved signature_scores_paired.png\n")

# heatmap of shrunken log2FC for every signature gene, all doses
lfc <- do.call(rbind, lapply(all_contrasts, function(cn) {
  d <- read.delim(file.path("results/deseq2", paste0(cn, ".tsv")), stringsAsFactors = FALSE)
  do.call(rbind, lapply(names(sets), function(s) {
    h <- d[d$gene_name %in% sets[[s]], ]
    if (!nrow(h)) return(NULL)
    data.frame(set = s, gene = h$gene_name, contrast = sub("_vs_.*", "", cn), log2FC = h$log2FC_shrunk, sig = !is.na(h$padj) & h$padj < 0.05)
  }))
}))
lfc$contrast <- factor(lfc$contrast, levels = sub("_vs_.*", "", all_contrasts))
lfc$set <- factor(lfc$set, levels = names(sets))
lim <- max(abs(lfc$log2FC), na.rm = TRUE)
p <- ggplot(lfc, aes(contrast, gene, fill = log2FC)) + geom_tile(colour = "white") +
  geom_point(data = lfc[lfc$sig, ], shape = 16, size = 0.8) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim)) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 8) +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Shrunken log2FC vs control (dot = padj < 0.05)")
ggsave(file.path(out, "signature_genes_heatmap.png"), p, width = 6, height = 13, dpi = 200)
cat("saved signature_genes_heatmap.png\n")

# ================================================================ C. ependymal / cilia genes one by one
cat("\n\n== C. Ependymal / cilia markers at dose 100 (does everything drop together?)\n")
cil <- sets$ependymal_cilia
per_rep <- function(g, sz) {
  id <- id_of(g); if (!length(id)) return(rep(NA, 4))
  sapply(names(rep_col), function(r) log2((norm[id, rownames(cd)[cd$replicate == r & cd$group == paste0(sz, "D")]] + 0.5) /
                                           (norm[id, rownames(cd)[cd$replicate == r & cd$group == paste0(sz, "A")]] + 0.5)))
}
tc <- do.call(rbind, lapply(cil, function(g) {
  a <- res[["50D_vs_50A"]]; b <- res[["500D_vs_500A"]]
  ra <- a[a$gene_name == g, ]; rb <- b[b$gene_name == g, ]
  if (!nrow(ra) || !nrow(rb)) return(NULL)
  data.frame(gene = g, baseMean = round(rb$baseMean[1]),
             l2fc_50D = round(ra$log2FC_shrunk[1], 2), padj_50D = signif(ra$padj[1], 2),
             l2fc_500D = round(rb$log2FC_shrunk[1], 2), padj_500D = signif(rb$padj[1], 2),
             reps_500D = paste(round(per_rep(g, "500"), 1), collapse = " "))
}))
if (is.null(tc)) cat("  (none of the markers were tested)\n") else print(tc, row.names = FALSE)
cat("(reps_500D = log2(500D / 500A) inside rep1 rep2 rep3 rep4)\n")
sink()
cat("done — see", out, "\n")
