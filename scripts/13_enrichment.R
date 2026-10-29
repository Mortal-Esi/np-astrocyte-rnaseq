# 13_enrichment.R — which biological processes do the DEGs point to?
#   Contrasts: 50D_vs_50A, 500D_vs_500A, response50_minus_response500 (interaction: + = stronger/more up in 50 nm)
#   1. ORA  (over-representation): are GO terms over-represented among up / down DEGs?
#           universe = genes DESeq2 actually tested (padj not NA), NOT the whole genome
#   2. GSEA (gene set enrichment): rank ALL tested genes by the Wald statistic (log2FC / lfcSE) and ask whether
#           a GO term's genes sit at the top or bottom of the ranking (no DEG cut-off needed)
#   3. KEGG (ORA + GSEA) — needs internet access to rest.kegg.jp; skipped with a message if unreachable
#   4. Theme check: best GSEA term for a few pre-specified themes (inflammation, interferon, complement, ...)
# Run from the project root:  source("scripts/13_enrichment.R", echo = FALSE)      (~5–15 min)
# Needs: results/deseq2/*.tsv (08). Packages: clusterProfiler, org.Rn.eg.db, enrichplot, ggplot2
# Outputs: results/enrichment/ (enrichment_log.txt, ora_*.tsv, gsea_*.tsv, kegg_*.tsv, *.png)

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Rn.eg.db); library(enrichplot); library(ggplot2)
})
out <- "results/enrichment"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
sink(file.path(out, "enrichment_log.txt"), split = TRUE)
set.seed(2024)

contrasts <- c("50D_vs_50A", "500D_vs_500A", "response50_minus_response500")
padj_cut <- 0.05; lfc_cut <- log2(1.5)                  # same DEG criterion as the authors (|FC| >= 1.5)
cat("versions: clusterProfiler", as.character(packageVersion("clusterProfiler")),
    "| org.Rn.eg.db", as.character(packageVersion("org.Rn.eg.db")), "\n")

# ---------------------------------------------------------------- Ensembl → Entrez (GO/KEGG use Entrez IDs)
# org.Rn.eg.db is built from NCBI; its Ensembl links can lag behind Ensembl 116 (GRCr8),
# so genes not found by Ensembl ID are tried by gene symbol.
ens_keys <- keys(org.Rn.eg.db, keytype = "ENSEMBL"); sym_keys <- keys(org.Rn.eg.db, keytype = "SYMBOL")
to_entrez <- function(ids, names) {
  e <- rep(NA_character_, length(ids))
  ok <- ids %in% ens_keys
  if (any(ok)) e[ok] <- mapIds(org.Rn.eg.db, keys = ids[ok], keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")
  ok2 <- is.na(e) & names %in% sym_keys
  if (any(ok2)) e[ok2] <- mapIds(org.Rn.eg.db, keys = names[ok2], keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  e
}

res <- list()
for (cn in contrasts) {
  d <- read.delim(file.path("results/deseq2", paste0(cn, ".tsv")), stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj), ]                               # genes DESeq2 really tested = the universe
  d$entrez <- to_entrez(d$gene_id, d$gene_name)
  d$stat <- d$log2FC / d$lfcSE
  lfc <- d$log2FC_shrunk                                 # interaction table: shrunk column = unshrunken value
  d$dir <- ifelse(d$padj < padj_cut & lfc >= lfc_cut, "up", ifelse(d$padj < padj_cut & lfc <= -lfc_cut, "down", "ns"))
  res[[cn]] <- d
}

cat("\n== 0. Gene ID mapping and DEG numbers\n")
for (cn in contrasts) {
  d <- res[[cn]]
  cat(sprintf("%-30s tested %5d | mapped to Entrez %5d (%.1f%%) | up %4d (mapped %4d) | down %4d (mapped %4d)\n", cn,
              nrow(d), sum(!is.na(d$entrez)), 100 * mean(!is.na(d$entrez)),
              sum(d$dir == "up"), sum(d$dir == "up" & !is.na(d$entrez)),
              sum(d$dir == "down"), sum(d$dir == "down" & !is.na(d$entrez))))
}
label_of <- function(cn, dir) {
  if (cn == "response50_minus_response500") return(if (dir == "up") "interaction: stronger in 50 nm" else "interaction: stronger in 500 nm")
  paste(sub("_vs_.*", "", cn), dir)
}

show_terms <- function(x, n = 12) {                     # compact printing of an enrichResult / gseaResult / data.frame
  df <- if (is.null(x)) data.frame() else as.data.frame(x)
  if (nrow(df) == 0) { cat("   (no significant terms)\n"); return(invisible()) }
  if ("NES" %in% names(df)) {
    df <- df[order(df$p.adjust), ]
    print(data.frame(term = substr(df$Description, 1, 60), size = df$setSize, NES = round(df$NES, 2),
                     padj = signif(df$p.adjust, 2))[seq_len(min(n, nrow(df))), ], row.names = FALSE)
  } else {
    fe <- sapply(strsplit(df$GeneRatio, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2])) /
          sapply(strsplit(df$BgRatio, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
    print(data.frame(term = substr(df$Description, 1, 60), genes = df$Count, fold_enrich = round(fe, 1),
                     padj = signif(df$p.adjust, 2))[seq_len(min(n, nrow(df))), ], row.names = FALSE)
  }
}
simplify_safe <- function(x) tryCatch(simplify(x, cutoff = 0.7, by = "p.adjust", select_fun = min),
                                      error = function(e) { cat("   (simplify failed:", conditionMessage(e), "— showing unsimplified)\n"); x })

# ---------------------------------------------------------------- 1. ORA (GO biological process)
cat("\n== 1. ORA — GO biological process (redundant terms merged with simplify, similarity cutoff 0.7)\n")
ora <- list()
for (cn in contrasts) for (dir in c("up", "down")) {
  d <- res[[cn]]; lab <- label_of(cn, dir)
  univ <- unique(na.omit(d$entrez)); genes <- unique(na.omit(d$entrez[d$dir == dir]))
  cat(sprintf("\n-- %s (%d genes)\n", lab, length(genes)))
  if (length(genes) < 10) { cat("   (fewer than 10 genes — skipped)\n"); next }
  eg <- enrichGO(genes, universe = univ, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "BP",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 15, maxGSSize = 500, readable = TRUE)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) { cat("   (no significant terms)\n"); next }
  write.table(as.data.frame(eg), file.path(out, sprintf("ora_GO_%s_%s.tsv", cn, dir)), sep = "\t", quote = FALSE, row.names = FALSE)
  egs <- simplify_safe(eg)
  cat(sprintf("   %d significant terms, %d after simplify\n", nrow(as.data.frame(eg)), nrow(as.data.frame(egs))))
  show_terms(egs); ora[[lab]] <- egs
}
if (length(ora)) {
  long <- do.call(rbind, lapply(names(ora), function(l) { df <- as.data.frame(ora[[l]]); df <- df[order(df$p.adjust), ][seq_len(min(6, nrow(df))), ]
    data.frame(set = l, term = substr(df$Description, 1, 55), padj = df$p.adjust, genes = df$Count) }))
  long$set <- factor(long$set, levels = names(ora))
  p <- ggplot(long, aes(set, term, size = genes, colour = -log10(padj))) + geom_point() +
       scale_colour_viridis_c() + theme_bw(base_size = 9) + labs(x = NULL, y = NULL, title = "ORA, GO BP: top 6 terms per gene set") +
       theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(out, "ora_GO_overview.png"), p, width = 9, height = max(4, 0.22 * length(unique(long$term)) + 2), dpi = 200)
}

# ---------------------------------------------------------------- 2. GSEA (GO biological process)
cat("\n\n== 2. GSEA — GO biological process, genes ranked by Wald statistic (positive NES = shifted up)\n")
rank_list <- function(d) {
  d <- d[!is.na(d$entrez) & is.finite(d$stat), ]
  d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$entrez), ]          # one value per Entrez gene
  sort(setNames(d$stat, d$entrez), decreasing = TRUE)
}
gsea <- list()
for (cn in contrasts) {
  gl <- rank_list(res[[cn]])
  cat(sprintf("\n-- %s (%d ranked genes)\n", cn, length(gl)))
  gg <- tryCatch(gseGO(gl, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "BP", minGSSize = 15, maxGSSize = 500,
                       pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE),
                 error = function(e) { cat("   gseGO failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(gg) || nrow(as.data.frame(gg)) == 0) { cat("   (no significant terms)\n"); next }
  gg <- setReadable(gg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(gg), file.path(out, sprintf("gsea_GO_%s.tsv", cn)), sep = "\t", quote = FALSE, row.names = FALSE)
  ggs <- simplify_safe(gg); gsea[[cn]] <- gg
  df <- as.data.frame(ggs)
  cat(sprintf("   %d significant terms (%d up, %d down); %d after simplify\n", nrow(as.data.frame(gg)),
              sum(as.data.frame(gg)$NES > 0), sum(as.data.frame(gg)$NES < 0), nrow(df)))
  cat("   top shifted UP:\n");   show_terms(df[df$NES > 0, ], 10)
  cat("   top shifted DOWN:\n"); show_terms(df[df$NES < 0, ], 10)
  p <- dotplot(ggs, showCategory = 10, split = ".sign", font.size = 8) + facet_grid(. ~ .sign) + ggtitle(paste("GSEA GO BP:", cn))
  ggsave(file.path(out, sprintf("gsea_GO_%s.png", cn)), p, width = 10, height = 7, dpi = 200)
}

# ---------------------------------------------------------------- 3. KEGG (needs internet)
cat("\n\n== 3. KEGG (rat, organism = 'rno')\n")
kegg_ok <- TRUE
for (cn in contrasts) {
  if (!kegg_ok) break
  d <- res[[cn]]; univ <- unique(na.omit(d$entrez))
  for (dir in c("up", "down")) {
    genes <- unique(na.omit(d$entrez[d$dir == dir])); if (length(genes) < 10) next
    ek <- tryCatch(enrichKEGG(genes, organism = "rno", universe = univ, pvalueCutoff = 0.05, minGSSize = 15, maxGSSize = 500),
                   error = function(e) { cat("KEGG unreachable (", conditionMessage(e), ") — skipped. Retry with VPN if needed.\n"); NULL })
    if (is.null(ek)) { kegg_ok <- FALSE; break }
    ek <- setReadable(ek, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    ek@result$Description <- sub(" - Rattus norvegicus \\(rat\\)$", "", ek@result$Description)
    cat(sprintf("\n-- KEGG ORA %s\n", label_of(cn, dir))); show_terms(ek, 10)
    write.table(as.data.frame(ek), file.path(out, sprintf("kegg_ora_%s_%s.tsv", cn, dir)), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!kegg_ok) break
  gk <- tryCatch(gseKEGG(rank_list(d), organism = "rno", minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE),
                 error = function(e) { cat("gseKEGG failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(gk) && nrow(as.data.frame(gk))) {
    gk@result$Description <- sub(" - Rattus norvegicus \\(rat\\)$", "", gk@result$Description)
    cat(sprintf("\n-- KEGG GSEA %s\n", cn)); show_terms(gk, 15)
    write.table(as.data.frame(gk), file.path(out, sprintf("kegg_gsea_%s.tsv", cn)), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# ---------------------------------------------------------------- 4. pre-specified themes (from GSEA, unsimplified)
cat("\n\n== 4. Theme check: best GSEA GO term per theme (NES / padj; '-' = no significant term)\n")
themes <- c(inflammation = "inflammatory response", interferon_virus = "interferon|virus", complement = "complement",
            chemotaxis = "chemotaxis", lysosome_autophagy = "lysosom|autophag", cholesterol = "cholesterol|sterol",
            oxidative_stress = "oxidative stress|reactive oxygen", cell_cycle = "cell cycle|mitotic|DNA replication",
            extracellular_matrix = "extracellular matrix", synapse = "synap", glutamate_GABA = "glutamate|gamma-aminobutyric",
            neuron_differentiation = "neuron differentiation|axonogenesis|neurogenesis")
tab <- sapply(names(gsea), function(cn) {
  df <- as.data.frame(gsea[[cn]])
  sapply(themes, function(pat) {
    h <- df[grepl(pat, df$Description, ignore.case = TRUE), ]
    if (!nrow(h)) return("-")
    h <- h[which.min(h$p.adjust), ]
    sprintf("%s (%+.2f / %.1e)", substr(h$Description, 1, 38), h$NES, h$p.adjust)
  })
})
if (length(gsea)) print(noquote(tab))
sink()
cat("done — see", out, "\n")
