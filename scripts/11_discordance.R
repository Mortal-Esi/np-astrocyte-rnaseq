# 11_discordance.R — two follow-up checks after 09 (validation) and 10 (robustness)
#   A. Where do we and the authors disagree at dose 100, and why?
#      (ID vs name matching, expression level, number of transcripts, authors' FC exactly 1 / padj exactly 1,
#       and whether every replicate supports OUR direction)
#   B. Is sample 500C of rep3 (SRR28006125) a genuine 500 nm sample with a strong response, or could it be mislabelled?
#      Each sample's response is compared with the 50-nm-specific and 500-nm-specific response pattern,
#      learned ONLY from the other three replicates (so the tested sample never helps define the pattern).
# Run from the project root:  source("scripts/11_discordance.R", echo = FALSE)
# Needs: results/deseq2/dds.rds and *_vs_*.tsv (08), data/authors/*.tsv.gz, results/tximport/*.
# Outputs: results/discordance/ (discordance_log.txt, discordant_genes_<size>.tsv, sample_signature_scores.tsv, signature_map.png)
# v2: section B uses a two-pattern regression (v1 correlation score failed its positive control)

suppressPackageStartupMessages(library(DESeq2))
out <- "results/discordance"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
sink(file.path(out, "discordance_log.txt"), split = TRUE)

dds <- readRDS("results/deseq2/dds.rds")
cd  <- as.data.frame(colData(dds))
cd$group <- as.character(cd$group); cd$replicate <- as.character(cd$replicate)
norm <- counts(dds, normalized = TRUE)
L <- log2(norm + 1)
ann_all <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)
reps <- sort(unique(cd$replicate))
sample_of <- function(rep, grp) rownames(cd)[cd$replicate == rep & cd$group == grp]

# number of transcripts per gene (from the GTF; skipped if the GTF is not there)
n_tx <- tryCatch({
  gl <- readLines(gzfile("data/ref/Rattus_norvegicus.GRCr8.116.gtf.gz"))
  gl <- gl[grepl("\ttranscript\t", gl, fixed = TRUE)]
  table(sub('.*gene_id "([^"]+)".*', "\\1", gl))
}, error = function(e) NULL)
if (exists("gl", inherits = FALSE)) rm(gl); invisible(gc())

# ================================================================ A. discordant genes
to_log2 <- function(fc) sign(fc) * log2(abs(fc))
read_authors <- function(f) { a <- read.delim(gzfile(f), stringsAsFactors = FALSE); a$EnsemblGene <- sub("\\.[0-9]+$", "", a$EnsemblGene); a }
auth <- list(`50`  = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_50.tsv.gz"),
             `500` = read_authors("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_500.tsv.gz"))
match_ids <- function(a) {               # same rule as 09_validate.R
  direct <- a$EnsemblGene %in% ann_all$gene
  uniq_ours <- ann_all$gene_name[!duplicated(ann_all$gene_name) & !duplicated(ann_all$gene_name, fromLast = TRUE)]
  uniq_auth <- a$Gene_name[!duplicated(a$Gene_name) & !duplicated(a$Gene_name, fromLast = TRUE)]
  by_name <- !direct & a$Gene_name %in% uniq_ours & a$Gene_name %in% uniq_auth & !grepl("^ENSRNOG|^AABR|^AC[0-9]", a$Gene_name)
  id <- ifelse(direct, a$EnsemblGene, NA)
  id[by_name] <- ann_all$gene[match(a$Gene_name[by_name], ann_all$gene_name)]
  # a name match must not land on a gene that another authors' row already matched by ID (e.g. Ogn)
  clash <- by_name & id %in% id[direct]
  id[clash] <- NA; by_name[clash] <- FALSE
  data.frame(our_id = id, how = ifelse(direct, "id", ifelse(by_name, "name", "none")), stringsAsFactors = FALSE)
}

cat("== A. Disagreements with the authors at dose 100\n")
cat("authors' table columns:", paste(names(auth[["50"]]), collapse = ", "), "\n")
for (sz in names(auth)) {
  a <- auth[[sz]]; mm <- match_ids(a); a$our_id <- mm$our_id; a$how <- mm$how
  o <- read.delim(file.path("results/deseq2", sprintf("%sD_vs_%sA.tsv", sz, sz)), stringsAsFactors = FALSE)
  m <- merge(data.frame(gene_id = a$our_id, authors_id = a$EnsemblGene, authors_name = a$Gene_name, how = a$how,
                        authors_FC = a$FoldChange_D_vs_A, authors_l2 = to_log2(a$FoldChange_D_vs_A), authors_padj = a$padj_D_vs_A),
             o[, c("gene_id", "gene_name", "baseMean", "log2FC", "log2FC_shrunk", "FC_shrunk_signed", "padj")], by = "gene_id")
  deg_a <- m$authors_padj < 0.05 & abs(m$authors_l2) >= log2(1.5)
  deg_o <- !is.na(m$padj) & m$padj < 0.05 & abs(m$FC_shrunk_signed) >= 1.5
  m$class <- ifelse(deg_a & deg_o, "1 DEG in both",
             ifelse(deg_o & !deg_a & abs(m$authors_l2) < 0.3, "2 ours only, authors ~0",
             ifelse(deg_a & !deg_o & abs(m$log2FC_shrunk) < 0.3, "3 authors only, ours ~0",
             ifelse(deg_a | deg_o, "4 DEG in one, partial agreement", "5 neither"))))
  # per-replicate support for OUR direction: log2 of (dose 100 / control) inside each replicate
  per_rep <- sapply(reps, function(r) log2((norm[m$gene_id, sample_of(r, paste0(sz, "D"))] + 0.5) /
                                           (norm[m$gene_id, sample_of(r, paste0(sz, "A"))] + 0.5)))
  m$reps_agree_with_ours <- rowSums(sign(per_rep) == sign(m$log2FC))
  m$n_transcripts <- if (is.null(n_tx)) NA else as.integer(n_tx[m$gene_id])
  m$name_differs <- m$how == "id" & m$authors_name != m$gene_name

  cat(sprintf("\n---- %s nm: %d shared genes\n", sz, nrow(m)))
  cat(sprintf("authors' table overall: FC exactly +/-1 in %d genes; padj exactly 1 in %d genes (of %d)\n",
              sum(abs(a$FoldChange_D_vs_A) == 1, na.rm = TRUE), sum(a$padj_D_vs_A == 1, na.rm = TRUE), nrow(a)))
  s <- do.call(rbind, lapply(split(m, m$class), function(x) data.frame(
    n = nrow(x),
    median_baseMean = round(median(x$baseMean)),
    pct_matched_by_name = round(100 * mean(x$how == "name"), 1),
    pct_name_differs = round(100 * mean(x$name_differs), 1),
    median_transcripts = median(x$n_transcripts),
    pct_authors_padj_1 = round(100 * mean(x$authors_padj == 1, na.rm = TRUE), 1),
    pct_4of4_reps_agree_ours = round(100 * mean(x$reps_agree_with_ours == 4), 1))))
  print(s)
  for (cl in c("2 ours only, authors ~0", "3 authors only, ours ~0")) {
    x <- m[m$class == cl, ]
    x <- x[order(-abs(x$log2FC_shrunk - x$authors_l2)), ]
    cat(sprintf("\ntop of '%s' (largest disagreement first):\n", cl))
    show <- head(x[, c("gene_name", "authors_name", "how", "baseMean", "n_transcripts", "authors_l2", "authors_padj",
                       "log2FC_shrunk", "padj", "reps_agree_with_ours")], 15)
    show$baseMean <- round(show$baseMean); show$authors_l2 <- round(show$authors_l2, 2)
    show$log2FC_shrunk <- round(show$log2FC_shrunk, 2); show$authors_padj <- signif(show$authors_padj, 2); show$padj <- signif(show$padj, 2)
    if (nrow(show)) print(show, row.names = FALSE) else cat("  (none)\n")
  }
  write.table(m[m$class != "5 neither", ], file.path(out, sprintf("discordant_genes_%s.tsv", sz)), sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\n---- Npy, Nes, Vim in the authors' tables (searched by name AND by our Ensembl ID)\n")
for (g in c("Npy", "Nes", "Vim")) {
  our_id <- ann_all$gene[ann_all$gene_name == g]
  cat(sprintf("\n%s — our ID(s): %s | baseMean in our data: %s\n", g, paste(our_id, collapse = ","),
              paste(round(rowMeans(norm[intersect(our_id, rownames(norm)), , drop = FALSE])), collapse = ",")))
  for (sz in names(auth)) {
    a <- auth[[sz]]; hit <- a[a$Gene_name == g | a$EnsemblGene %in% our_id, c("EnsemblGene", "Gene_name", "FoldChange_D_vs_A", "padj_D_vs_A")]
    if (nrow(hit)) { cat(sz, "nm authors:\n"); print(hit, row.names = FALSE) } else cat(sz, "nm authors: not in their table\n")
  }
}

# ================================================================ B. sample signature scores
# For each replicate r, learn from the OTHER three replicates:
#   R50  = average response of 50D vs 50A,  R500 = average response of 500D vs 500A   (log2 scale)
#   genes used: expressed (mean normalised count >= 50) and |R50| >= 1 or |R500| >= 1
# Then for every treated sample s of replicate r:  d = log2(s) - log2(its own control in r), and fit
#   d = a * R50 + b * R500 + noise          (one linear regression per sample, no intercept)
#   dose_strength = a + b    (1 = as strong as a typical dose 100, 0 = no response)
#   frac_50       = a / (a + b)   (1 = pure 50 nm pattern, 0 = pure 500 nm pattern; NA if almost no response)
# Version 1 of this section used cor(d, R50 - R500). It failed its positive control: the 500 nm response is
# larger and points the same way as the 50 nm one, so ANY response correlated negatively with R50 - R500.
# The regression separates "how much response" (a + b) from "which kind" (a vs b).
cat("\n\n== B. Does each sample look like 50 nm or 500 nm? (patterns learned from the other 3 replicates)\n")
expressed <- rowMeans(norm) >= 50
resp <- function(r, grp, ctrl) L[expressed, sample_of(r, grp)] - L[expressed, sample_of(r, ctrl)]
rows <- list()
for (r in reps) {
  others <- setdiff(reps, r)
  R50  <- rowMeans(sapply(others, resp, grp = "50D",  ctrl = "50A"))
  R500 <- rowMeans(sapply(others, resp, grp = "500D", ctrl = "500A"))
  use <- abs(R50) >= 1 | abs(R500) >= 1
  for (sz in c("50", "500")) for (dose in c("B", "C", "D")) {
    d <- resp(r, paste0(sz, dose), paste0(sz, "A"))
    cf <- coef(lm(d[use] ~ 0 + R50[use] + R500[use]))
    a <- unname(cf[1]); b <- unname(cf[2])
    rows[[length(rows) + 1]] <- data.frame(
      sample = sample_of(r, paste0(sz, dose)), replicate = r, group = paste0(sz, dose),
      a_50 = round(a, 2), b_500 = round(b, 2), dose_strength = round(a + b, 2),
      frac_50 = if (a + b >= 0.2) round(a / (a + b), 2) else NA,
      n_genes = sum(use), cor_R50_R500 = round(cor(R50[use], R500[use]), 2))
  }
}
sc <- do.call(rbind, rows)
print(sc, row.names = FALSE)
write.table(sc, file.path(out, "sample_signature_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(out, "signature_map.png"), width = 1600, height = 1500, res = 200)
col <- ifelse(grepl("^500", sc$group), "#e66101", "#5e3c99")
pch <- c(B = 1, C = 17, D = 16)[substr(sc$group, nchar(sc$group), nchar(sc$group))]
plot(sc$a_50, sc$b_500, pch = pch, col = col, cex = 1.4, asp = 1,
     xlab = "a: amount of the 50 nm pattern", ylab = "b: amount of the 500 nm pattern",
     main = "Each treated sample = a x (50 nm pattern) + b x (500 nm pattern)\n(patterns learned from the other 3 replicates)", cex.main = 0.9)
abline(h = 0, v = 0, col = "grey80"); abline(0, 1, lty = 2, col = "grey60")
text(sc$a_50, sc$b_500, paste0(sc$group, ".", sub("rep", "r", sc$replicate)), pos = 4, cex = 0.55, col = col)
legend("topright", bty = "n", cex = 0.8, pch = c(16, 16, 1, 17, 16), col = c("#5e3c99", "#e66101", "grey30", "grey30", "grey30"),
       legend = c("labelled 50 nm", "labelled 500 nm", "dose 1", "dose 10", "dose 100"))
invisible(dev.off())
cat("\nsaved signature_map.png\n")
sink()
cat("done — see", out, "\n")
