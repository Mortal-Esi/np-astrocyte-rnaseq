# 06_tximport.R — Salmon transcript estimates → gene-level count matrix (tximport)
#
# Run from the project root. In RStudio: open np-astrocyte-rnaseq.Rproj, then
#   source("scripts/06_tximport.R", echo = TRUE)
#
# Inputs : samples.csv, results/salmon/<run>/quant.sf, the Ensembl GTF
# Outputs: results/tximport/txi.rds             (tximport object → input for DESeq2)
#          results/tximport/gene_counts.csv.gz  (rounded counts, genes × samples)
#          results/tximport/gene_annotation.tsv (gene_id, gene_name, gene_biotype)
#          results/tximport/tximport_checks.txt (all checks printed below)

suppressPackageStartupMessages(library(tximport))

out_dir <- "results/tximport"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "tximport_checks.txt")
sink(log_file, split = TRUE)   # everything printed below also goes into the log file

cat("== 1. Samples and Salmon files\n")
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
files <- file.path("results/salmon", samples$run, "quant.sf")
names(files) <- samples$run
if (!all(file.exists(files))) {
  stop("missing quant.sf for: ", paste(samples$run[!file.exists(files)], collapse = ", "))
}
cat(length(files), "samples found\n")

hdr <- names(read.delim(files[1], nrows = 1))
need <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
cat("quant.sf columns:", paste(hdr, collapse = ", "), "\n")
if (!all(need %in% hdr)) stop("quant.sf does not have the columns tximport expects: ", paste(setdiff(need, hdr), collapse = ", "))

cat("\n== 2. Transcript → gene table from the GTF\n")
gtf <- "data/ref/Rattus_norvegicus.GRCr8.116.gtf.gz"
lines <- readLines(gzfile(gtf))
lines <- lines[!startsWith(lines, "#")]
fields3 <- sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*$", "\\1", lines)      # column 3 = feature type
tx_lines <- lines[fields3 == "transcript"]
rm(lines, fields3); invisible(gc())
attrs <- sub("^([^\t]*\t){8}", "", tx_lines)                          # column 9 = attributes
get_attr <- function(key) {
  m <- regmatches(attrs, regexec(paste0("(^|; )", key, ' "([^"]*)"'), attrs))
  vapply(m, function(x) if (length(x)) x[3] else NA_character_, character(1))
}
tx_table <- data.frame(
  tx           = get_attr("transcript_id"),
  gene         = get_attr("gene_id"),
  gene_name    = get_attr("gene_name"),
  gene_biotype = get_attr("gene_biotype"),
  stringsAsFactors = FALSE
)
tx_table$gene_name[is.na(tx_table$gene_name)] <- tx_table$gene[is.na(tx_table$gene_name)]
cat(nrow(tx_table), "transcripts,", length(unique(tx_table$gene)), "genes in the GTF\n")

cat("\n== 3. Do Salmon's transcript names match the GTF? (version suffix check)\n")
q_names <- read.delim(files[1])$Name
cat("example Salmon name:", q_names[1], "| example GTF id:", tx_table$tx[1], "\n")
q_nover <- sub("\\.[0-9]+$", "", q_names)
missing <- setdiff(q_nover, tx_table$tx)
cat(length(q_names), "transcripts in quant.sf;", length(missing), "not found in the GTF after removing versions\n")
if (length(missing) > 0) cat("first missing:", head(missing, 5), "\n")

cat("\n== 4. tximport (gene level)\n")
txi <- tximport(files, type = "salmon",
                tx2gene = tx_table[, c("tx", "gene")],
                ignoreTxVersion = TRUE,
                dropInfReps = TRUE)        # no bootstraps were made (--numBootstraps 0)
counts <- txi$counts
cat("count matrix:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("countsFromAbundance:", txi$countsFromAbundance, "(no = raw estimated counts + length offset for DESeq2)\n")

cat("\n== 5. Reads per sample: tximport total vs Salmon 'mapped'\n")
summ_file <- "results/salmon/salmon_summary.tsv"
if (file.exists(summ_file)) {
  summ <- read.delim(summ_file)
  cmp <- data.frame(sample = colnames(counts),
                    tximport_total = round(colSums(counts)),
                    salmon_mapped = summ$mapped[match(colnames(counts), summ$sample)])
  cmp$ratio <- round(cmp$tximport_total / cmp$salmon_mapped, 4)
  print(cmp, row.names = FALSE)
}

cat("\n== 6. How many genes are actually expressed?\n")
zero <- sum(rowSums(counts) == 0)
expressed <- rowSums(counts >= 10) >= 4      # ≥10 reads in ≥4 samples (smallest group size)
cat("genes with 0 reads in all samples:", zero, "\n")
cat("genes with ≥10 reads in ≥4 samples:", sum(expressed), "\n")
gene_ann <- unique(tx_table[, c("gene", "gene_name", "gene_biotype")])
gene_ann <- gene_ann[!duplicated(gene_ann$gene), ]
bt <- gene_ann$gene_biotype[match(rownames(counts)[expressed], gene_ann$gene)]
cat("biotypes of expressed genes:\n"); print(sort(table(bt), decreasing = TRUE)[1:min(8, length(unique(bt)))])

cat("\n== 7. Cross-check with Salmon's own gene-level file (first sample)\n")
gfile <- file.path("results/salmon", samples$run[1], "quant.genes.sf")
if (file.exists(gfile)) {
  g <- read.delim(gfile)
  common <- intersect(rownames(counts), g$Name)
  a <- counts[common, 1]; b <- g$NumReads[match(common, g$Name)]
  cat(length(common), "genes in common; correlation =", round(cor(a, b), 6),
      "; max |difference| =", round(max(abs(a - b)), 3), "reads\n")
} else cat("quant.genes.sf not found — skipped\n")

cat("\n== 8. Biological sanity check: are these astrocytes?\n")
markers <- c(Gfap = "astrocyte", Aqp4 = "astrocyte", Slc1a3 = "astrocyte", Aldh1l1 = "astrocyte",
             Rbfox3 = "neuron", Syp = "neuron", Mbp = "oligodendrocyte", Plp1 = "oligodendrocyte",
             Aif1 = "microglia", Ptprc = "immune")
ctrl <- samples$run[samples$dose_ugml == 0]
mk <- data.frame(gene = names(markers), cell_type = markers,
                 mean_counts_controls = sapply(names(markers), function(s) {
                   id <- gene_ann$gene[match(s, gene_ann$gene_name)]
                   if (is.na(id) || !id %in% rownames(counts)) NA else round(mean(counts[id, ctrl]))
                 }))
print(mk, row.names = FALSE)

cat("\n== 9. Save\n")
saveRDS(txi, file.path(out_dir, "txi.rds"))
cnt_out <- data.frame(gene_id = rownames(counts), round(counts), check.names = FALSE)
gz <- gzfile(file.path(out_dir, "gene_counts.csv.gz"), "w"); write.csv(cnt_out, gz, row.names = FALSE); close(gz)
write.table(gene_ann, file.path(out_dir, "gene_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("saved txi.rds, gene_counts.csv.gz, gene_annotation.tsv in", out_dir, "\n")
sink()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("log:", log_file, "\n")
