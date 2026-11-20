# High-dose competitive matched-background analysis on full-precision DESeq2 contrasts.
# Run:  Rscript --vanilla scripts/22_amplitude_matched_sets.R . results/amplitude
# Needs: results/deseq2/dds.rds, results/tximport/gene_annotation.tsv,
#        results/interaction/{interaction_genes.tsv.gz,set_members.tsv}

suppressPackageStartupMessages({library(DESeq2)})
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
out  <- if (length(args) > 1) args[2] else file.path(root, "results/amplitude")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
setwd(root)

N_DRAWS <- 20000L
WINDOW  <- 50L          # +/- this many ranks when matching on the 500 nm effect
SEED    <- 20261115L

dds <- readRDS("results/deseq2/dds.rds")
ann <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors = FALSE)

# ---- within-series responses at 100 ug/mL, on the same gene universe -------------
res50  <- results(dds, contrast = c("group", "50D",  "50A"),  alpha = 0.05)
res500 <- results(dds, contrast = c("group", "500D", "500A"), alpha = 0.05)

# ---- the between-series contrast, as implemented in scripts 17-19 ---------------
inter <- results(dds,
  contrast = list(c("group_50D_vs_50A", "group_500A_vs_50A"), "group_500D_vs_50A"),
  alpha = 0.05)

stopifnot(identical(rownames(res50), rownames(inter)),
          identical(rownames(res500), rownames(inter)))

sym <- ann$gene_name[match(rownames(inter), ann$gene)]
ok  <- is.finite(inter$log2FoldChange) & is.finite(res500$log2FoldChange) &
       !is.na(sym) & nzchar(sym)

tab <- data.frame(gene = sym[ok],
                  l500 = res500$log2FoldChange[ok],
                  l50  = res50$log2FoldChange[ok],
                  diff = inter$log2FoldChange[ok],
                  stringsAsFactors = FALSE)

# ---- describe the global amplitude relationship ---------------------------------
slope <- sum(tab$l50 * tab$l500) / sum(tab$l500^2)
global <- data.frame(
  n_genes          = nrow(tab),
  median_abs_l50   = median(abs(tab$l50)),
  median_abs_l500  = median(abs(tab$l500)),
  slope_50_on_500  = slope)
write.table(global, file.path(out, "global_amplitude.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- sets ------------------------------------------------------------------------
members <- read.delim("results/interaction/set_members.tsv", stringsAsFactors = FALSE)
sets <- split(as.character(members$values), members$ind)

# Per-gene estimates underlying the high-dose complement means. Adjusted P values
# retain their original genome-wide test family in each contrast.
comp <- which(ok & sym %in% sets[["complement_activation"]])
stopifnot(length(comp) > 0L, !anyDuplicated(sym[comp]))
comp <- comp[order(sym[comp])]
complement_genes <- data.frame(
  gene_id = rownames(inter)[comp], gene_symbol = sym[comp],
  log2FC_50nm = res50$log2FoldChange[comp], padj_50nm = res50$padj[comp],
  log2FC_500nm = res500$log2FoldChange[comp], padj_500nm = res500$padj[comp],
  log2FC_difference = inter$log2FoldChange[comp], padj_difference = inter$padj[comp]
)
write.table(complement_genes, file.path(out, "complement_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- matched-background permutation ---------------------------------------------
# Null: sets of the same size whose members are drawn from genes with a similar
# 500 nm log2 fold change. This holds amplitude fixed and asks whether the set
# still differs between series more than its neighbours do.
ord  <- order(tab$l500)
diff_sorted <- tab$diff[ord]
rank_of <- setNames(seq_len(nrow(tab)), tab$gene[ord])
n <- length(diff_sorted)

set.seed(SEED)
rows <- lapply(names(sets), function(s) {
  stopifnot(!any(sets[[s]] %in% tab$gene[duplicated(tab$gene)]))
  idx <- rank_of[intersect(sets[[s]], names(rank_of))]
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(NULL)
  obs <- mean(diff_sorted[idx])
  draws <- vapply(seq_len(N_DRAWS), function(i) {
    pick <- vapply(idx, function(j)
      sample(seq(max(1L, j - WINDOW), min(n, j + WINDOW)), 1L), integer(1))
    mean(diff_sorted[pick])
  }, numeric(1))
  data.frame(set = s, n_represented = length(idx),
             observed_mean_diff = obs,
             matched_null_mean  = mean(draws),
             matched_null_sd    = sd(draws),
             excess             = obs - mean(draws),
             p_one_sided        = (sum(draws >= obs) + 1) / (N_DRAWS + 1))
})
result <- do.call(rbind, rows)
result$BH_three_sets <- p.adjust(result$p_one_sided, "BH")
write.table(result, file.path(out, "set_excess_over_amplitude.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(global); print(result)

# Same matched null for the independently sourced identity sets.
source_file <- "data/marker_sources/sets_source.tsv"
if (file.exists(source_file)) {
  definitions <- read.delim(source_file)
  defs <- subset(definitions, origin=="external" & grepl("^Panglao_",set))
  keep <- which(ok)
  pick <- order(-abs(inter$stat[keep]),rownames(inter)[keep],na.last=TRUE)
  key <- toupper(tab$gene[pick])
  chosen <- pick[!duplicated(key)]
  rank_index <- integer(n)
  rank_index[ord] <- seq_len(n)
  lookup <- setNames(rank_index[chosen],toupper(tab$gene[chosen]))
  set.seed(SEED)
  ext <- lapply(split(toupper(defs$symbol),defs$set),function(gs) {
    idx <- unname(lookup[intersect(gs,names(lookup))])
    observed <- mean(diff_sorted[idx])
    means <- numeric(N_DRAWS)
    for (j in idx) {
      neighbours <- seq(max(1L,j-WINDOW),min(n,j+WINDOW))
      means <- means+diff_sorted[sample(neighbours,N_DRAWS,replace=TRUE)]/length(idx)
    }
    data.frame(n=length(idx),observed=observed,expected=mean(means),
      excess=observed-mean(means),p=(1+sum(means>=observed))/(1+N_DRAWS))
  })
  ext <- do.call(rbind,ext)
  ext$BH_two_sets <- p.adjust(ext$p,"BH")
  write.table(data.frame(set=rownames(ext),ext),file.path(out,"external_matched.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

source("scripts/review_amplitude_descriptives.R",local=TRUE)
