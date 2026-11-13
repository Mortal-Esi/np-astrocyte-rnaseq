# Rscript --vanilla compare_counts.R PROJECT OUTPUT_DIRECTORY
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
out <- normalizePath(args[2])
setwd(args[1])
# Match published genes by the same ID / unique-name rule as step 12.
ann_all <- read.delim("results/tximport/gene_annotation.tsv", stringsAsFactors=FALSE)
published <- function(sz) {
  file <- sprintf("data/authors/GSE256038_DifferentialExpression_MicroPlastics_Astrocytes_Gill_%s.tsv.gz", sz)
  x <- read.delim(file, stringsAsFactors=FALSE)
  x$EnsemblGene <- sub("\\.[0-9]+$", "", x$EnsemblGene)
  x
}
unique_names <- function(x) !duplicated(x) & !duplicated(x,fromLast=TRUE)
match_ids <- function(a) {
  direct <- a$EnsemblGene %in% ann_all$gene
  id <- ifelse(direct,a$EnsemblGene,NA_character_)
  eligible <- !direct & unique_names(a$Gene_name) &
    a$Gene_name %in% ann_all$gene_name[unique_names(ann_all$gene_name)] &
    !grepl("^ENSRNOG|^AABR|^AC[0-9]",a$Gene_name)
  id[eligible] <- ann_all$gene[match(a$Gene_name[eligible],ann_all$gene_name)]
  id[eligible & id %in% id[direct]] <- NA_character_
  id
}
dds <- readRDS("results/deseq2/dds.rds")
plain <- dds
design(plain) <- ~group
plain <- DESeq(plain, quiet=TRUE)
rows <- list()
for (model in c("unblocked","inferred")) {
 d <- if (model=="unblocked") plain else dds
 for (sz in c("50","500")) {
  if (sz=="500") {
   d$group <- relevel(d$group,"500A")
   d <- nbinomWaldTest(d,quiet=TRUE)
  }
  coef <- paste0("group_",sz,"D_vs_",sz,"A")
  r <- results(d,name=coef,alpha=.05)
  shr <- lfcShrink(d,coef=coef,type="apeglm",quiet=TRUE)
  shared <- rownames(r) %in% match_ids(published(sz))
  for (universe in c("all","shared")) {
   use <- if (universe=="all") rep(TRUE,nrow(r)) else shared
   sig <- use & !is.na(r$padj) & r$padj<.05
   rows[[length(rows)+1]] <- data.frame(model,size_nm=sz,universe,genes=sum(use),with_padj=sum(use & !is.na(r$padj)),padj05=sum(sig),padj05_apeglm_FC1.5=sum(sig & abs(shr$log2FoldChange)>=log2(1.5),na.rm=TRUE))
  }
 }
}
write.table(do.call(rbind,rows),file.path(out,"model_counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# The step-08 summary rounds before counting. This audit applies both cutoffs
# at full precision. The separate permutation audit uses contrast= for all fits.
