# Run: Rscript 19_dose_sensitivity.R PROJECT ABSOLUTE_OUTPUT_DIRECTORY
suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Rn.eg.db)
  library(clusterProfiler)
})
args <- commandArgs(TRUE)
if (length(args)) setwd(args[1])
out <- if (length(args)>1) args[2] else "results/dose_sensitivity"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_table <- function(x,file) write.table(x,file.path(out,file),sep="\t",quote=FALSE,row.names=FALSE)
s <- read.csv("samples.csv",row.names=1)
s$group <- relevel(factor(s$group),"50A")
s$replicate <- factor(s$replicate)
stopifnot(all(table(s$group,s$replicate)==1))
txi <- readRDS("results/tximport/txi.rds")
stopifnot(identical(colnames(txi$counts),rownames(s)))
full <- readRDS("results/deseq2/dds.rds")
stopifnot(identical(colnames(full),rownames(s)),all(full$group==as.character(s$group)),
          all(full$replicate==as.character(s$replicate)),identical(all.vars(design(full)),c("replicate","group")))
ids <- rownames(full)
ann <- read.delim("results/tximport/gene_annotation.tsv")
genes <- ann$gene_name[match(ids,ann$gene)]
sets_file <- read.delim("results/interaction/set_members.tsv")
sets <- split(as.character(sets_file$values),sets_file$ind)[c("ependymal_cilia","complement_activation")]
mapping <- stack(sets)
evaluate <- function(d,label) {
  set.seed(20261113)
  rows <- lapply(seq_along(c("B","C","D")),function(i) {
    tag <- c("B","C","D")[i]
    contrast <- list(c(paste0("group_50",tag,"_vs_50A"),"group_500A_vs_50A"),paste0("group_500",tag,"_vs_50A"))
    r <- results(d,contrast=contrast,alpha=.05)
    ok <- which(is.finite(r$stat) & !is.na(genes) & nzchar(genes))
    ok <- ok[order(-abs(r$stat[ok]),ids[ok])]
    ok <- ok[!duplicated(genes[ok])]
    ranks <- sort(setNames(r$stat[ok],genes[ok]),decreasing=TRUE)
    g <- as.data.frame(GSEA(ranks,TERM2GENE=data.frame(term=mapping$ind,gene=mapping$values),
      minGSSize=4,maxGSSize=500,pvalueCutoff=1,eps=0,seed=TRUE,verbose=FALSE))
    stopifnot(setequal(g$ID,names(sets)))
    data.frame(fit=label,dose=c(1,10,100)[i],g[,c("ID","setSize","NES","pvalue","core_enrichment")])
  })
  tab <- do.call(rbind,rows)
  tab$FDR_six_tests <- p.adjust(tab$pvalue,"BH")
  tab
}
refit <- function(use,formula,label) {
  x <- txi
  for (k in c("counts","length","abundance")) x[[k]] <- txi[[k]][ids,use,drop=FALSE]
  d <- DESeqDataSetFromTximport(x,droplevels(s[use,,drop=FALSE]),formula)
  evaluate(DESeq(d,quiet=TRUE),label)
}
fits <- list(evaluate(full,"full"),refit(rep(TRUE,nrow(s)),~group,"unblocked"))
for (b in levels(s$replicate)) fits[[paste0("omit_",b)]] <- refit(s$replicate!=b,~replicate+group,paste0("omit_",b))
save_table(do.call(rbind,fits),"pathway_stability.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
