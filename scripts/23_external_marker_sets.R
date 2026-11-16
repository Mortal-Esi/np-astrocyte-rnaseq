# Run from the project root, after step 08. Fixed memberships in data/marker_sources/sets_source.tsv.
suppressPackageStartupMessages({ library(DESeq2);library(clusterProfiler) })
members <- read.delim("data/marker_sources/sets_source.tsv")
external <- subset(members,origin=="external")
term <- data.frame(id=external$set,gene=toupper(external$symbol))
dds <- readRDS("results/deseq2/dds.rds")
ann <- read.delim("results/tximport/gene_annotation.tsv")
ids <- rownames(dds)
symbols <- toupper(ann$gene_name[match(ids,ann$gene)])
tests <- list(dose100_50=c("group","50D","50A"), dose100_500=c("group","500D","500A"),
              difference=list(c("group_50D_vs_50A","group_500A_vs_50A"),"group_500D_vs_50A"))
rows <- lapply(names(tests),function(name) {
 r <- results(dds,contrast=tests[[name]],alpha=.05)
 i <- which(!is.na(symbols)&nzchar(symbols)&is.finite(r$stat))
 i <- i[order(-abs(r$stat[i]),ids[i])]
 i <- i[!duplicated(symbols[i])]
 ranks <- sort(setNames(r$stat[i],symbols[i]),decreasing=TRUE)
 set.seed(20261116)
 g <- as.data.frame(GSEA(ranks,TERM2GENE=term,minGSSize=4,maxGSSize=500,pvalueCutoff=1,eps=0,seed=TRUE,verbose=FALSE))
 stopifnot(setequal(g$ID,unique(term$id)))
 data.frame(contrast=name,g[,c("ID","setSize","NES","pvalue","p.adjust","core_enrichment")])
})
out <- "results/marker_sensitivity"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(do.call(rbind,rows),file.path(out,"external_gsea.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
