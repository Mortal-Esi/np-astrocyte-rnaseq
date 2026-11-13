# From the repo root: Rscript scripts/17_interaction_robustness.R
# Optional arguments: project directory, output directory
suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Rn.eg.db)
  library(clusterProfiler)
})
args <- commandArgs(TRUE)
if (length(args)) setwd(args[1])
out <- if (length(args)>1) args[2] else "results/interaction"
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
contrast <- list(c("group_50D_vs_50A","group_500A_vs_50A"),"group_500D_vs_50A")
extract <- function(d) as.data.frame(results(d,contrast=contrast,alpha=.05))
# Rebuild from tximport so that each refit estimates its own normalization and dispersions.
refit <- function(use,formula) {
  x <- txi
  for (k in c("counts","length","abundance")) x[[k]] <- txi[[k]][ids,use,drop=FALSE]
  d <- DESeqDataSetFromTximport(x,droplevels(s[use,,drop=FALSE]),formula)
  extract(DESeq(d,quiet=TRUE))
}
blocks <- levels(s$replicate)
fits <- list(full=extract(full), unblocked=refit(rep(TRUE,nrow(s)),~group))
for (b in blocks) fits[[paste0("omit_",b)]] <- refit(s$replicate!=b,~replicate+group)
stopifnot(all(vapply(fits,function(r) identical(rownames(r),ids),logical(1))))
est <- fits$full
tab <- data.frame(gene_id=ids,gene=genes,baseMean=est$baseMean)
for (k in names(fits)) {
  r <- fits[[k]]
  tab[[paste0(k,"_log2")]] <- r$log2FoldChange
  tab[[paste0(k,"_padj")]] <- r$padj
}
tab$ci95_low <- est$log2FoldChange-qnorm(.975)*est$lfcSE
tab$ci95_high <- est$log2FoldChange+qnorm(.975)*est$lfcSE
L <- log2(counts(full,normalized=TRUE)+.5)
for (b in blocks) {
  v <- function(g) L[,s$replicate==b & s$group==g]
  tab[[paste0(b,"_observed_log2")]] <- (v("50D")-v("50A"))-(v("500D")-v("500A"))
}
sig <- function(p) !is.na(p) & p<.05
direction <- function(x) is.finite(x) & is.finite(tab$full_log2) & x*tab$full_log2>0
tab$loo_same_sign <- Reduce("&",lapply(fits[-(1:2)],function(r) direction(r$log2FoldChange)))
tab$loo_all_significant <- Reduce("&",lapply(fits[-(1:2)],function(r) sig(r$padj)))
tab$blocks_same_sign <- Reduce("&",lapply(tab[paste0(blocks,"_observed_log2")],direction))
tab$unblocked_same_sign <- direction(tab$unblocked_log2)
stable <- sig(tab$full_padj) & tab$loo_same_sign & tab$blocks_same_sign
tab$direction_stable <- stable
tab$strict_stable <- stable & tab$loo_all_significant
# Restrict reporting after genome-wide adjustment; do not recalculate FDR on selected genes.
summary <- do.call(rbind,lapply(c(0,10,20,50),function(t) {
 keep <- tab$baseMean>=t
 n <- function(x) sum(keep & x,na.rm=TRUE)
 data.frame(min_baseMean=t,retained=sum(keep),with_padj=n(!is.na(tab$full_padj)),
 full_hits=n(sig(tab$full_padj)),unblocked_hits=n(sig(tab$unblocked_padj)),
 shared_hits_same_direction=n(sig(tab$full_padj)&sig(tab$unblocked_padj)&tab$unblocked_same_sign),
 full_hits_same_sign_loo=n(sig(tab$full_padj)&tab$loo_same_sign),
 full_hits_same_sign_blocks=n(sig(tab$full_padj)&tab$blocks_same_sign),
 direction_stable=n(stable),strict_stable=n(tab$strict_stable))
}))
save_table(summary,"sensitivity.tsv")
f <- gzfile(file.path(out,"interaction_genes.tsv.gz"),"wt")
write.table(tab,f,sep="\t",row.names=FALSE,quote=FALSE)
close(f)
key <- c("Aqp4","Glul","Gja1","Atp1a2","Foxj1","C3","Ccl2","Gpnmb","Npy","Snca","Pink1","Uchl1","Dnajc6")
save_table(tab[match(key,tab$gene),],"key_genes.tsv")
# Targeted follow-up: run the same three prespecified sets in every fit.
# BH family = these three sets per fit. Original genome-wide GO FDR has a different denominator.
sets <- list(homeostatic_astro=c("Aqp4","Slc1a2","Slc1a3","Glul","Kcnj10","Aldh1l1","Gja1","Sparcl1","Atp1a2","Slc4a4"),
 ependymal_cilia=c("Foxj1","Dnah5","Dnah9","Rsph1","Rsph4a","Ccdc153","Tmem212","Dynlrb2","Spag6","Sntn"),
 complement_activation=unique(na.omit(AnnotationDbi::select(org.Rn.eg.db,keys="GO:0006956",keytype="GOALL",columns="SYMBOL")$SYMBOL)))
save_table(stack(sets),"set_members.tsv")
set.seed(20261113)
path_results <- lapply(names(fits),function(k) {
  r <- fits[[k]]
  keep <- is.finite(r$stat) & !is.na(genes) & nzchar(genes)
  stats <- sort(setNames(r$stat[keep],genes[keep]),decreasing=TRUE)
  stats <- stats[order(-abs(stats),names(stats))]
  stats <- sort(stats[!duplicated(names(stats))],decreasing=TRUE)
  mapping <- stack(sets)
  g <- as.data.frame(GSEA(stats, TERM2GENE=data.frame(term=mapping$ind,gene=mapping$values), minGSSize=4,maxGSSize=500,pvalueCutoff=1,eps=0,seed=TRUE,verbose=FALSE))
  stopifnot(setequal(g$ID,names(sets)))
  data.frame(fit=k,pathway=g$ID,size=g$setSize,NES=g$NES,pval=g$pvalue,padj=g$p.adjust)
})
save_table(do.call(rbind,path_results),"targeted_sets.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

# Model-based, pointwise CIs and the four descriptive within-block estimates.
show <- tab[match(c("Aqp4","Glul","Gja1","Atp1a2","Foxj1","C3","Gpnmb","Npy"),tab$gene),]
draw <- function() {
  par(mar=c(4.4,5.2,4,1))
  y <- rev(seq_len(nrow(show)))
  block_est <- as.matrix(show[paste0(blocks,"_observed_log2")])
  lim <- range(show$ci95_low, show$ci95_high,block_est,finite=TRUE)
  plot(NA,xlim=lim,ylim=c(0.3,max(y)+.8),yaxt="n",ylab="",xlab="Difference in log2 response (50 nm minus 500 nm)",
       main="High-dose response at equal mass concentration",bty="n")
  axis(2,at=y,labels=show$gene,las=1,tick=FALSE)
  abline(v=0,lty=2,col="grey55")
  segments(show$ci95_low,y+.12,show$ci95_high,y+.12,lwd=2,col="#243B53")
  points(show$full_log2,y+.12,pch=18,col="#243B53")
  for (i in seq_along(blocks)) points(block_est[,i],y-.22+(i-1)*.06,pch=i,col="#A94C13")
  legend("top",legend=c("Model (95% Wald CI)",blocks),pch=c(18,1:4),col=c("#243B53",rep("#A94C13",4)),horiz=TRUE,bty="n",cex=.78)
}
pdf(file.path(out,"selected_interactions.pdf"),width=8,height=5.5)
draw()
invisible(dev.off())
