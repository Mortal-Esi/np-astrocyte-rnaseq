# Rscript --vanilla 18_dose_response.R PROJECT OUTPUT_DIR
suppressPackageStartupMessages({library(DESeq2);library(clusterProfiler)})
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
save_tab <- function(x,n) write.table(x,file.path(out,n),sep="\t",row.names=FALSE,quote=FALSE)
d <- readRDS(file.path(root,"results/deseq2/dds.rds"))
d$series <- factor(sub("[ABCD]$","",d$group))
d$dose <- factor(d$dose_ugml,levels=c(0,1,10,100))
stopifnot(all(table(d$series,d$dose,d$replicate)==1))
ann <- read.delim(file.path(root,"results/tximport/gene_annotation.tsv"))
sym <- ann$gene_name[match(rownames(d),ann$gene)]
sets_file <- read.delim(file.path(root,"results/interaction/set_members.tsv"))
target <- c("ependymal_cilia","complement_activation")
sets <- split(as.character(sets_file$values),sets_file$ind)[target]
# group spans all eight cells; the reduced design retains the two main effects.
factorial <- d
design(factorial) <- ~replicate+series*dose
lrt <- nbinomLRT(factorial,reduced=~replicate+series+dose,quiet=TRUE)
omni <- as.data.frame(results(lrt,alpha=.05))
genes <- data.frame(id=rownames(d),symbol=sym,LRT_3df=omni$stat,omnibus_p=omni$pvalue,omnibus_FDR=omni$padj)
norm <- log2(counts(d,normalized=TRUE)+.5)
score_rows <- list()
member_rows <- list()
gsea_rows <- list()
jack <- NULL
set.seed(20261113)
run_gsea <- function(gl,sets) {
  m <- stack(sets)
  g <- as.data.frame(GSEA(gl,TERM2GENE=data.frame(term=m$ind,gene=m$values),
      minGSSize=4,maxGSSize=500,pvalueCutoff=1,eps=0,seed=TRUE,verbose=FALSE))
  stopifnot(setequal(g$ID,names(sets)))
  g[,c("ID","setSize","NES","pvalue","core_enrichment")]
}
for (i in seq_along(c("B","C","D"))) {
  tag <- c("B","C","D")[i]
  dose <- c(1,10,100)[i]
  terms <- list(c(paste0("group_50",tag,"_vs_50A"),"group_500A_vs_50A"),paste0("group_500",tag,"_vs_50A"))
  r <- results(d,contrast=terms,alpha=.05)
  genes[[paste0("dose",dose,"_log2")]] <- r$log2FoldChange
  genes[[paste0("dose",dose,"_FDR")]] <- r$padj
  ok <- which(is.finite(r$stat) & !is.na(sym) & nzchar(sym))
  ok <- ok[order(-abs(r$stat[ok]),rownames(d)[ok])]
  ok <- ok[!duplicated(sym[ok])]
  gl <- sort(setNames(r$stat[ok],sym[ok]),decreasing=TRUE)
  gsea_rows[[i]] <- data.frame(dose,run_gsea(gl,sets))
  if (dose==100) {
    variants <- list()
    for (k in names(sets)) for (g in intersect(sets[[k]],names(gl)))
      variants[[paste(k,g,sep=":")]] <- setdiff(sets[[k]],g)
    jack <- run_gsea(gl,variants)
    jack$FDR_across_omissions <- p.adjust(jack$pvalue,"BH")
  }
  for (k in names(sets)) {
    ids <- ok[sym[ok] %in% sets[[k]]]
    for (b in levels(d$replicate)) {
      get_col <- function(gr) norm[ids,d$group==gr & d$replicate==b]
      delta <- (get_col(paste0("50",tag))-get_col("50A"))-(get_col(paste0("500",tag))-get_col("500A"))
      member_rows[[length(member_rows)+1]] <- data.frame(set=k,dose,block=b,gene=sym[ids],delta)
      score_rows[[length(score_rows)+1]] <- data.frame(set=k,dose,block=b,mean_log2_difference=mean(delta))
    }
  }
}
g <- do.call(rbind,gsea_rows)
g$FDR_six_tests <- p.adjust(g$pvalue,"BH")
save_tab(g,"pathway_doses.tsv")
save_tab(jack,"omit_one_gene.tsv")
members <- do.call(rbind,member_rows)
save_tab(members,"member_block_effects.tsv")
scores <- do.call(rbind,score_rows)
save_tab(scores,"block_scores.tsv")
con <- gzfile(file.path(out,"gene_tests.tsv.gz"),"wt")
write.table(genes,con,sep="\t",row.names=FALSE,quote=FALSE)
close(con)
pdf(file.path(out,"dose_response.pdf"),width=10,height=7)
par(mfrow=c(2,2),mar=c(4,4.3,3,1))
for (k in names(sets)) {
  x <- subset(scores,set==k)
  plot(1:3,rep(0,3),type="n",ylim=range(c(0,x$mean_log2_difference)),xaxt="n",
       xlab="Dose (ug/mL)",ylab="Mean log2 response difference",main=gsub("_"," ",k))
  axis(1,1:3,c(1,10,100)); abline(h=0,lty=3)
  cols <- c("#176D74","#77548F","#AF5724","#475C9F")
  for (i in 1:4) lines(1:3,subset(x,block==levels(d$replicate)[i])$mean_log2_difference,type="b",col=cols[i],pch=i)
  legend("topleft",levels(d$replicate),col=cols,pch=1:4,bty="n",cex=.72)
  m <- subset(members,set==k & dose==100)
  avg <- sort(tapply(m$delta,m$gene,mean))
  par(mar=c(4,5.8,3,1))
  plot(avg,seq_along(avg),xlim=range(c(0,m$delta)),yaxt="n",ylab="",pch=18,
       xlab="50 nm minus 500 nm (log2)",main="Individual genes at 100 ug/mL")
  axis(2,seq_along(avg),names(avg),las=1,cex.axis=.68)
  abline(v=0,lty=3)
  for(i in 1:4) {z<-subset(m,block==levels(d$replicate)[i]);points(z$delta,match(z$gene,names(avg)),col=cols[i],pch=i,cex=.68)}
  par(mar=c(4,4.3,3,1))
}
invisible(dev.off())
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
