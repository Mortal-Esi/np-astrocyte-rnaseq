# Run with Rscript --vanilla block_control.R PROJECT OUTPUT_DIRECTORY
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
out <- normalizePath(args[2])
setwd(root)
samples <- read.csv("samples.csv", row.names=1)
samples$group <- factor(samples$group, levels=c("50A","50B","50C","50D","500A","500B","500C","500D"))
samples$replicate <- factor(samples$replicate)
txi <- readRDS("results/tximport/txi.rds")
stopifnot(identical(colnames(txi$counts), rownames(samples)))
base <- DESeqDataSetFromTximport(txi, colData=samples, design=~group)
base <- base[rowSums(counts(base)>=10)>=4,]
fit <- function(d) DESeq(d,quiet=TRUE)
sizes <- c("50","500")
count_sig <- function(d) sapply(sizes,function(sz) {
 r <- results(d,contrast=c("group",paste0(sz,c("D","A"))),alpha=.05)
 sum(r$padj<.05,na.rm=TRUE)
})
write_tab <- function(x,name) write.table(x,file.path(out,name),sep="\t",quote=FALSE,row.names=FALSE)
original <- base
design(original) <- ~replicate+group
paired <- fit(original)
plain <- fit(base)
counts <- rbind(unblocked=count_sig(plain), inferred=count_sig(paired))
write_tab(data.frame(model=rownames(counts),counts,check.names=FALSE),"control_baselines.tsv")
set.seed(20261113)
B <- 39L
draws <- replicate(B,as.character(unsplit(lapply(split(samples$replicate,samples$group),sample),samples$group)))
write_tab(data.frame(run=rownames(samples),group=samples$group,original=samples$replicate,draws),"permuted_labels.tsv")
# For DE counts the two sizes use a contrast under the same parameterization throughout.
for (i in seq_len(B)) {
 d <- original
 d$replicate <- factor(draws[,i],levels=levels(samples$replicate))
 stopifnot(all(table(d$replicate,d$group)==1))
 n <- count_sig(fit(d))
 row <- data.frame(draw=i,n50=n[1],n500=n[2])
 write.table(row,file.path(out,"permuted_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=i==1,append=i>1)
}
