# Run from the project root. Prefer fresh step-08 tables. Release archives may store them compressed.
sz <- c("50","500")
a <- lapply(sz, function(s) {
  path <- file.path("results/deseq2", sprintf("%sD_vs_%sA.tsv", s, s))
  if (!file.exists(path)) path <- paste0(path, ".gz")
  read.delim(path)
})
d <- merge(a[[1]],a[[2]],by="gene_id",suffixes=c("_50","_500"))
x <- d$log2FC_shrunk_500
y <- d$log2FC_shrunk_50
a50 <- !is.na(d$padj_50) & d$padj_50 < .05
a500 <- !is.na(d$padj_500) & d$padj_500 < .05
subsets <- list(all=rep(TRUE,nrow(d)),sig500=a500,large_sig500=a500 & abs(x)>.5,either=a50|a500,both=a50&a500)
z <- do.call(rbind,lapply(names(subsets),function(name) {
 k <- subsets[[name]]
 data.frame(subset=name,n=sum(k),median_abs_50=median(abs(y[k])),median_abs_500=median(abs(x[k])),slope=sum(x[k]*y[k])/sum(x[k]^2))
}))
if (!exists("out",inherits=FALSE)) out <- "results/amplitude"
dir.create(out,showWarnings=FALSE,recursive=TRUE)
write.table(z,file.path(out,"shrunken_amplitude.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
