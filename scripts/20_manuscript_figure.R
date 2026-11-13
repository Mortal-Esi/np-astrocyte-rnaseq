library(ggplot2)
library(patchwork)
a <- read.delim("results/dose_response/block_scores.tsv")
b <- read.delim("results/dose_sensitivity/pathway_stability.tsv")
labs_set <- c(ependymal_cilia="Ciliary / ependymal",complement_activation="Complement activation")
a$path <- factor(labs_set[a$set],levels=unname(labs_set))
b$path <- factor(labs_set[b$ID],levels=unname(labs_set))
a$dose <- factor(a$dose,levels=c(1,10,100))
b$dose <- factor(b$dose,levels=c(1,10,100))
b$model <- factor(b$fit,levels=rev(c("full","unblocked",paste0("omit_rep",1:4))),labels=rev(c("Full","Group only",paste("Omit",1:4))))
style <- theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),plot.title=element_text(face="bold"),legend.position="bottom")
p <- ggplot(a,aes(dose,mean_log2_difference,group=block,color=block,shape=block))+
 geom_hline(yintercept=0,linetype=2,color="grey60")+geom_line(linewidth=.55)+geom_point(size=2.5)+facet_wrap(~path)+
 scale_color_manual(values=c("#007469","#6649A0","#A4491C","#285887"))+
 labs(title="A    Expression differences within inferred blocks",x="Nominal mass concentration (µg/mL)",y="Mean difference in log2 response",color="Block",shape="Block")+style
q <- ggplot(b,aes(dose,model,fill=NES))+geom_tile(colour="white",linewidth=.8)+
 geom_point(data=subset(b,FDR_six_tests<.05),shape=16,size=2)+facet_wrap(~path)+
 scale_fill_gradient2(low="#ab6e58",mid="#f8f8f4",high="#63aea4",midpoint=0,limits=c(-2.2,2.2))+
 labs(title="B    Enrichment across model fits",subtitle="Dot = BH-adjusted P < 0.05 across six tests per fit",x="Nominal mass concentration (µg/mL)",y=NULL,fill="NES")+style+theme(panel.grid=element_blank())
fig <- p/q+plot_layout(heights=c(1,1.07))
ggsave("manuscript/Figure1.pdf",fig,width=8,height=8,device="pdf")
ggsave("manuscript/Figure1.png",fig,width=8,height=8,dpi=220)
