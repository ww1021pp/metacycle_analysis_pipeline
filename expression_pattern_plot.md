# expression pattern plot with R for 3 replicants for different condition and error bar in the plots 


## for the inout file are homer expression table 
the format is (transcripts*sample) with transcripts annotation information

```
rm(list=ls())
setwd("~/Documents/DY_RNA_Seq_4_8_19/4_B6_129_NC_HFD_RNAseq-206142945/circadian_analysis1")

#raw <- read.delim("../4_homer/normlization.txt", sep="\t", fill=T, h=T)
raw <- read.delim("../4_homer/quatile.txt",fill=T, h=T)

colnames(raw)=gsub("\\.reads.*","",colnames(raw))
colnames(raw)[1]="Transcript_ID"
tags<-c("NC\\.129","HF\\.129","NC.B6","HF.B6")
raw$Symbol<-sapply(strsplit(as.character(raw$Annotation), "\\|"), "[", 1)
raw<-raw[which(raw$Symbol !="0.000"),]
head(raw)
dim(raw)

library(ggpubr)
library(stringr)
library(reshape2)
library(matrixStats)
raw_data= raw
raw=raw[,-56]

sd=data.frame(
  
  ZT2_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT2\\.",colnames(raw))]),
  ZT6_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT6\\.",colnames(raw))]),
  ZT10_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT10\\.",colnames(raw))]),
  ZT14_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT14\\.",colnames(raw))]),
  ZT18_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT18\\.",colnames(raw))]),
  ZT22_NC_B6 =rowMeans(raw[,grepl(pattern = "NC\\.B6\\.ZT22\\.",colnames(raw))]),
  ZT2_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT2\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT2\\.",colnames(raw)))),
  ZT6_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT6\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT6\\.",colnames(raw)))),
  ZT10_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT10\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT10\\.",colnames(raw)))),
  ZT14_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT14\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT14\\.",colnames(raw)))),
  ZT18_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT18\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT18\\.",colnames(raw)))),
  ZT22_NC_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.B6\\.ZT22\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.B6\\.ZT22\\.",colnames(raw)))),
  ZT2_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT2\\.",colnames(raw))]),
  ZT6_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT6\\.",colnames(raw))]),
  ZT10_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT10\\.",colnames(raw))]),
  ZT14_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT14\\.",colnames(raw))]),
  ZT18_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT18\\.",colnames(raw))]),
  ZT22_HF_B6 =rowMeans(raw[,grepl(pattern = "HF\\.B6\\.ZT22\\.",colnames(raw))]),
  ZT2_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT2\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT2\\.",colnames(raw)))),
  ZT6_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT6\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT6\\.",colnames(raw)))),
  ZT10_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT10\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT10\\.",colnames(raw)))),
  ZT14_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT14\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT14\\.",colnames(raw)))),
  ZT18_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT18\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT18\\.",colnames(raw)))),
  ZT22_HF_B6_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.B6\\.ZT22\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.B6\\.ZT22\\.",colnames(raw)))),
  ZT2_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT2\\.",colnames(raw))]),
  ZT6_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT6\\.",colnames(raw))]),
  ZT10_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT10\\.",colnames(raw))]),
  ZT14_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT14\\.",colnames(raw))]),
  ZT18_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT18\\.",colnames(raw))]),
  ZT22_NC_129 =rowMeans(raw[,grepl(pattern = "NC\\.129\\.ZT22\\.",colnames(raw))]),
  ZT2_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT2\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT2\\.",colnames(raw)))),
  ZT6_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT6\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT6\\.",colnames(raw)))),
  ZT10_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT10\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT10\\.",colnames(raw)))),
  ZT14_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT14\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT14\\.",colnames(raw)))),
  ZT18_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT18\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT18\\.",colnames(raw)))),
  ZT22_NC_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "NC\\.129\\.ZT22\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "NC\\.129\\.ZT22\\.",colnames(raw)))),
  ZT2_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT2\\.",colnames(raw))]),
  ZT6_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT6\\.",colnames(raw))]),
  ZT10_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT10\\.",colnames(raw))]),
  ZT14_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT14\\.",colnames(raw))]),
  ZT18_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT18\\.",colnames(raw))]),
  ZT22_HF_129 =rowMeans(raw[,grepl(pattern = "HF\\.129\\.ZT22\\.",colnames(raw))]),
  ZT2_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT2\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT2\\.",colnames(raw)))),
  ZT6_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT6\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT6\\.",colnames(raw)))),
  ZT10_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT10\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT10\\.",colnames(raw)))),
  ZT14_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT14\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT14\\.",colnames(raw)))),
  ZT18_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT18\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT18\\.",colnames(raw)))),
  ZT22_HF_129_sd =rowSds(as.matrix(raw[,grepl(pattern = "HF\\.129\\.ZT22\\.",colnames(raw))]))/sqrt(sum(grepl(pattern = "HF\\.129\\.ZT22\\.",colnames(raw))))
   )
head(sd)
dim(sd)
sd=cbind(raw$Transcript_ID,raw$Symbol,sd)
colnames(sd)[1:2]=c("Transcript_ID","Symbol")



#####fileter low expression genes by condition avarage expression gt>1, AT LEAST  half condition expression >filter_value#####
#raw[,"count"] <- rowSums(sd[,]>3)
#keep=raw$count>12
#raw <- raw[keep,];dim(raw)
#sd<-sd[keep,];head(sd);dim(sd)

library(R.utils)

core_genes<-read.csv("17clocklist.txt",stringsAsFactors = FALSE,header = F)

core_genes$mouse<-tolower(core_genes$V1)
tar_genes<-capitalize(core_genes$mouse)
tar_genes<-c("Sntg2","Taok3","Atf1")
expressed_genes<-intersect(tar_genes,raw$Symbol)
genes_matrix<-raw[which(raw$Symbol %in% expressed_genes),]
genes_matrix<-genes_matrix[order(genes_matrix$Symbol,genes_matrix$Length,decreasing=TRUE),]
Trans_ID=genes_matrix$Transcript_ID

library(reshape2)

rownames(genes_matrix) <- genes_matrix[,1]
name <- colnames(genes_matrix)
name_exp<-str_detect(name,'\\.ZT[0-9].*')
dat_exp <- genes_matrix[,name_exp]
##########################gene_expression line plot ############

dat_exp <- cbind(genes_matrix[,c("Transcript_ID","Symbol")],dat_exp)
colnames(dat_exp)[1:2]=c("Transcript_ID","Symbol")
dat_plot <- reshape2::melt(dat_exp,id  = c("Transcript_ID","Symbol"))
head(dat_plot)
dat_plot$group<-factor(sapply(strsplit(as.character(dat_plot$variable), "\\.ZT"), "[", 1),levels = c("NC.B6","HF.B6","NC.129","HF.129"),
                       labels = c("NC_B6","HF_B6","NC_129","HF_129"))


dat_plot$time<-factor(gsub(".*ZT","",gsub("\\.R.*","",dat_plot$variable)),levels = c("2","6","10","14","18","22"))
head(dat_plot)
dat_plot$Transcript_ID=as.factor(dat_plot$Transcript_ID)
########################################line_plot for each genes and time####
#############################################################################
mean_Exp<-sd[which(raw$Symbol %in% tar_genes),grepl(pattern = "^ZT.*[69]$",colnames(sd))]
mean_Exp<-cbind(sd[which(raw$Symbol %in% tar_genes),c("Transcript_ID","Symbol")],mean_Exp)
mean_plot <- reshape2::melt(mean_Exp,id  = c("Transcript_ID","Symbol"),value.name = "mean")
mean_plot$group=factor(gsub("ZT[0-9]+_","",mean_plot$variable),levels = c("NC_B6","HF_B6","NC_129","HF_129"))
mean_plot$time=factor(gsub("ZT","",gsub("_.*","",mean_plot$variable)),levels = c("2","6","10","14","18","22"))
head(mean_plot)

Exp_sd<-sd[which(raw$Symbol %in% tar_genes),grepl(pattern = "sd",colnames(sd))]
Exp_sd<-cbind(sd[which(raw$Symbol %in% tar_genes),c("Transcript_ID","Symbol")],Exp_sd)
sd_plot <- reshape2::melt(Exp_sd,id  = c("Transcript_ID","Symbol"),value.name = "sd");head(sd_plot)
sd_plot$group=factor(gsub("ZT[0-9]+_","",gsub("_sd","",sd_plot$variable)),levels = c("NC_B6","HF_B6","NC_129","HF_129"))
sd_plot$time=factor(gsub("ZT","",gsub("_.*","",mean_plot$variable)),levels = c("2","6","10","14","18","22"))
head(sd_plot)

df2<-merge(mean_plot,sd_plot,by=c("Transcript_ID","Symbol","group","time"))
df2$col=factor(df2$group,levels = c("NC_B6","HF_B6","NC_129","HF_129"),labels = c("green","red","blue","darkred"))
Figure_tpm=list()
for (geid in 1:length(Trans_ID)){
  ge=Trans_ID[geid]
  dat_plot_gene=dat_plot[which(dat_plot$Transcript_ID == ge),]
  symbol=unique(dat_plot_gene$Symbol)
  df2_line=df2[which(df2$Transcript_ID == ge),]
  p1 <-ggplot() +
    geom_point(data=dat_plot_gene, aes(x = time, y = value,shape=group))+
    geom_line(data=df2_line,aes(x=time,y=mean,group=group,color=group))+ 
    geom_errorbar(data=df2_line,aes(x=time,ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05))+
    scale_x_discrete(limits = c("0","2","6","10","14","18","22"))+
    scale_color_manual(values = c("#008000",
                                  "#FF0000",
                                  "#0000FF",
                                  "#FFA500"))
  
  Figure_tpm[[geid]] <- p1 +theme_bw()+theme(legend.position="none")+
    labs(title = paste0(symbol," ",ge),x = "Time", y = "rpkm expressed value") +
    #labs(title = paste0(symbol," ",ge),x = "Time", y = "norm1e7 expressed value") +
    theme(axis.text.x=element_text(face="bold", size=10, angle = 45, vjust = 0.5),
          axis.text.y=element_text(face="bold", size=6),
          axis.title.y=element_text(size=8, face="bold"),
          strip.text.y = element_text(size = 8, face ="bold", angle = 0 ),
          legend.position="top")
}

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(data.table)
library(tibble)
library(grid)
arrange_pred=ggarrange(plotlist=Figure_tpm, nrow=6, ncol=6,  common.legend = TRUE)
arrange_pred=annotate_figure(arrange_pred, left = textGrob("quantie expresseion value", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                             bottom = textGrob("zt (hr)", gp = gpar(cex = 1.3)))
ggsave(paste0("interested_genes_quan.pdf"), arrange_pred, width = 15, height = 15)

```
