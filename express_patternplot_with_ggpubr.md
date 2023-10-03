```
rm(list=ls())
wd="/workspace/rsrch2/panpanliu/project/RNA-seq/RNA-seq_082523/M_Jh_circadian_analysis"  
setwd(wd)

options(stringsAsFactors = FALSE)
library(data.table)
library(MetaCycle)
library(dplyr)
library(ggplot2)
library(ggpubr)



data_raw<-read.delim("../4_homer/MZT_Jh_quatile.txt", header = TRUE,stringsAsFactors = FALSE) ## read expression data ###
colnames(data_raw)[1]="Transcript_ID"

rownames(data_raw)=data_raw$Transcript_ID
colnames(data_raw)<-sapply(colnames(data_raw),function(x) gsub(".*M_ZT","M_ZT",as.character(x)))
colnames(data_raw)<-sapply(colnames(data_raw),function(x) gsub("\\.reads.*","",as.character(x)))
data_raw$Symbol<-sapply(strsplit(data_raw$Annotation.Divergence,"\\|"),getElement, 1)


#######################remove the transcripts without gene annotation(As we focused on the known protein)#####
data_raw<-data_raw %>% filter(Symbol != "0.000") 
print(dim(data_raw))


cir_genes <- read.delim("./meta_out.M_rawcountsgt3.txt",header = T,stringsAsFactors = F) %>%   #### filter circadian genes 
  filter(meta2d_pvalue < 0.01 & Amp_C> 2 ) 

cir_uniq_tmp<-cir_genes[order(cir_genes[,"Symbol"],cir_genes[,"Length"],decreasing=TRUE),]
cir_uniq <- cir_uniq_tmp[match(unique(cir_uniq_tmp$Symbol), cir_uniq_tmp$Symbol),]

## convert long dataframe to wide dataframe with library(reshape2)

library(reshape2)


data_plot <- data_raw %>% filter(Transcript_ID %in% cir_uniq$Transcript_ID) %>% select(c(25,9:24)) 
rownames(data_plot)=NULL

library(tidyverse)
data_plot <- column_to_rownames(data_plot, var="Symbol")
data_plot <- t(data_plot)
```
![image](https://github.com/ww1021pp/metacycle_analysis_pipeline/assets/60449311/7265acfc-3989-4608-9cd6-01b25a9b7cfe)



```
data_plot <- tibble::rownames_to_column(as.data.frame(data_plot), "Sample") 

exp_data<- data_plot %>% mutate(condition = as.factor(sapply(strsplit(Sample,"_"),getElement, 1)),
                                  time = factor(sapply(strsplit(Sample,"_"),getElement, 2),levels = c("ZT4","ZT10","ZT16","ZT22")),
                                  drugTime = sapply(strsplit(Sample,"_"),getElement, 3))
core_clock<-c("ESRRG","ARNTL","NPAS2","CRY1","NFIL3","CLOCK","RORC","NR1D1","CRY2","BHLHE41","PER2","DBP","TEF","PER1","NR1D2","PER3")
core_clock <- unlist(sapply(core_clock,str_to_title))

genes<-intersect(core_clock,cir_genes$Symbol)

ggsave("coreclock_gene_expression.pdf",width = 8,)
figure_list <-list()
for (i in 1:length(genes)){
figure_list[[i]] <- ggline(exp_data, x = "time",
       y = genes[i],
       combine = TRUE,
       color = "Black",
       add = c("mean_sd","point"),
       add.params = list(color = "grey", size = 0.5)
            )
figure_list[[i]]<-ggpar(figure_list[[i]],
      main = genes[i],
      xlab = FALSE,
      ylab = FALSE,
      font.main=c(14,"bold"),
      font.xtickslab=c(14,"plain","black"), 
      font.ytickslab=c(14,"plain","black"),
     # font.x = c(14,"bold"),
     #  font.y = c(14,"bold")
  #    xlab = "ZT time", ylab = "quantie expression"
  )
}

library(grid)
arrange_pred=ggarrange(plotlist=figure_list, nrow=3, ncol=4,  common.legend = TRUE)
arrange_pred=annotate_figure(arrange_pred, left = textGrob("quantie expresseion value", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                             bottom = textGrob("zt (hr)", gp = gpar(cex = 1.3)))
ggsave(paste0("interested_genes_quan.pdf"), arrange_pred, width = 15, height = 10)
```

![image](https://github.com/ww1021pp/metacycle_analysis_pipeline/assets/60449311/2d159e4f-0773-43a5-a707-8d02d79221a1)

