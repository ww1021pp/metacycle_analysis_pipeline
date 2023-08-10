rm(list=ls())

setwd("/workspace/rsrch2/panpanliu/project/RNA-seq/20_RNAseq.12W.HFD.EsrrgF.129Male/20_21_circadian_analysis/all_quatile/")

##read GFP B6

GFP_B6 <- read.delim("meta_out.GFP_B6.txt",stringsAsFactors = F)
rownames(GFP_B6)=GFP_B6$Transcript_ID
GFP_B6[, "GFPB6_max"] <- apply(GFP_B6[, 45:50], 1, max)
GFP_B6[, "GFPB6_min"] <- apply(GFP_B6[, 45:50], 1, min)
GFP_B6 <- transform(GFP_B6, GFPB6_Amp_C = (GFP_B6$GFPB6_max+0.01)/(GFP_B6$GFPB6_min+0.01))
GFP_B6_meta_above_2 <- GFP_B6[(GFP_B6$GFPB6_Amp_C>2)&(GFP_B6$meta2d_pvalue<0.01),];


##read GFP 129
GFP_129 <- read.delim("meta_out.GFP_129.txt",stringsAsFactors = F)
rownames(GFP_129)=GFP_129$Transcript_ID
GFP_129[, "GFP129_max"] <- apply(GFP_129[, 45:50], 1, max)
GFP_129[, "GFP129_min"] <- apply(GFP_129[, 45:50], 1, min)
GFP_129 <- transform(GFP_129, GFP129_Amp_C = (GFP_129$GFP129_max+0.01)/(GFP_129$GFP129_min+0.01))
GFP_129_meta_above_2 <- GFP_129[(GFP_129$GFP129_Amp_C>2)&(GFP_129$meta2d_pvalue<0.01),]

GFP129_spe <- setdiff(GFP_129_meta_above_2$Transcript_ID,GFP_B6_meta_above_2$Transcript_ID)

intersect_spe <- intersect(GFP129_spe,intersect(GFP_B6$Transcript_ID, GFP_129$Transcript_ID))
heatmap_data<-cbind(GFP_129[intersect_spe,45:50]/GFP_129[intersect_spe,"GFP129_max"],
                    GFP_B6[intersect_spe,45:50]/GFP_129[intersect_spe,"GFP129_max"])
colnames(heatmap_data,paste0("GFP_129.",colnames(GFP_129)[45:50]),
                      paste0("GFP_B6.",colnames(GFP_B6)[45:50]))

write.table(heatmap_data,"GFP_129_diff_B6_circadian_gene.txt", sep = "\t",
            quote = F,na="NA",row.names = F)

## sort the heat data by each time points and selected top n genes as following steps:
# step 1: sort the firt time pints with descending order and selected top n genes.
# step 2: sort the retained genes by second time point and  selected top genes.
# step 3: repeat the above steps untill sort all time points and all transcript
 after sort data by above steps, we can use pheatmap get a figure



heat_data <- read.delim("GFP_129_diff_B6_circadian_gene_sorted_byGFP129.txt",
                        stringsAsFactors = F,header = T)

library(pheatmap)
library(RColorBrewer)
library(gplots)

palette.breaks <- seq(0,1,0.02)
color.palette  <- colorRampPalette(c("blue","black", "yellow"))(length(palette.breaks) - 1)

pheatmap(heat_data,breaks =palette.breaks,color = color.palette,cluster_rows = F,cluster_cols = F,show_rownames = F)




