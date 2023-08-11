# metacycle analysis results summary and plot venny plots for different condition to see the overlap of catogary

## metacycle statistics from meta2d out files

metacycle_result_statics.R
### this scripts use the files of meta2d() function, and merge avaerage expression for each time points and get the peak_through (max_mean/mean_max) for one condition in different time series.
```
rm(list=ls())
setwd("~/Documents/DY_RNA_Seq_4_8_19/4_B6_129_NC_HFD_RNAseq-206142945/circadian_analysis1")
NC_129_meta <- read.delim("meta_out.NC_129.txt",header = T)
NC_129_meta<-NC_129_meta[which(NC_129_meta$raw_uniq.Symbol !="0.000"),]
head(NC_129_meta);dim(NC_129_meta)
NC_129_meta[, "max"] <- apply(NC_129_meta[, 45:50], 1, max)
NC_129_meta[, "min"] <- apply(NC_129_meta[, 45:50], 1, min)
NC_129_meta <- transform(NC_129_meta, Amp_C = (NC_129_meta$max+0.01)/(NC_129_meta$min+0.01))
NC_129_meta_above_1.5 <- NC_129_meta[(NC_129_meta$Amp_C>1.5)&(NC_129_meta$meta2d_pvalue<0.01),];head(NC_129_meta);dim(NC_129_meta_above_1.5)
write.table(NC_129_meta_above_1.5,file="NC_129_cir3592.txt",row.names=F,col.names=T, na="",sep="\t")

#############################################3
###############GFP_129###############
#########################
HF_129_meta <- read.delim("meta_out.HF_129.txt",header = T)
HF_129_meta=HF_129_meta[which(HF_129_meta$raw_uniq.Symbol != "0.000"),]
head(HF_129_meta);dim(HF_129_meta)
HF_129_meta[, "max"] <- apply(HF_129_meta[, 45:50], 1, max)
HF_129_meta[, "min"] <- apply(HF_129_meta[, 45:50], 1, min)
HF_129_meta <- transform(HF_129_meta, Amp_C = (HF_129_meta$max+0.01)/(HF_129_meta$min+0.01))
HF_129_meta_above_1.5 <- HF_129_meta[(HF_129_meta$Amp_C>1.5)&(HF_129_meta$meta2d_pvalue<0.01),];head(HF_129_meta_above_1.5);dim(HF_129_meta_above_1.5)
write.table(HF_129_meta_above_1.5,file="HF_129_cir2470.txt",row.names=F,col.names=T, na="",sep="\t")
###########Cre_B6############################
###############################33
#############################################
NC_B6_meta <- read.delim("meta_out.NC_B6.txt",header = T)
NC_B6_meta=NC_B6_meta[which(NC_B6_meta$raw_uniq.Symbol != "0.000"),]
head(NC_B6_meta);dim(NC_B6_meta)
NC_B6_meta[, "max"] <- apply(NC_B6_meta[, 45:50], 1, max)
NC_B6_meta[, "min"] <- apply(NC_B6_meta[, 45:50], 1, min)
NC_B6_meta <- transform(NC_B6_meta, Amp_C = (NC_B6_meta$max+0.01)/(NC_B6_meta$min+0.01))

NC_B6_meta_above_1.5 <- NC_B6_meta[(NC_B6_meta$Amp_C>1.5)&(NC_B6_meta$meta2d_pvalue<0.01),];head(NC_B6_meta_above_1.5);dim(NC_B6_meta_above_1.5)
write.table(NC_B6_meta_above_1.5,file="NC_B6_cir971.txt",row.names=F,col.names=T, na="",sep="\t")

############################################################
###############################GFP_B6#####################
############################################################
HF_B6_meta <- read.delim("meta_out.HF_B6.txt",header = T)
HF_B6_meta=HF_B6_meta[which(HF_B6_meta$raw_uniq.Symbol != ""),]
head(HF_B6_meta);dim(HF_B6_meta)
HF_B6_meta[, "max"] <- apply(HF_B6_meta[, 45:50], 1, max)
HF_B6_meta[, "min"] <- apply(HF_B6_meta[, 45:50], 1, min)
HF_B6_meta <- transform(HF_B6_meta, Amp_C = (HF_B6_meta$max+0.01)/(HF_B6_meta$min+0.01))
HF_B6_meta_above_1.5 <- HF_B6_meta[(HF_B6_meta$Amp_C>1.5)&(HF_B6_meta$meta2d_pvalue<0.01),];head(HF_B6_meta_above_1.5);dim(HF_B6_meta_above_1.5)
write.table(HF_B6_meta_above_1.5,file="HF_B6_cir741.txt",row.names=F,col.names=T, na="",sep="\t")

```


## venny plot for 4 groups 

```
library(VennDiagram);
venn.diagram(
  x = list(
    NC.B6= NC_B6_Transcript,  ## the list for each group we are interested
    HF.129 = HF_129_Transcript, 
    HF.B6 = HF_B6_Transcript,
    NC.129 = NC_129_Transcript
  ),
  filename = "NC_HFD_B6_129_FC1.5_group.tiff",
  col = "black",
  lty = "solid",
  lwd = 4,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  #cex=0, ##for num label
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.1,
  ##cat.cex=0, ###for catogory
  cat.fontfamily = "serif"
)
```
![image](https://github.com/ww1021pp/metacycle_analysis_pipeline/assets/60449311/7138ef9e-39ac-48ca-adb9-1b6782af1f3f)

