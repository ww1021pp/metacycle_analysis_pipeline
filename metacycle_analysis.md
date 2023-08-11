# metacycle analysis to predict circadian geness pipiline for time-series data with R

### metacycle analysis for time series expression data to predict circadian genes/transcripts
### This is a script used to detect circadian genes from quantile expression value got from Homer analysisRepeat.pl
![image](https://github.com/ww1021pp/metacycle_analysis_pipeline/assets/60449311/97627478-57fa-4469-9569-41fdf68048f8)

![image](https://github.com/ww1021pp/metacycle_analysis_pipeline/assets/60449311/754f0f21-48a5-44ce-8287-20c935606150)


```
rm(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/RNA-seq/20_RNAseq.12W.HFD.EsrrgF.129Male/circadian_analysis") ## set work directory

library(MetaCycle)            ## import libarary 
GFP_129_raw<-read.delim("../4_homer/ErrgGFP_129_quatile.txt")  ## import the counts table 
Cre_129_raw<-read.delim("../4_homer/ErrgCre_129_quatile.txt")
colnames(GFP_129_raw)[1]="Transcript_ID"
colnames(Cre_129_raw)[1]="Transcript_ID"   ## names each column with our own defination

rownames(GFP_129_raw)=GFP_129_raw$Transcript_ID   ## name each row with transcript ID to extract the interested gene 

rownames(Cre_129_raw)=Cre_129_raw$Transcript_ID ## names each row with transcript/gene ID
count_GFP_129=GFP_129_raw[,9:26]         ## extract count table

raw<-merge(Cre_129_raw,count_GFP_129,by="row.names")[,-1]

colnames(raw)=gsub("\\.reads.*","",colnames(raw))  ## substitue colnames of count table
colnames(raw)[1]="Transcript_ID"    
tags<-c("ErrgCreHFD\\.129","ErrgGFPHFD\\.129")    ## the tags we are interested
time_points<-as.numeric(rep(c("2","6","10","14","18","22"),each=3))  ### define the time points and marked them as numberic
raw$Symbol<-sapply(strsplit(as.character(raw$Annotation.Divergence), "\\|"), "[", 1)  ## define gene symbol for each transcript
raw=raw[which(raw$Annotation.Divergence != "0.000"),] ## filter the transcripts without gene name
head(raw)
dim(raw)
write.table(raw,file="../4_homer/merged_Cre_GFP_129_quatile.txt",row.names=F,quote = F,sep = "\t") ## write out the file we are interested

raw$ZT2_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT2\\."),colnames(raw))])  ## get the row of each condition for each time point
raw$ZT6_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT6\\."),colnames(raw))]) 
raw$ZT10_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT10\\."),colnames(raw))])
raw$ZT14_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT14\\."),colnames(raw))])
raw$ZT18_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT18\\."),colnames(raw))])
raw$ZT22_GFP <- rowMeans(raw[,grepl(pattern = paste0("ErrgGFPHFD\\.129","\\.ZT22\\."),colnames(raw))])

raw$ZT2_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT2\\."),colnames(raw))])
raw$ZT6_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT6\\."),colnames(raw))])
raw$ZT10_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT10\\."),colnames(raw))])
raw$ZT14_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT14\\."),colnames(raw))])
raw$ZT18_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT18\\."),colnames(raw))])
raw$ZT22_Cre <- rowMeans(raw[,grepl(pattern = paste0("ErrgCreHFD\\.129","\\.ZT22\\."),colnames(raw))])

head(raw)
dim(raw)



#####fileter low expression genes by condition avarage expression gt>1, AT LEAST  half condition expression >filter_value#####
raw[,"count"] <- rowSums(raw[,46:57]>5);head(raw);dim(raw) 
raw <- raw[raw$count>6,];dim(raw)
#raw[, "max"] <- apply(raw[, 46:57], 1, max)
#raw[, "min"] <- apply(raw[, 46:57], 1, min)
#a<-as.matrix(raw[,46:57]/raw$max)
#raw_uniq_tmp<-raw[order(raw[,9],raw[,6],decreasing=TRUE),]
#raw_uniq <- raw_uniq_tmp[match(unique(raw_uniq_tmp$Symbol), raw_uniq_tmp$Symbol),]
raw$annot_cobind=paste0(raw$Transcript_ID,"##",raw$Symbol)

for(con in tags){
  exp_data=cbind(raw[,"annot_cobind"],
                 raw[,grepl(pattern = paste0(con,"\\.ZT2\\."),colnames(raw))],
                 raw[,grepl(pattern = paste0(con,"\\.ZT6\\."),colnames(raw))],
                 raw[,grepl(pattern = paste0(con,"\\.ZT10\\."),colnames(raw))],
                 raw[,grepl(pattern = paste0(con,"\\.ZT14\\."),colnames(raw))],
                 raw[,grepl(pattern = paste0(con,"\\.ZT18\\."),colnames(raw))],
                 raw[,grepl(pattern = paste0(con,"\\.ZT22\\."),colnames(raw))]
                 )  ## cbind expression table for each condition for the following metacycle analysis
  colnames(exp_data)[1]="Transcript_ID" 
  tmp=gsub(".*ZT","",colnames(exp_data)[-1])  
  time_points=as.numeric(gsub("\\..*","",tmp)) 
  write.csv(exp_data, file=paste0(con,".csv"), row.names=FALSE)
  meta2d(infile=paste0(con,".csv"), filestyle="csv",
         outdir=paste0(con,"_metacycle_res"), timepoints=time_points,outSymbol=con)  ## outdir is the results outdirectory
 }

options(stringsAsFactors = FALSE)
Cre_129=read.csv("./ErrgCreHFD\\.129_metacycle_res/ErrgCreHFD\\.129meta2d_.129.csv",header=TRUE);head(Cre_129);dim(Cre_129)
Cre_129_Exp=read.csv("./ErrgCreHFD\\.129.csv")
head(Cre_129_Exp);dim(Cre_129_Exp)
Cre_129_results=cbind(raw[,1:8],raw$Symbol,Cre_129[,-1],Cre_129_Exp[,-1],raw[,grepl("ZT[0-9]+_Cre",colnames(raw))])
save(Cre_129_results,file=paste("meta_out","Cre_129","rda",sep="."))
write.table(Cre_129_results,file=paste("meta_out","Cre_129","txt",sep="."),row.names=F,quote = F,sep = "\t",na="")

GFP_129=read.csv("./ErrgGFPHFD\\.129_metacycle_res/ErrgGFPHFD\\.129meta2d_.129.csv",header=TRUE)
GFP_129_Exp=read.csv("./ErrgGFPHFD\\.129.csv")
GFP_129_results=cbind(raw[,1:8],raw$Symbol,GFP_129[,-1],GFP_129_Exp[,-1],raw[,grepl("ZT[0-9]+_GFP",colnames(raw))])
save(GFP_129_results,file=paste("meta_out","GFP_129","rda",sep="."))
write.table(GFP_129_results,file=paste("meta_out","GFP_129","txt",sep="."),row.names=F,quote=F,sep = "\t",na="")

```
