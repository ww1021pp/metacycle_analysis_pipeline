# venny plot for 4 groups 


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
