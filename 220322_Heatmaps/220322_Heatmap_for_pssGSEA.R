#Heatmap for pssGSEA
library(tidyverse)
library(pheatmap)
setwd("<folder where your tab delimited txt file is")
#"tab delimited txt file" contains consolidated pssGSEA Normalized Enrichment Scores (NES)
#Example txt file for epithelial cells (Figure 3) uploaded as "Example pssGSEA EPI.txt"
PJH_Heatmap <- read_tsv(file = "tab delimited txt file")
PJH_Heatmap <- as.data.frame(PJH_Heatmap)
PJH_Heatmap <- column_to_rownames(PJH_Heatmap, var = "Cluster")
p <- pheatmap(PJH_Heatmap)