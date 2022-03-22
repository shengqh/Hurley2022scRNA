#FeaturePlot
library(ggplot2)
finalList<-readRDS("<filepath to PH_scRNA.final.rds>")
obj<-finalList$obj
rm(finalList)
all_cells=colnames(obj)
obj$group="Tumor"
obj$group[grepl("Benign", all_cells)] = "Benign"
clusters<-read.csv("<filepath to PH_scRNA.rename_celltype.csv>", stringsAsFactors = F)
rownames(clusters)<-clusters$seurat_clusters
obj$cell_type=clusters[as.character(obj$seurat_clusters), "renamed_cellactivity_clusters"]
obj$seurat_cellactivity_clusters=clusters[as.character(obj$seurat_clusters), "seurat_renamed_cellactivity_clusters"]
#add split.by='group' to separate ICC/IDC-enriched vs benign-enriched
FeaturePlot(obj, features="<gene of interest>")