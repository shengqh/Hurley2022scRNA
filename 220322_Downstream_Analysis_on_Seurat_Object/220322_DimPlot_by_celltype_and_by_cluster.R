#DimPlot (by cell type and by cluster)
library(Seurat)
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
prefix="PH_scRNA"
#DimPlot by cell type
#add split.by='group' to separate ICC/IDC-enriched vs benign-enriched
png(file=paste0(prefix, ".celltype.png"), width=3300, height=3000, res=300)
g=DimPlot(obj, group.by="cell_type", label=T) + guides(colour = guide_legend(override.aes = list(size = 6), ncol=1))
print(g)
dev.off()
#DimPlot by clusters
#add split.by='group' to separate ICC/IDC-enriched vs benign-enriched
png(file=paste0(prefix, ".seurat_cellactivity_clusters.png"), width=3300, height=3000, res=300)
g=DimPlot(obj, group.by="seurat_cellactivity_clusters", label=T) + guides(colour = guide_legend(override.aes = list(size = 6), ncol=1))
print(g)
dev.off()