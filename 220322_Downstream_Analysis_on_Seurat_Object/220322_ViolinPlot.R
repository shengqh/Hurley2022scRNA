#ViolinPlot
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
sort_cell_type<-function(cts, sort_column){
  result<-cts[order(cts$seurat_clusters),]
  ct_uniq<-result[!duplicated(result[,sort_column]),]
  result[,sort_column]=factor(result[,sort_column], levels=ct_uniq[,sort_column])
  result<-result[order(result[,sort_column], result[,"seurat_clusters"]),]
  return(result)
}
gene="<gene of interest>"
clusters<-sort_cell_type(clusters, "seurat_clusters")
obj$seurat_cellactivity_clusters<-factor(obj$seurat_cellactivity_clusters, levels=clusters$seurat_cellactivity_clusters)
png(file=paste0(prefix, ".", gene, ".violin.1.png"), width=4000, height=2000, res=300)
#add split.by='group' to separate ICC/IDC-enriched vs benign-enriched
g=VlnPlot(obj, gene, group.by="seurat_cellactivity_clusters") + theme(axis.text.x = element_text(angle = 90)) + guides(colour = guide_legend(override.aes = list(size = 2), ncol=1)) + xlab("")
print(g)
dev.off()
