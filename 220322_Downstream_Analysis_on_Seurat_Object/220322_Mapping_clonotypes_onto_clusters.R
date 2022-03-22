#Mapping clonotypes onto clusters
rm(list=ls()) 
outFile='PH_scRNA'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parSampleFile4=''
parFile1='<filepath to PH_scRNA.final.rds>'
parFile2='<filepath to PH_scRNA.cluster.csv>'
parFile3='<filepath to clonotypes.csv>'
outputDirectory='.'
setwd('<folder you want your output files to be saved in>')

read_cell_cluster_file<-function(fileName, sort_cluster_name="seurat_clusters"){
  result<-read.csv(fileName, stringsAsFactors = F, row.names = 1)
  
  display_sort_cluster_name = paste0("display_", sort_cluster_name)
  result[,display_sort_cluster_name] = paste0("Cluster ", result[,sort_cluster_name])
  
  cluster_names=colnames(result)[grepl("_clusters", colnames(result))]
  
  sort_clusters_num = length(unique(result[,sort_cluster_name]))
  for(cluster_name in cluster_names){
    cluster_num = length(unique(result[,cluster_name]))
    if(cluster_name == sort_cluster_name){
      next
    }
    
    if (cluster_num != sort_clusters_num) {
      next
    }
    
    cf<-unique(result[,c(sort_cluster_name, cluster_name)])
    if(nrow(cf) != sort_clusters_num){
      next
    }
    
    cf<-cf[order(as.numeric(cf[,sort_cluster_name]), decreasing = T),]
    cf_levels=cf[,cluster_name]
    result[,cluster_name] = factor(result[,cluster_name], levels=cf_levels)
  }
  return(result)
}

library(Seurat)
library(ggplot2)
library(ggpubr)
finalList<-readRDS(parFile1)
obj<-finalList$obj
cell_df<-read_cell_cluster_file(parFile2)
clonotypes<-read.csv(parFile3)
clonotypes<-clonotypes[c(1:min(10, nrow(clonotypes))),,drop=F]
clonocells<-unlist(strsplit(clonotypes$cells, split = ";"))
cell_df_clono<-subset(cell_df, rownames(cell_df) %in% clonocells)
cell_df_clono<-cell_df_clono[order(cell_df_clono$seurat_clusters),]
clono_clusters<-unique(cell_df_clono$seurat_clusters)
display_clusters<-as.character(unique(cell_df_clono$seurat_cellactivity_clusters))
valid_cell_df<-subset(cell_df, cell_df$seurat_clusters %in% clono_clusters)
valid_obj=subset(obj, cells=rownames(valid_cell_df))
rm(obj)
valid_obj[["final_seurat_clusters"]]=valid_cell_df$seurat_clusters
gcell<-DimPlot(valid_obj, group.by="final_seurat_clusters", reduction = "umap", label=T) + theme(legend.position = "none") + ggtitle("")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

idx=1
for(idx in c(1:nrow(clonotypes))) {
  curname=clonotypes$clonotype_id[idx]
  cdr3s_aa=clonotypes$cdr3s_aa[idx]
  cdr3s_cells=strsplit(clonotypes$cells[idx], split = ";")
  names(cdr3s_cells)<-cdr3s_aa
  cg<-DimPlot(valid_obj, reduction = "umap", cells.highlight = cdr3s_cells, label=F) + theme(legend.position = "none") + 
    scale_color_manual(values = c("lightgrey", "red")) + ggtitle(cdr3s_aa)
  
  cobj=subset(valid_obj, cells=unlist(cdr3s_cells))
  sc<-data.frame(cell=colnames(cobj), sample=cobj$orig.ident)
  samplelist<-tapply(sc$cell,sc$sample,list)
  sg<-DimPlot(valid_obj, reduction = "umap", cells.highlight = samplelist, cols.highlight=gg_color_hue(length(samplelist)) , label=F) + theme(legend.position = "top") + scale_color_manual(values = c("lightgrey", '#FF0000', '#0070C0'))
  
  pdf(file=paste0(curname, ".umap.pdf"), width=21, height=7)
  p<-ggarrange(gcell, cg, sg, ncol=3)
  print(p)
  dev.off()
}
