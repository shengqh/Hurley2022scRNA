rm(list=ls()) 
outFile='PH_scRNA'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct/result/PH_scRNA.final.rds'
parFile2='/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype_rename/result/PH_scRNA.rename_cluster.csv'
parFile3='/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/hto_clonotype_4_convert/result/clonotypes.csv'
outputDirectory='.'


setwd('/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/hto_clonotype_5_vis/result')

source("scRNA_func.r")
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

valid_obj[["final_seurat_clusters"]]=valid_cell_df$seurat_cellactivity_clusters

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
  sg<-DimPlot(valid_obj, reduction = "umap", cells.highlight = samplelist, cols.highlight=gg_color_hue(length(samplelist)) , label=F) + theme(legend.position = "top")
  
  pdf(file=paste0(curname, ".umap.pdf"), width=21, height=7)
  p<-ggarrange(gcell, cg, sg, ncol=3)
  print(p)
  dev.off()
}
