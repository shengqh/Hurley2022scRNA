#Pseudotime Analysis on T Cells with Monocle
library(Seurat)
library(SeuratWrappers)
finalList<-readRDS("<file path to PH_scRNA.final.rds>")
obj<-finalList$obj
rm(finalList)
clusters<-read.csv("<filepath to PH_scRNA.rename_celltype.csv>", stringsAsFactors = F)
rownames(clusters)<-clusters$seurat_clusters
obj$cell_type=clusters[as.character(obj$seurat_clusters), "renamed_cellactivity_clusters"]
obj$seurat_cellactivity_clusters=clusters[as.character(obj$seurat_clusters), "seurat_renamed_cellactivity_clusters"]
all_cells=colnames(obj)
tcell_clusters=subset(clusters, clusters$renamed_cellactivity_clusters == "T cells")
tcell_cells=all_cells[obj$seurat_clusters %in% tcell_clusters$seurat_clusters]
tcell_obj=subset(obj, cells=tcell_cells)
rm(obj)
tcell.cds<-as.cell_data_set(tcell_obj)
plot_cells(tcell.cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")
tcell.cds <- cluster_cells(tcell.cds, reduction_method = c("UMAP"))
plot_cells(tcell.cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")
tcell.cds <- learn_graph(tcell.cds, use_partition = TRUE)
plot_cells(tcell.cds, label_groups_by_cluster=FALSE,  color_cells_by = "partition", graph_label_size=3)
cds <- order_cells(tcell.cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           show_trajectory_graph = FALSE)