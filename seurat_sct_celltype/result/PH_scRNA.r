rm(list=ls()) 
outFile='PH_scRNA'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct/result/PH_scRNA.cluster.normByUpQuantile.csv'
parFile2='/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct/result/PH_scRNA.cluster.csv'
parFile3=''


setwd('/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype/result')

source("scRNA_func.r")

library(data.table)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file

remove_subtype_of=ifelse(annotate_tcell, "T cells", "")
remove_subtype_of=paste0(remove_subtype_of, ";", remove_subtype)

data.norm=read.csv(parFile1, row.names=1,check.names = F)

cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file)
if("curated_markers_file" %in% names(myoptions)){
  curated_markerfile<-myoptions$curated_markers_file
  curated_markers_df<-read.table(curated_markerfile, sep="\t", header=F, stringsAsFactors=F)
  curated_markers_celltype<-split(curated_markers_df$V2, curated_markers_df$V1)
  cellType=cell_activity_database$cellType
  for(cmct in names(curated_markers_celltype)){
    cellType[[cmct]]=curated_markers_celltype[[cmct]]
  }
  weight=calc_weight(cellType)
  cell_activity_database=list(cellType=cellType, weight=weight)
}

predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)

cta_ora_mat = get_cta_ora_mat(predict_celltype)
cell_activity_database$predicted<-predict_celltype
cell_activity_database$cta_mat=cta_ora_mat$cta_mat
cell_activity_database$ora_mat=cta_ora_mat$ora_mat

cell_type=list(cell_activity_database=cell_activity_database)

new.cluster.ids<-names(predict_celltype$max_cta)
names(new.cluster.ids) <- colnames(data.norm)

id_tbl=data.frame("seurat_clusters"=names(new.cluster.ids), "cell_type"=new.cluster.ids)
rownames(id_tbl)=id_tbl$seurat_clusters

if (annotate_tcell){
  tcell_activity_database<-get_selfdefinedCelltype(tcell_markers_file)
  tcell_clusters<-names(new.cluster.ids)[new.cluster.ids=="T cells"]
  tcell_data.norm<-data.norm[,tcell_clusters,drop=F]
  tcell_predict_celltype<-ORA_celltype(tcell_data.norm,tcell_activity_database$cellType,tcell_activity_database$weight)

  cta_ora_mat = get_cta_ora_mat(tcell_predict_celltype)
  tcell_activity_database$predicted<-tcell_predict_celltype
  tcell_activity_database$cta_mat=cta_ora_mat$cta_mat
  tcell_activity_database$ora_mat=cta_ora_mat$ora_mat
  
  tcell_new.cluster.ids<-names(tcell_predict_celltype$max_cta)
  names(tcell_new.cluster.ids) <- tcell_clusters
  
  id_tbl$tcell_type=""
  id_tbl[tcell_clusters,"tcell_type"]=tcell_new.cluster.ids
  
  new.cluster.ids[tcell_clusters] = tcell_new.cluster.ids
  
  cell_type$tcell_activity_database=tcell_activity_database
}

id_tbl$cellactivity_clusters=new.cluster.ids
write.csv(id_tbl, file=paste0(outFile, ".celltype.csv"), row.names=F)
saveRDS(cell_type, file=paste0(outFile, ".celltype.rds"))

clusters=read.csv(parFile2, header=T, stringsAsFactors = F, row.names=1)
clusters$cellactivity_clusters=new.cluster.ids[as.character(clusters$seurat_clusters)]
clusters$seurat_cellactivity_clusters=paste0(clusters$seurat_clusters, " : ", clusters$cellactivity_clusters)
write.csv(clusters, file=paste0(outFile, ".celltype_cluster.csv"))
