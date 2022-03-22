#rm(list=ls()) 
outFile='PH_scRNA'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct/result/PH_scRNA.final.rds'
parFile2='C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype_rename/result/PH_scRNA.rename_cluster.csv'
parFile3=''
pvalue=0.05;useRawPvalue=0;foldChange=1.5;DE_by_cell=1;filter_minTPM=1;filter_cellPercentage=0.2;bBetweenCluster=1;cluster_name='seurat_renamed_cellactivity_clusters'

setwd('C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype_rename_edgeR_betweenCluster_byCell/result')


library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(statmod)

sumcount<-function(ct_count, names, sample_df){
  result<-lapply(names, function(x){
    res<-ct_count[,sample_df$Cell[sample_df$Sample ==x],drop=F]
    apply(res, 1, sum)
  })
  rescount<-do.call(cbind, result)
  colnames(rescount)<-names
  return(rescount)
}

if(!exists('obj')){
  finalList<-readRDS(parFile1)
  obj<-finalList$obj
}

clusterDf<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
clusterDf<-clusterDf[colnames(obj),]
clusterDf$sample<-obj$orig.ident[rownames(clusterDf)]

obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]

sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
colnames(sampleGroups)<-c("Sample","Group")

comparisons<-data.frame(Cluster=c("12", "5", "12", "6", "12", "11", "21", "5", "21", "6", "21", "11"),
                                Group=c("Benign", "Tumor", "Benign", "Tumor", "Benign", "Tumor", "Benign", "Tumor", "Benign", "Tumor", "Benign", "Tumor"),
                                Comparison=c("Tumor_5_vs_Benign_12", "Tumor_5_vs_Benign_12", 
                                             "Tumor_6_vs_Benign_12", "Tumor_6_vs_Benign_12", 
                                             "Tumor_11_vs_Benign_12", "Tumor_11_vs_Benign_12", 
                                             "Tumor_5_vs_Benign_21", "Tumor_5_vs_Benign_21", 
                                             "Tumor_6_vs_Benign_21", "Tumor_6_vs_Benign_21", 
                                             "Tumor_11_vs_Benign_21", "Tumor_11_vs_Benign_21"))
comparisonNames<-unique(comparisons$Comparison)

comp <-comparisonNames[1]

designMatrix<-NULL
for (comp in comparisonNames){
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  
  cur_groups=unique(comp_groups$Group)
  
  controlGroup<-cur_groups[1]
  sampleGroup<-cur_groups[2]
  
  control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
  sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
  
  control_clusters<-comp_groups$Cluster[comp_groups$Group==controlGroup]
  sample_clusters<-comp_groups$Cluster[comp_groups$Group==sampleGroup]
  
  control_name=paste0(controlGroup, ".", paste(control_clusters, collapse = "."))
  sample_name=paste0(sampleGroup, ".", paste(sample_clusters, collapse = "."))
  
  prefix<-paste0(outFile, ".", comp, ".edgeR")
  
  control_cells<-rownames(clusterDf)[clusterDf$seurat_clusters %in% control_clusters & clusterDf$sample %in% control_names]  
  sample_cells<-rownames(clusterDf)[clusterDf$seurat_clusters %in% sample_clusters & clusterDf$sample %in% sample_names]  

  all_cells<-c(control_cells, sample_cells)

  de_obj<-subset(obj, cells=all_cells)
  de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
  de_obj$DisplayGroup<-c(rep(control_name, length(control_cells)), rep(sample_name, length(sample_cells)))
    
  designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
  rm(de_obj)
  
  designdata$Sample=gsub("_Benign|_Tumor","",designdata$Sample)
  
  designfile<-paste0(prefix, ".design")
  write.csv(designdata, file=designfile, row.names=F, quote=F)
  
  curdf<-data.frame(prefix=prefix, cellType="", comparison=comp, sample_in_formula=1, design=designfile, stringsAsFactors = F)
  if (is.null(designMatrix)){
    designMatrix = curdf
  }else{
    designMatrix = rbind(designMatrix, curdf)
  }
}

result<-NULL
idx<-1
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  cat(prefix, "\n")
  
  designfile=designMatrix[idx, "design"]
  cellType=designMatrix[idx, "cellType"]
  comp=designMatrix[idx, "comparison"]
  sample_in_formula=designMatrix[idx, "sample_in_formula"]
  
  dge_filename <-paste0(prefix, ".csv")
  sigFile<-paste0(prefix, ".sig.csv")
  sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=bBetweenCluster, "sample_in_formula"=sample_in_formula, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
  
  #if(file.exists(gseaFile)){
  #  next
  #}
  
  designdata<-read.csv(designfile, stringsAsFactors = F)
  groups<-designdata$Group
  
  de_obj<-subset(obj, cells=designdata$Cell)
  cells<-as.matrix(de_obj[["RNA"]]@counts)
    
  #filter genes with zero count
  cells<-cells[rowSums(cells)>0,]

  #filter genes by tpm
  tpm = sweep(cells, 2, colSums(cells)/1e6, "/")
  min_sample<-filter_cellPercentage * ncol(cells)
  valid_genes=rowSums(tpm > filter_minTPM)
  keep_rows <- valid_genes >= min_sample
  rm(tpm)
  
  cat("get GSEA rank file using most detected genes\n")
  if(length(valid_genes) > 10000){
    valid_genes<-valid_genes[order(valid_genes, decreasing = T)]
    gsea_cells<-cells[names(valid_genes[1:10000]),]
  }else{
    gsea_cells<-cells
  }
  
  cdr <- scale(colMeans(gsea_cells > 0))
  
  if(sample_in_formula){
    samples<-designdata$Sample
    design <- model.matrix(~ cdr + samples + groups)
    gs<-table(designdata$Group, designdata$Sample)
    notpaired_samples<-colnames(gs)[colSums(gs> 0)==1]
    design<-design[, !(colnames(design) %in% paste0("samples", notpaired_samples))]
  }else{
    design <- model.matrix(~ cdr + groups)
  }
  
  rownames(design)<-colnames(gsea_cells)
  dge<-DGEList(gsea_cells, group=groups)
  cat("  calcNormFactors", "\n")
  dge<-calcNormFactors(dge)
  
  cat("  estimateDisp", "\n")
  dge<-estimateDisp(dge,design=design)
  
  cat("  glmQLFit", "\n")
  fitqlf<-glmQLFit(dge,design=design,robust=TRUE)
  qlf<-glmQLFTest(fitqlf)
  out<-topTags(qlf, n=Inf)
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * out$table$F)
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  cat("now perform DE with filtered genes, filter genes by tpm\n")
  
  cells<-cells[keep_rows,]
  cdr <- scale(colMeans(cells > 0))
  
  if(sample_in_formula){
    samples<-designdata$Sample
    design <- model.matrix(~ cdr + samples + groups)
    gs<-table(designdata$Group, designdata$Sample)
    notpaired_samples<-colnames(gs)[colSums(gs> 0)==1]
    design<-design[, !(colnames(design) %in% paste0("samples", notpaired_samples))]
  }else{
    design <- model.matrix(~ cdr + groups)
  }
  
  rownames(design)<-colnames(cells)
  write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
  
  dge<-DGEList(cells, group=groups)
  cat("  calcNormFactors", "\n")
  dge<-calcNormFactors(dge)
  
  cat("  estimateDisp", "\n")
  dge<-estimateDisp(dge,design=design)
  
  cat("  glmQLFit", "\n")
  fitqlf<-glmQLFit(dge,design=design,robust=TRUE)
  qlf<-glmQLFTest(fitqlf)
  out<-topTags(qlf, n=Inf)
  dge_filename <-paste0(prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  sigFile<-paste0(prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)
  
  siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)
