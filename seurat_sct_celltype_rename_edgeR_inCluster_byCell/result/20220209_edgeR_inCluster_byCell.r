#rm(list=ls()) 
outFile='PH_scRNA'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='c:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct/result/PH_scRNA.final.rds'
parFile2='c:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype_rename/result/PH_scRNA.rename_cluster.csv'
parFile3=''
pvalue=0.05;useRawPvalue=0;foldChange=1.5;DE_by_cell=1;filter_minTPM=1;filter_cellPercentage=0.2;bBetweenCluster=0;cluster_name='seurat_renamed_cellactivity_clusters'

setwd('c:/projects/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/seurat_sct_celltype_rename_edgeR_inCluster_byCell/result')


library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)

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
obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]

sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
colnames(sampleGroups)<-c("Sample","Group")

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
colnames(comparisons)<-c("Group", "Comparison")
comparisonNames<-unique(comparisons$Comparison)

comp <-comparisonNames[1]

designMatrix<-NULL
for (comp in comparisonNames){
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  
  controlGroup<-comp_groups$Group[1]
  sampleGroup<-comp_groups$Group[2]
  
  control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
  sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
  
  if(bBetweenCluster){
    prefix<-paste0(outFile, ".", comp, ".edgeR")
    
    control_cells<-rownames(clusterDf)[clusterDf[,cluster_name] %in% control_names]  
    sample_cells<-rownames(clusterDf)[clusterDf[,cluster_name] %in% sample_names]  
  
    all_cells<-c(control_cells, sample_cells)

    de_obj<-subset(obj, cells=all_cells)
    de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
    de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
      
    designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
    designfile<-paste0(prefix, ".design")
    write.csv(designdata, file=designfile, row.names=F, quote=F)
    
    curdf<-data.frame(prefix=prefix, cellType="", comparison=comp, sampleInGroup=1, design=designfile, stringsAsFactors = F)
    if (is.null(designMatrix)){
      designMatrix = curdf
    }else{
      designMatrix = rbind(designMatrix, curdf)
    }
  }else{
    cts = unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name])
    prefixList<-gsub(" ", "_", cts)
    prefixList<-gsub(":", "_", prefixList)
    prefixList<-gsub("_+", "_", prefixList)
    
    idx<-1
    for (idx in c(1:length(cts))){
      ct = cts[idx]
      prefix = paste0(outFile, ".", prefixList[idx], ".", comp, ".edgeR")
      
      clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
      de_obj<-subset(obj, cells=rownames(clusterCt))
      clusterCt$sample=de_obj$orig.ident

      invalid_control_names= control_names[!(control_names %in% unique(clusterCt$sample))]
      invalid_sample_names= sample_names[!(sample_names %in% unique(clusterCt$sample))]

      if (length(invalid_control_names) == length(control_names)){
        stop(paste0("There were no control ", paste0(invalid_control_names, collapse=","), " found in object sample names!"))
      }
      
      if (length(invalid_sample_names)  == length(sample_names)){
        stop(paste0("There were no sample ", paste0(invalid_sample_names, collapse=","), " found in object sample names!"))
      }
      
      control_cells<-rownames(clusterCt)[clusterCt$sample %in% control_names]  
      sample_cells<-rownames(clusterCt)[clusterCt$sample %in% sample_names]  
      
      all_cells<-c(control_cells, sample_cells)
      
      de_obj<-subset(obj, cells=all_cells)
      de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
      de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
      
      designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
      designdata$Sample<-gsub("_Benign|_Tumor", "", designdata$Sample)
      designfile<-paste0(prefix, ".design")
      write.csv(designdata, file=designfile, row.names=F, quote=F)
      
      curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, sampleInGroup=1, design=designfile, stringsAsFactors = F)
      if (is.null(designMatrix)){
        designMatrix = curdf
      }else{
        designMatrix = rbind(designMatrix, curdf)
      }
    }
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
  sampleInGroup=designMatrix[idx, "sampleInGroup"]
  
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
  dge_filename <-paste0(prefix, ".csv")
  sigFile<-paste0(prefix, ".sig.csv")
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=bBetweenCluster, "sampleInGroup"=sampleInGroup, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }

  if(file.exists(sigGenenameFile)){
    next
  }

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
  
  if(sampleInGroup){
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
  
  if(sampleInGroup){
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
  write.csv(out$table, file=dge_filename, quote=F)

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  write.csv(sigout, file=sigFile, quote=F)
  
  siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)
