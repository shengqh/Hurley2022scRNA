min_clonotype_cell=10

setwd("/scratch/cqs/paula_hurley_projects/20210303_scRNA_human/hto_clonotype_6_fishertest")

clonotypes<-read.csv("../hto_clonotype_4_convert/result/clonotypes.csv", header=T)
rownames(clonotypes)<-clonotypes$clonotype_id
cells<-read.csv("../hto_samples_cutoff_summary/result/PH_scRNA.HTO.summary.csv", header=T)
rownames(cells)<-cells$Sample

#clonotypes<-clonotypes[c(1:10),]  

sample_columns=c(paste0(cells$Sample, "_Benign"),paste0(cells$Sample, "_Tumor"))
sample_data=clonotypes[,sample_columns]
max_column<-colnames(sample_data)[max.col(sample_data, ties.method = "first")]
max_sample<-gsub("_Benign|_Tumor","", max_column)

clonotypes$MaxSample=max_sample

benign_names=paste0(max_sample, "_Benign")
tumor_names=paste0(max_sample, "_Tumor")

clonotypes$MaxSample_Clonotype_Benign=unlist(lapply(c(1:length(max_sample)), function(x){
  return(clonotypes[x,benign_names[x]])
}))
clonotypes$MaxSample_Clonotype_Tumor=unlist(lapply(c(1:length(max_sample)), function(x){
  return(clonotypes[x,tumor_names[x]])
}))

clonotypes$MaxSample_Cell_Benign=cells[max_sample, "Benign"]
clonotypes$MaxSample_Cell_Tumor=cells[max_sample, "Tumor"]

fortest<-which(clonotypes$MaxSample_Clonotype_Benign >= min_clonotype_cell | clonotypes$MaxSample_Clonotype_Tumor >= min_clonotype_cell)
testdata<-clonotypes[fortest,]

fts<-apply(testdata, 1, function(x){
  clono_benign=as.numeric(x[["MaxSample_Clonotype_Benign"]])
  clono_tumor=as.numeric(x[["MaxSample_Clonotype_Tumor"]])
  cell_benign=as.numeric(x[["MaxSample_Cell_Benign"]])
  cell_tumor=as.numeric(x[["MaxSample_Cell_Tumor"]])
  clono_test<-matrix(c(clono_tumor, clono_benign,  cell_tumor-clono_tumor, cell_benign-clono_benign),
                     nrow = 2,
                     dimnames = list(Cell = c("Tumor", "Benign"),
                                     Clonotype = c("Yes", "No")))
  ft=fisher.test(clono_test)
  return(data.frame("pvalue"=ft$p.value,"odds_ratio"=ft$estimate))
})
ftsvalues<-do.call(rbind, fts)
ftsvalues$clonotype_id<-testdata$clonotype_id
ftsvalues$adjustp<-p.adjust(ftsvalues$pvalue, method="fdr")
ftsvalues$significant<-ftsvalues$adjustp<0.01

final<-merge(clonotypes,ftsvalues,by = "clonotype_id", all.x = TRUE)
final<-final[order(final$frequency,decreasing = T),]

write.csv(final, file="PH_scRNA.clonotype.fishertest.csv", row.names=F, na = "")
