#Survival curve for 4-gene signature (TCGA, PanCancer Atlas)
#Data from: https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
setwd("<folder with all relevant unzipped TCGA files>")
read.delim("<filepath to prad_tcga_clinical_data.tsv>", sep="\t")->data
newdata <- data[-c(1, 2, 3, 4), ]
newdata$altid<-gsub("-", ".", newdata$X.Patient.Identifier)
newdata<-separate(data = newdata, col = Progression.Free.Status, into = c("Status", "ProgressionFree"), sep = ":")
newdata$Status<-as.numeric(newdata$Status)
Zscores<-read.delim("<filpath to data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt>")
Zscores$Entrez_Gene_Id<-NULL
colnames(Zscores) <-gsub(".01","", colnames(Zscores))
meta<-newdata[newdata$altid %in% colnames(Zscores)[2:ncol(Zscores)],]
Zscores2<-Zscores[,meta$altid]
Zscores2$gene<-Zscores$Hugo_Symbol
Geneset <- Zscores2[Zscores2$gene %in% c("ASPN","CTHRC1","ENG","FAP"), ]
Geneset$gene<-NULL
Geneset_Mean<-as.data.frame(colMeans(Geneset))
colnames(Geneset_Mean)<-"Zscore Mean"
meta$`Zscore Mean`<-Geneset_Mean$`Zscore Mean`
meta$`Zscore Mean`<- as.numeric(meta$`Zscore Mean`)
summary(meta$`Zscore Mean`)
quantile(meta$`Zscore Mean`)
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`> 0.6468875, "4", "3"))
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`< -0.0163875, "2", Cat))
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`< -0.6621563, "1", Cat))
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`>-0.0163875, "High", "Low"))
meta$Progress.Free.Survival..Months.<-as.numeric(meta$Progress.Free.Survival..Months.)
fit<-survfit(Surv(Progress.Free.Survival..Months., Status) ~ Cat, data=meta)
print(fit)
ggsurvplot(fit,
           data=meta,
           pval = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           ggtheme = theme_bw(),
           font.y=15, font.x=15, font.tickslab="black",
           palette = c("#E7B800", "#2E9FDF", "#86AA00", "#FB8072"))