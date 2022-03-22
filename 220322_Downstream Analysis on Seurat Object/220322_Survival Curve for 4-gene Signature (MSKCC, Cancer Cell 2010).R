#Survival Curve for 4-gene Signature (MSKCC, Cancer Cell 2010)
#Data from: https://www.cbioportal.org/study/summary?id=prad_mskcc
setwd("<folder with all relevant unzipped MSKCC files>")
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
read.delim("<filepath to data_clinical_combined_patient_sample_MSKCC2010.txt>", sep="\t")->data
newdata <- data[-c(1, 2, 3, 4), ] 
newdata<-separate(data = newdata, col = Disease.Free.Status, into = c("Status", "DiseaseFree"), sep = ":")
newdata$Status<-as.numeric(newdata$Status)
Zscores<-read.delim("<filepath to data_mrna_agilent_microarray_zscores_ref_diploid_samples_WHY.txt>")
meta<-newdata[newdata$X.Patient.Identifier %in% colnames(Zscores)[2:ncol(Zscores)],]
Zscores2<-Zscores[,meta$X.Patient.Identifier]
Zscores2$gene<-Zscores$Entrez_Gene_Id
#Using Entrez Gene IDs instead of Hugo Gene Symbols
#ASPN = 54829
#CTHRC1 = 115908
#ENG = 2022
#FAP = 2191
Geneset <- Zscores2[Zscores2$gene %in% c("54829","115908","2022","2191"), ]
Geneset$gene<-NULL
Geneset_Mean<-as.data.frame(colMeans(Geneset))
colnames(Geneset_Mean)<-"Zscore Mean"
meta$`Zscore Mean`<-Geneset_Mean$`Zscore Mean`
meta$`Zscore Mean`<- as.numeric(meta$`Zscore Mean`)
summary(meta$`Zscore Mean`)
quantile(meta$`Zscore Mean`)

#IF only median (2-lined curve) use next line
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`>0.02649987, "High", "Low"))
#IF only median (2-lined curve) use line above

#IF using quartiles (4-lined curve) use next 3 lines
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`> 0.7423112, "4", "3"))
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`< 0.0264998, "2", Cat))
meta<-mutate(meta, Cat=ifelse(meta$`Zscore Mean`< -0.3653848, "1", Cat))
#IF using quartiles (4-lined curve) use 3 lines above

meta$Disease.Free..Months.<- as.numeric(meta$Disease.Free..Months.)
fit<-survfit(Surv(Disease.Free..Months., Status) ~ Cat, data=meta)
print(fit)
ggsurvplot(fit,
           data=meta,
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           ggtheme = theme_bw(),
           font.y=15, font.x=15, font.tickslab="black",
           palette = c("#E7B800", "#2E9FDF", "#86AA00", "#FB8072"))