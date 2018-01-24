###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:5790)]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

library(glmnet)
library(survcomp)
library(caret)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]]
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

#clinical variables 
add_clin = fread("OV_clinical_core.txt", data.table=F)
clin = canc_data[,2359:ncol(canc_data)]

#clinical variables 
add_clin = fread("OV_clinical_core.txt", data.table=F)
add_clin = add_clin[,-5]
colnames(add_clin)[1] = "patient"
add_clin = merge(add_clin, clin, by="patient")
rownames(add_clin) = add_clin$patient
z = which(add_clin$stage == "[Not Available]")
add_clin = add_clin[-z,]
z = which(add_clin$grade %in% c("G4", "GB", "GX"))
add_clin = add_clin[-z,]

add_clin$stage[add_clin$stage == "IIA"] = 2
add_clin$stage[add_clin$stage == "IIB"] = 2
add_clin$stage[add_clin$stage == "IIC"] = 2
add_clin$stage[add_clin$stage == "IIIA"] = 3
add_clin$stage[add_clin$stage == "IIIB"] = 3
add_clin$stage[add_clin$stage == "IIIC"] = 3
add_clin$stage[add_clin$stage == "IV"] = 4

add_clin$grade[add_clin$grade == "G1"] = 1
add_clin$grade[add_clin$grade == "G2"] = 2
add_clin$grade[add_clin$grade == "G3"] = 3

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% add_clin$patient)
add_clin = add_clin[,-c(6,7)]

gene_data = t(canc_data[,1:(ncol(rna)-5)])
	
#1. remove any genes that have 0 counts within cancer
sums = apply(canc_data[,1:(ncol(rna)-5)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(canc_data) %in% zeroes)
if(!(length(z)==0)){
	canc_data = canc_data[,-z]
}

###---------------------------------------------------------------
###Results - using just 202 patients 
###---------------------------------------------------------------

#1. 100 CV just clinical variables 
##cindeces
m2 = readRDS("MODEL2_updated_code_CLINICAL_predictors_cindes_CV_100timesJan17.RDS")
m2 = as.data.frame(m2)
colnames(m2) = "cindex"
m2$Predictor = "Clinical"
m2$Model = 1

##cindeces
m8A = readRDS("MODEL8A_202_OV_pats_3bestGenes_CV_100timesJan19.RDS")
m8A = as.data.frame(m8A)
colnames(m8A) = "cindex"
m8A$Predictor = "AC018647.1"
m8A$Model = 2

##cindeces
m9A = readRDS("MODEL9A_202_OV_pats_3bestGenes_CV_100timesJan19.RDS")
m9A = as.data.frame(m9A)
colnames(m9A) = "cindex"
m9A$Predictor = "LINC00263"
m9A$Model = 3

##cindeces
m10A = readRDS("MODEL10A_202_OV_pats_3bestGenes_CV_100timesJan19.RDS")
m10A = as.data.frame(m10A)
colnames(m10A) = "cindex"
m10A$Predictor = "RP11-443B7.3"
m10A$Model = 4

##cindeces
m6C = readRDS("MODEL6C_202_OV_pats_3bestGenes_CV_100timesJan19.RDS")
m6C = as.data.frame(m6C)
colnames(m6C) = "cindex"
m6C$Predictor = "Combined_lncRNAs"
m6C$Model = 5

##cindeces
m6 = readRDS("MODEL6_202_OV_pats_3bestGenesANDclinicalVariables_CV_100timesJan19.RDS")
m6 = as.data.frame(m6)
colnames(m6) = "cindex"
m6$Predictor = "Combined_lncRNAs_and_Clinical"
m6$Model = 6

cinds_all = rbind(m2, m8A, m9A, m10A, m6C, m6)
cinds_all$Predictor = factor(cinds_all$Predictor, levels=c("Clinical",
	"AC018647.1", "LINC00263","RP11-443B7.3","Combined_lncRNAs","Combined_lncRNAs_and_Clinical"))

###---------------------------------------------------------------
###Plot C-Indeces for the different models 
###---------------------------------------------------------------

pdf("202_ovarian_cancer_patients_comparing_model_performanceJan19.pdf", width=12, height=12)
g = ggboxplot(cinds_all, x="Model", y= "cindex", fill="Predictor", palette=mypal, 
	notch=T, main="202 Ovarian Cancer Patients, 100 Cross Validations", 
	order=c("1", "2", "3", "4", "5", "6"), 
	legend="right", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15)+
		stat_compare_means(method = "wilcox.test",
                     ref.group = "1", label = "p.format")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")
ggpar(g, xlab ="Model", ylab="C-index")
print(g)
dev.off()

###---------------------------------------------------------------
###Results - using all 305 patients 
###---------------------------------------------------------------

m8B = readRDS("MODEL8B_305_OV_pats_bestGenes_CV_100timesJan19.RDS")
m8B = as.data.frame(m8B)
colnames(m8B) = "cindex"
m8B$Predictor = "AC018647.1"
m8B$Model = 1

m9B = readRDS("MODEL9B_305_OV_pats_bestGenes_CV_100timesJan19.RDS")
m9B = as.data.frame(m9B)
colnames(m9B) = "cindex"
m9B$Predictor = "LINC00263"
m9B$Model = 2

m10B = readRDS("MODEL10B_305_OV_pats_bestGenes_CV_100timesJan19.RDS")
m10B = as.data.frame(m10B)
colnames(m10B) = "cindex"
m10B$Predictor = "RP11-443B7.3"
m10B$Model = 3

m6B = readRDS("MODEL6B_305_OV_pats_3bestGenes_CV_100timesJan19.RDS")
m6B = as.data.frame(m6B)
colnames(m6B) = "cindex"
m6B$Predictor = "Combined_lncRNAs"
m6B$Model = 4

cinds_all = rbind(m8B, m9B, m10B, m6B)
cinds_all$Predictor = factor(cinds_all$Predictor, levels=c("AC018647.1", 
	"LINC00263","RP11-443B7.3","Combined_lncRNAs"))

pdf("305_ovarian_cancer_patients_comparing_model_performanceJan19.pdf", width=12, height=12)
g = ggboxplot(cinds_all, x="Model", y= "cindex", fill="Predictor", palette=mypal, 
	notch=T, main="305 Ovarian Cancer Patients, 100 Cross Validations", 
	order=c("1", "2", "3", "4"), 
	legend="right", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15)+
		stat_compare_means(method = "wilcox.test",
                     ref.group = "4", label = "p.format")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")
ggpar(g, xlab ="Model", ylab="C-index")
print(g)
dev.off()













































