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
###Results 
###---------------------------------------------------------------

#1. 100 CV just clinical variables 
##cindeces
c1 = readRDS("updated_code_CLINICAL_predictors_cindes_CV_100timesJan17.RDS")
c1 = as.data.frame(c1)
colnames(c1) = "cindex"
c1$test = "clinical"

#2. 100 CV just lncRNAs (using all 305 patients)
##genes 
genes2 = readRDS("ALL_OV_pats_updated_code_binary_predictors_list_of_sig_genes_CVJan17.rds")
genes2 = unlist(genes2)
genes2 = as.data.table(table(genes2))
genes2 = genes2[order(N)]
##cindeces 
c2 = readRDS("ALL_OV_pats_updated_code_binary_predictors_cindes_CV_100timesJan17.RDS")
c2 = as.data.frame(c2)
colnames(c2) = "cindex"
c2$test = "lncRNAs_all305patients"

#3. 100 CV just lncRNAs (using just 205 patients)
##genes 
genes3 = readRDS("only_pats_wclinical_updated_code_binary_predictors_list_of_sig_genes_CVJan17.rds")
genes3 = unlist(genes3)
genes3 = as.data.table(table(genes3))
genes3 = genes3[order(N)]
##cindeces 
c3 = readRDS("only_pats_wclinical_updated_code_binary_predictors_cindes_CV_100timesJan17.RDS")
c3 = as.data.frame(c3)
colnames(c3) = "cindex"
c3$test = "lncRNAs_just205wclinicalPatients"

#4. 100 CV lncRNAs + clinical 
##genes
genes4 = readRDS("lncRNAS_AND_wclinical_updated_code_binary_predictors_list_of_sig_genes_CVJan16.rds")
##cindices 
c4 = readRDS("lncRNAS_AND_wclinical_updated_code_binary_predictors_cindes_CV_100timesJan16.RDS")
c4 = as.data.frame(c4)
colnames(c4) = "cindex"
c4$test = "combined"

cinds_all = rbind(c1, c2, c3, c4)
















