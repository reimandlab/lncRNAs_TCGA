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
library(stringr)

###---------------------------------------------------------------
###Results - using just 202 patients 
###---------------------------------------------------------------

#1. 100 CV just clinical variables 
##cindeces
genes = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch1.rds")
genes = unlist(genes)
genes = as.data.table(table(genes))
genes = genes[order(N)]

genes = subset(genes, N >=20)







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













































