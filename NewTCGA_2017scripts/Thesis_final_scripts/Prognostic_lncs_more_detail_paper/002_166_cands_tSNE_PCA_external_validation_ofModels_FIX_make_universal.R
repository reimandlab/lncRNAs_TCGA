library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library("FactoMineR")


#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands

#-------------------------------------------------------------------
#-----------------PCA using just 166 lncRNA candidates--------------
#-------------------------------------------------------------------

#subset to cancers,  n = 23
z = which(rna$Cancer %in% cands$cancer)

#subset to lncRNAs, keep cancer type 
rna = rna[z,]

z1 = which(colnames(rna) %in% cands$gene)
z2 = which(colnames(rna) %in% c("Cancer", "patient", "type"))

rna = rna[,c(z1, z2)]
rownames(rna) = rna$patient
rna$patient = NULL
rna$Cancer = NULL

#log transform 
#rna[,1:(ncol(rna)-1)] = log1p(rna[,1:(ncol(rna)-1)]) 

#calculate PCA
rna.pca <- PCA(rna[,-ncol(rna)], graph = FALSE)

pdf("166_cands_PCA_plots_June15.pdf", width=12)


fviz_pca_ind(rna.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = rna$type, # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
             ) + scale_shape_manual(values=seq(0,23)) + theme_bw()


dev.off()

