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
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

dim(rna)
rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% "type")
rna = rna[,c(z1, z2)]

#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, sum)
z = which(sums ==0)
if(!(length(z)==0)){
rm = names(sums)[z]}

#2. remove those with MAD < 0? 
#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, mad)
z = which(sums <= 0)
rna = rna[,-z]

#2. cluster cancer types, ie -> each dot should be cancer type 
library("FactoMineR")
z1 = which(str_detect(colnames(rna), "ENSG"))
logged_rna = rna
logged_rna[,z1] = log1p(logged_rna[,z1])

#logged - DONE! 
#pdf("logged_nonMAD0lncRNAs_cands_PCA_plots_June25.pdf", width=12)
#autoplot(prcomp(logged_rna[,z1]), data = logged_rna, colour = 'type')
#dev.off()


#NOT RUN ---------------------------------------------------------------
#not logged & scaled 
#pdf("scaled_nonMAD0lncRNAs_cands_PCA_plots_June25.pdf", width=12)
#z1 = which(str_detect(colnames(rna), "ENSG"))
#autoplot(prcomp(rna[,z1], scale. = TRUE), data = rna, colour = 'type')
#dev.off()
#PCA--------------------------------------------------------------------
#NOT RUN END------------------------------------------------------------


#-------------------------------------------------------------------
#-----------------PCA using mean/lncRNA/cancer to get --------------
#-----------------		32 points in the end          --------------
#-------------------------------------------------------------------





















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

