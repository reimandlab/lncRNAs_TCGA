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

#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands


#--------This script ------------------------------------------------

#build models using all of TCGA for each lncRNA candidate
#conduct bootstrapping on PCAWG sample 
#and evluate perforamnce of TCGA candidate model 

#--------------------------------------------------------------------

#1, Need PCAWG data 
#2, integrate it into the internal cross-validation but instead of testing on test set -> test 
#on PCAWG for which there is data 











