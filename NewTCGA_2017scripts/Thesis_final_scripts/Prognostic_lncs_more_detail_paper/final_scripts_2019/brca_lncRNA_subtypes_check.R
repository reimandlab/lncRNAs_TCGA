library(survAUC)
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(EnvStats)

source("check_lnc_exp_cancers.R")

#------FEATURES-----------------------------------------------------

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

library(TCGAbiolinks)

#--------This script ------------------------------------------------

#include additional clinical variables from more detailed
#files for each cancer type 
#fit multivariate models using those variables for each candidate
#foresplots?
#correlation plots?

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

cancers = unique(allCands$canc)
get_canc_dat = function(canc){
  canc_d = subset(rna, Cancer == canc)
  return(canc_d)
}
cancer_data = llply(cancers, get_canc_dat)

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "residual_tumor", "margin_status", "PFI", "PFI.time"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)
brca_full = filtered_data[[17]]
brca_full = brca_full[,c("patient", "PFI", "PFI.time")]

#subtypes available from biolinks
subtypes_data = toupper(c("acc", "brca", "coad", "gbm", "hnsc", "kich", "kirp", 
  "kirc", "lgg", "luad", "lusc", "prad", "pancan", "read", "skcm", "stad", "thca", "ucec"))

#--------ADD CLINICAL VARIABLES----------------------------------------

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

brca = clin[[10]]
brca = merge(brca, brca_full)

brca = brca[,which(colnames(brca) %in% c("patient", "ENSG00000243926", "OS", "OS.time", "PAM50.mRNA"))]
med = median(brca$ENSG00000243926)
#median = 0
if(med==0){
brca$med = ""
brca$med[which(brca$ENSG00000243926 > 0)] = "High"
brca$med[which(brca$ENSG00000243926 == 0)] = "Low"
}
if(!(med==0)){
brca$med = ""
brca$med[which(brca$ENSG00000243926 >= med)] = "High"
brca$med[which(brca$ENSG00000243926 < med)] = "Low"
}

z = which(is.na(brca$PAM50.mRNA))
if(!(length(z)==0)){
brca = brca[-z,]}

pdf("brca_pam50_status_vs_TIPARP_AS1_km_plots.pdf", width=9)

  brca$OS = as.numeric(brca$OS)
  brca$OS.time = as.numeric(brca$OS.time)
  brca$med = factor(brca$med, levels = c("Low","High"))
 
    #lnc model
    lnc_model = coxph(Surv(OS.time, OS) ~ med, data = brca)
    hr_lnc = round(summary(lnc_model)$coefficients[2], digits=3)

    #idh model
    idh_model = coxph(Surv(OS.time, OS) ~ PAM50.mRNA, data = brca)
    hr_idh = round(summary(idh_model)$coefficients[2], digits=3)

    #lnc + idh model
    both = coxph(Surv(OS.time, OS) ~ med + PAM50.mRNA, data = brca)

    #lnc vs both
    lnc_vs_both = anova(lnc_model, both)[2,4]

    #idh vs both
    idh_vs_both = anova(idh_model, both)[2,4]
  
    brca$OS.time = brca$OS.time/365

  fit <- survfit(Surv(OS.time, OS) ~ med, data = brca)
s <- ggsurvplot(
          fit, 
          data = brca,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          facet.by ="PAM50.mRNA")
          print(s)

dev.off()

 



