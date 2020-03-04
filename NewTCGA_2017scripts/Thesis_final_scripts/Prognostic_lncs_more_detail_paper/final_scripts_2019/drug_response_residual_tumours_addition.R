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

#final script used to generate clinical analyiss data now figure 4

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

library(TCGAbiolinks)

check=function(pat){
  pat_new=paste(unlist(strsplit(pat, "\\."))[1:3], collapse="-")
  test=unlist(strsplit(pat, "\\."))[4]
  print(test)
  return(pat_new)
}

#--------This script ------------------------------------------------

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

#clean up extra new columns 
rna$new_tumor_event_type[rna$new_tumor_event_type == "#N/A" ] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "[Not Available]"] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "[Unknown]"] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "[Not Applicable]"] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "[Discrepancy]"] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "#N/A"] = "NA"
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "[Not Evaluated]"] = "NA"
rna$residual_tumor[rna$residual_tumor == "#N/A"] = "NA"
rna$residual_tumor[rna$residual_tumor == "[Not Available]"] = "NA"
rna$residual_tumor[rna$residual_tumor == "[Unknown]"] = "NA"
rna$margin_status[rna$margin_status == "#N/A"] = "NA"
rna$margin_status[rna$margin_status == "[Unknown]"] = "NA"
rna$margin_status[rna$margin_status == "[Not Available]"] = "NA"

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

get_clin_lnc_cors = function(dtt){
  canc = dtt$Cancer[1]
  print(canc)
  print(dim(dtt))
  #get lncs
  z = which(str_detect(colnames(dtt), "ENSG")) 
  lncs = colnames(dtt)[z]
  
  #look at individual lncRNAs 
  get_cor = function(lnc){
    z = which((str_detect(colnames(dtt), "ENSG") & !(colnames(dtt) %in% lnc)))
    new_dat = dtt
    if(length(z) > 0){
    new_dat = dtt[,-z]}
    #add 0/1 labels 
    new_dat$lncRNA_tag = ""
    med = median(new_dat[,which(colnames(new_dat) %in% lnc)])
    k = which(colnames(new_dat) %in% lnc)
    if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(new_dat[,k] > 0)
        l2 = which(new_dat[,k] ==0)
        new_dat$lncRNA_tag[l1] = 1
        new_dat$lncRNA_tag[l2] = 0
        }

        if(!(med ==0)){
        l1 = which(new_dat[,k] >= med)
        l2 = which(new_dat[,k] < med)
        new_dat$lncRNA_tag[l1] = 1
         new_dat$lncRNA_tag[l2] = 0
        }
    #get risk type 
    z = as.numeric(which((allCands$cancer %in% canc) & (allCands$gene %in% lnc) & (allCands$data == "TCGA")))
    hr = as.numeric(allCands$HR[z])
    new_dat$risk = ""
    if(hr >1){new_dat$risk = "HighExp"}
    if(hr <1){new_dat$risk = "LowExp"}
    
    hr = as.numeric(allCands$HR[allCands$gene == lnc])
    if(hr > 1){
      hr="high_exp_bad"
    }
    if(hr < 1){
      hr="low_exp_bad"
    }

    #look association between lncRNA and new variable 
    variables_check = c("treatment_outcome_first_course", 
      "residual_tumor")

    get_assoc = function(check){
      z1 = which(colnames(new_dat) == check)
      z2 = which(colnames(new_dat) %in% c("lncRNA_tag", "OS", "OS.time", "PFI", "PFI.time"))
      check_dat = new_dat[,c(z1,z2)]
      z = which(check_dat[,1] == "NA")
      if(!(length(z)==0)){
        check_dat=check_dat[-z,]
      }
      tb=table(check_dat[,1], check_dat$lncRNA_tag)
      if(!(dim(tb)[1]==0)){
      chisq_pval = as.numeric(tidy(chisq.test(tb))[2])
      answer = c(lnc, canc, check, chisq_pval, hr)
      return(answer)
      }
    }
    variables_checked = as.data.table(ldply(llply(variables_check, get_assoc)))
    return(variables_checked)
    }

    canc_results = as.data.table(ldply(llply(lncs, get_cor)))
    return(canc_results)
  }

all_cancers_cell_types = as.data.table(ldply(llply(filtered_data, get_clin_lnc_cors)))
colnames(all_cancers_cell_types) = c("lnc", "canc", "check", "chisq_pval", "hr")
all_cancers_cell_types$chisq_pval = as.numeric(all_cancers_cell_types$chisq_pval)
all_cancers_cell_types$chisq_fdr = p.adjust(all_cancers_cell_types$chisq_pval, method="fdr")
write.csv(all_cancers_cell_types, file="cands_treatment_outcome_residual_tumor_analysis.csv", row.names=F, quote=F)




