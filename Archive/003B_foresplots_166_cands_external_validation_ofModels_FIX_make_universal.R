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
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
#cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
#cands$Cancer = NULL
all_cands = cands

#--------This script ------------------------------------------------

#for each cancer type 
#look at respective lncrna cands 
#evluate their relative survival association using 
#foresplots 

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% all_cands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(all_cands$gene[all_cands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

add_tags = function(dtt){
  print(dtt$Cancer[1])
  rownames(dtt) = dtt$patient
  dtt$patient = NULL

  #log1p 
  z = which(str_detect(colnames(dtt), "ENSG"))
  if(length(z)>1){
  medians = apply(dtt[,z], 2, median)}
  if(length(z)==1){
    medians = median(dtt[,z])
  }
  #add high low tag
  for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dtt[,k] > 0)
    l2 = which(dtt[,k] ==0)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(dtt[,k] >= med)
    l2 = which(dtt[,k] < med)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }
  }  
  return(dtt)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")

get_survival_models = function(dtt){
  results_cox1 <- as.data.frame(matrix(ncol=7)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer")

  dat = dtt
  dat$Cancer = NULL
  dat$new_tumor_event_type = NULL
  dat$treatment_outcome_first_course = NULL

  num_genes = which(str_detect(colnames(dtt), "ENSG"))
  check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
  }
  keep = unlist(apply(dat, 2, check_contrasts))
  
  dat = dat[,c(which(colnames(dat) %in% names(keep)))]

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)
  dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)
  z = which(colnames(dat) == "age_at_initial_pathologic_diagnosis")
  colnames(dat)[z] = "Age"

  #z = which((!(str_detect(colnames(dat), "ENSG"))) & (!(colnames(dat) %in% c("Age", "OS", "OS.time"))))
  #if(!(length(z)==0)){
   # for(i in 1:length(z)){
    #  dat[,z[i]] = as.factor(dat[,z[i]])
    #}
  #}

  #how to summarize this 
  model <- coxph( Surv(OS.time, OS) ~ ., data = dat)
  print(colnames(dat)[1])
  print(ggforest(model, main = paste(rna$type[which(rna$Cancer %in% dtt$Cancer[1])], "Forest Plot")[1]))

return(results_cox1)

}

pdf("166_cands_foreplots.pdf", width=10, height=11)
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()




