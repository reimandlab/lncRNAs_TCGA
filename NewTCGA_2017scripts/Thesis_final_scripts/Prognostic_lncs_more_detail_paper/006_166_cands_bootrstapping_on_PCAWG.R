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

#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_May2.rds")
pcawg_data$combo = paste(pcawg_data$canc, pcawg_data$histo)

pcawg_data$canc[pcawg_data$combo == "Uterus Adenocarcinoma, endometrioid"] = "Uterine Corpus Endometrial Carcinoma" 
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, clear cell type"] = "Kidney renal clear cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Breast Infiltrating duct carcinoma"] = "Breast invasive carcinoma"
pcawg_data$canc[pcawg_data$combo == "Ovary Serous cystadenocarcinoma"] = "Ovarian serous cystadenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "Pancreas Pancreatic ductal carcinoma"] = "Pancreatic adenocarcinoma" 
pcawg_data$canc[pcawg_data$combo == "Liver Hepatocellular carcinoma"] = "Liver hepatocellular carcinoma"
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, papillary type"] = "Kidney renal papillary cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Squamous cell carcinoma"] = "Lung squamous cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Adenocarcinoma, invasive"] = "Lung adenocarcinoma"


cancers_tests = as.list(unique(cands$cancer[which(cands$cancer %in% pcawg_data$canc)]))

get_matched_data = function(cancer){
    dtt = subset(pcawg_data, canc == cancer)
    z = which(colnames(dtt) %in% c(as.character(all_cands$gene[all_cands$cancer == dtt$canc[1]]), "canc", 
    "histo", "time", "status", "sex"))
    dtt = dtt[,z]
    return(dtt)
}

filtered_data = llply(cancers_tests, get_matched_data)

add_tags = function(dtt){

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
  dat = dtt
  dat$canc = NULL
  dat$histo = NULL
  dat$sex = NULL
  dat$status[dat$status=="alive"] =0
  dat$status[dat$status=="deceased"] =1

  z = which(str_detect(colnames(dat), "ENSG"))
  check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
  }
  keep = unlist(apply(dat, 2, check_contrasts))
  dat = dat[,c(which(colnames(dat) %in% names(keep)))]
  dat$status = as.numeric(dat$status)
  dat$time = as.numeric(dat$time)
  num_genes = which(str_detect(colnames(dat), "ENSG"))
  return(dat)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")
 
#--------------------------------------------------------------------

#1, Need PCAWG data 
#2, integrate it into the internal cross-validation but instead of testing on test set -> test 
#on PCAWG for which there is data 

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
z <- which(fantom$CAT_geneName %in% rm)
rm <- fantom$CAT_geneName[z]
fantom <- fantom[-z,]


#-------------------------------------------------------------------
#------TCGA DATA---------------------------------------------------
#-------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------
#1 - get cancer data 
#-----------------------------------------------------------------------------------------------------------------------------------

z = which(cancers %in% unlist(cancers_tests))
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(cands$gene[cands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data)

#-----------------------------------------------------------------------------------------------------------------------------------
#2 - set up for cross validation 
#-----------------------------------------------------------------------------------------------------------------------------------

set_up_cv = function(dtt){
	dat = dtt
	genes_results = list()
	cinds_clin = c()
	cinds_justlncs = c()
	cinds_combined = c()
	#set.seed(101) 

	stage = as.data.table(table(dat$clinical_stage))
	rm = stage$V1[which(stage$N <10)]
	z = which(dat$clinical_stage %in% rm)
	if(!(length(z)==0)){
	dat = dat[-z,]}

	grade = as.data.table(table(dat$histological_grade))
	rm = grade$V1[which(grade$N <10)]
	z = which(dat$histological_grade %in% rm)
	if(!(length(z)==0)){
	dat = dat[-z,]}

	return(dat)
}

setup_data = llply(filtered_data, set_up_cv)










