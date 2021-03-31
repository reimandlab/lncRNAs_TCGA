#------------------------------------------------------------------------------

library(corrplot)

set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FUNCTIONS-----------------------------------------------------

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}


#-------------------ANALYSIS--------------------------------------------

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#pcgs = c("IDH1", "IDH2", "MGMT", "TERT", "ERBB2", "ESR1", "ATRX", "PGR",
 # "CDKN2A", "SETD2", "BAP1", "PBRM1", "PIK3CA", "ARID1A")

lncs = unique(allCands$gene)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

#library(TCGAbiolinks)

check=function(pat){
  pat_new=paste(unlist(strsplit(pat, "\\."))[1:3], collapse="-")
  test=unlist(strsplit(pat, "\\."))[4]
  print(test)
  return(pat_new)
}

#--------This script ------------------------------------------------

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

rna$residual_tumor[rna$residual_tumor == "R0"] = 1
rna$residual_tumor[rna$residual_tumor == "R1"] = 2
rna$residual_tumor[rna$residual_tumor == "R2"] = 3
rna$residual_tumor[rna$residual_tumor == "RX"] = "NA"

rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Complete Remission/Response"] = 1
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Progressive Disease"] = 4
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Stable Disease"] = 3
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Partial Remission/Response"] = 2
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "No Measureable Tumor or Tumor Markers"] = 1
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Persistent Disease"] = 4
rna$treatment_outcome_first_course[rna$treatment_outcome_first_course == "Normalization of Tumor Markers, but Residual Tumor Mass"] = 3

cancers = unique(allCands$cancer)
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
  dtt = dtt[,..z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

get_clin_lnc_cors = function(dtt){
  canc = dtt$Cancer[1]
  print(canc)
  print(dim(dtt))
  dtt=as.data.frame(dtt)
  #get lncs
  z = which(str_detect(colnames(dtt), "ENSG"))
  lncs = colnames(dtt)[z]

  #look at individual lncRNAs
  get_cor = function(lnc){
    z = which((str_detect(colnames(dtt), "ENSG") & !(colnames(dtt) %in% lnc)))
    new_dat = as.data.frame(dtt)
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

      if(!(dim(table(check_dat$lncRNA_tag))==1)){
      check_dat[,1] = as.numeric(check_dat[,1])
      check_dat$lncRNA_tag = as.numeric(check_dat$lncRNA_tag)

      tb=table(check_dat[,1], check_dat$lncRNA_tag)
      if(!(dim(tb)[1]==0)){
      chisq_pval = as.numeric(tidy(chisq.test(tb))[2])
      colnames(check_dat)[1]="variable"
      answer = c(lnc, canc, check, chisq_pval, hr)
      return(answer)
      }
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
write.csv(all_cancers_cell_types, file="/u/kisaev/cands_treatment_outcome_residual_tumor_analysis.csv", row.names=F, quote=F)
