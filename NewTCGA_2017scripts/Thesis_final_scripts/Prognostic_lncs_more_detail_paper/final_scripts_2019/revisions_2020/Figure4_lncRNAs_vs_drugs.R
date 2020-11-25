source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#library(TCGAbiolinks)

#--------This script ------------------------------------------------

#include additional clinical variables from more detailed
#files for each cancer type 
#fit multivariate models using those variables for each candidate
#foresplots?
#correlation plots?

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group


drug_data = readRDS("all_cancers_drug_info.rds")
drug_data = as.data.table(ldply(drug_data))
colnames(drug_data)[1] = "patient"

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
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,..z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)


#--------ADD CLINICAL VARIABLES----------------------------------------

add_clin_vars = function(dtt){
  canc = dtt$Cancer[1]
  canc = rna$type[rna$Cancer == canc][1]
  dtt$cancer = canc

  dtt_drug = merge(dtt, drug_data, by=c("patient", "cancer"))
  #Check if TCGA has 
  z = which(drug_data$patient %in% dtt$patient)
  if(!(length(z)==0)){
  clin_subtypes <- TCGAquery_subtype(tumor = canc)

  clin = clin_subtypes
  z = which(clin$patient %in% dtt$patient)
  clin = clin[z,]

    #check which columns have enoguh contrasts
    #remove columns where # of NAs is greater than 50% of patietnt cohort
    
    check_nas = function(col){
      check = length(which((col == "[Not Applicable]") | (col == "[Not Available]") | (col == "Unknown")))
        if(check < (dim(clin)[1] *0.5)){
          return("keep")
        }
      }
    
    keep_cols1 = unlist(apply(clin, 2, check_nas))
    clin = clin[,which(colnames(clin) %in% names(keep_cols1))]

    check_contrasts = function(col){
      check = dim(table(col))
        if(check >1){
          return("keep")
        }
      }
    
    keep_cols2 = unlist(apply(clin, 2, check_contrasts))
    clin = clin[,which(colnames(clin) %in% names(keep_cols2))]

    cols = colnames(clin)[which(colnames(clin) %in% colnames(dtt))]

    dtt = merge(dtt, clin, by=cols)
    return(dtt)

    } #end add_clin_vars 

  #if not in molecular profiles subset of biolinks
  #just look at whatever clinical variables are available
  #if(length(z)==0){
  #  clinical <- GDCquery_clinic(project = paste("TCGA-", canc, sep=""), type = "clinical")
  #}

}

clin_data_lncs = llply(filtered_data, add_clin_vars)

#remove Nulls
clin_data_lncs = Filter(Negate(is.null), clin_data_lncs)
#saved file --- below
saveRDS(clin_data_lncs, file="clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")




