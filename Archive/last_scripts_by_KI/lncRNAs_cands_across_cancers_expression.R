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
  ids = which(!(str_detect(colnames(dtt), "ENSG")))
  lncs = which((str_detect(colnames(dtt), "ENSG")))

  dtt.m1 = melt(dtt, id.vars = ids,
                measure.vars = lncs)
  dtt.m1$lnc = sapply(dtt.m1$variable, get_name)
  dtt.m1 = as.data.table(dtt.m1)
  return(dtt.m1)
}

filtered_data = as.data.table(ldply(llply(cancer_data, get_canc_data_for_plot)))
filtered_data$exp = log1p(filtered_data$value)
filtered_data = merge(filtered_data, canc_conv)

#plot lncRNAs from lowest median expression to highest median expression fill by cancer type
#x=lnc
#y=log1p exp
#fill=cancer

lncs = unique(filtered_data$lnc)

get_med = function(lncrna){
  dat = as.data.table(filter(filtered_data, lnc==lncrna))
  med = median(dat$exp)
  return(c(lncrna, med))
}

meds = as.data.table(ldply(llply(lncs, get_med)))
meds = meds[order(-V2)]
filtered_data$lnc = factor(filtered_data$lnc, levels=meds$V1)

#plot

pdf("/u/kisaev/figureSuppX_lncRNA_cands_expression.pdf", width=10, height=6)
g=ggerrorplot(filtered_data, x="lnc", y="exp", fill="type", x.text.angle = 90,  add = c("median", "median_mad"), error.plot = "errorbar")
ggpar(g, font.xtickslab=c(5,"plain", "black"))+theme_bw()
dev.off()



