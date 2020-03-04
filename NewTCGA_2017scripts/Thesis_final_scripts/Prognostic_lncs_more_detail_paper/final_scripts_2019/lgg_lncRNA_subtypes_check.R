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
lgg_full = filtered_data[[1]]
lgg_full = lgg_full[,c("patient", "PFI", "PFI.time")]

#subtypes available from biolinks
subtypes_data = toupper(c("acc", "lgg", "coad", "gbm", "hnsc", "kich", "kirp", 
  "kirc", "lgg", "luad", "lusc", "prad", "pancan", "read", "skcm", "stad", "thca", "ucec"))

#--------ADD CLINICAL VARIABLES----------------------------------------

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

lgg = clin[[1]]
lgg = merge(lgg, lgg_full)

lgg = lgg[,which(colnames(lgg) %in% c("patient", "ENSG00000253187", "PFI", "PFI.time", "IDH.status", "Chr.19.20.co.gain", "X1p.19q.codeletion"))]
med = median(lgg$ENSG00000253187)
#median = 0
if(med==0){
lgg$med = ""
lgg$med[which(lgg$ENSG00000253187 > 0)] = "High"
lgg$med[which(lgg$ENSG00000253187 == 0)] = "Low"
}
if(!(med==0)){
lgg$med = ""
lgg$med[which(lgg$ENSG00000253187 >= med)] = "High"
lgg$med[which(lgg$ENSG00000253187 < med)] = "Low"
}

lgg$PFI = as.numeric(lgg$PFI)
lgg$PFI.time = as.numeric(lgg$PFI.time)
lgg$med = factor(lgg$med, levels = c("Low","High"))
#write.csv(lgg, "lgg_summary_hoxa10as_molecular_features.csv", quote=F, row.names=F) 

z1=which(is.na(lgg$IDH.status))
z2=which(is.na(lgg$Chr.19.20.co.gain))
z3=which(is.na(lgg$X1p.19q.codeletion))

lgg = lgg[-unique(c(z1,z2,z3)),]

#lnc model
lnc_model = coxph(Surv(PFI.time, PFI) ~ med , data = lgg)
hr_lnc = round(summary(lnc_model)$coefficients[2], digits=3)

#idh model
idh_model = coxph(Surv(PFI.time, PFI) ~ IDH.status + Chr.19.20.co.gain + X1p.19q.codeletion, data = lgg)

#lnc + idh model
both = coxph(Surv(PFI.time, PFI) ~ med + IDH.status + Chr.19.20.co.gain + X1p.19q.codeletion, data = lgg)

#lnc vs both
lnc_vs_both = anova(lnc_model, both)[2,4]

#idh vs both
idh_vs_both = anova(idh_model, both)[2,4]
  
lgg$PFI.time = lgg$PFI.time/365
lgg$lnc_idh = paste(lgg$med, lgg$IDH.status)
lgg$Chr.19.20.co.gain = paste(lgg$med, lgg$Chr.19.20.co.gain)
lgg$X1p.19q.codeletion = paste(lgg$med, lgg$X1p.19q.codeletion)

pdf("LGG_HOXA10AS_molecular_variates.pdf", width=12, height=10)

fit <- survfit(Surv(PFI.time, PFI) ~ lnc_idh, data = lgg)

          s <- ggsurvplot(
          title = paste("LGG", "HOXA10-AS expression", "IDH status"),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

fit <- survfit(Surv(PFI.time, PFI) ~ Chr.19.20.co.gain, data = lgg)
 s <- ggsurvplot(
          title = paste("LGG", "HOXA10-AS expression", "Chr.19.20.co.gain"),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

fit <- survfit(Surv(PFI.time, PFI) ~ X1p.19q.codeletion, data = lgg)
 s <- ggsurvplot(
          title = paste("LGG", "HOXA10-AS expression", "X1p.19q.codeletion"),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

dev.off()

 



