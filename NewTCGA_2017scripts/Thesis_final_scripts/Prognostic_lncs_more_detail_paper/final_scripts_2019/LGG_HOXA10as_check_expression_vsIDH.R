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

z = which(colnames(all) %in% c("ENSG00000253187", "patient", "Cancer", "PFI", "OS", "PFI.time", "OS.time"))
lnc = all[,z]
lnc=subset(lnc, Cancer %in% c("Brain Lower Grade Glioma"))
lnc$exp_log1p = log1p(lnc$ENSG00000253187)

#boxplot

pdf("/u/kisaev/lgg_hoxa10as_expresion.pdf", width=10, height=10)

#--------ADD CLINICAL VARIABLES----------------------------------------

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

lgg = clin[[1]]
lgg = merge(lgg, lnc, by="patient")

lgg = lgg[,which(colnames(lgg) %in% c("patient", "ENSG00000253187.y", "OS.y", "OS.time.y", "PFI", "PFI.time", "IDH.status", "Chr.19.20.co.gain", "X1p.19q.codeletion"))]
lgg = subset(lgg, !(is.na(IDH.status)))
colnames(lgg)[c(5:6, 9)] = c("OS", "OS.time", "ENSG00000253187")
lgg$exp_log1p = log1p(lgg$ENSG00000253187)

  ggboxplot(lgg, x="IDH.status", y="ENSG00000253187", fill="IDH.status") + stat_compare_means() + stat_n_text()
  ggboxplot(lgg, x="IDH.status", y="exp_log1p", fill="IDH.status") + stat_compare_means() + stat_n_text()
   
  lgg = subset(lgg, IDH.status == "WT")

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

  lgg$OS = as.numeric(lgg$OS)
  lgg$OS.time = as.numeric(lgg$OS.time)
  lgg$med = factor(lgg$med, levels = c("Low","High"))

  lgg$OS.time = lgg$OS.time/365
  #KM plot
  fit <- survfit(Surv(OS.time, OS) ~ med, data = lgg)

            s <- ggsurvplot(
            title = paste("LGG IDH WT", "HOXA10-AS expression"),
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
            palette = "npg" 
            #ggtheme = theme_bw(), # customize plot and risk table with a theme.
            #risk.table.y.text.col = T, # colour risk table text annotations.
            #risk.table.y.text = FALSE # show bars instead of names in text annotations
                              # in legend of risk table
            )
            print(s)

  ggboxplot(lgg, x="med", y="ENSG00000253187", fill="med", title="IDH WT LGGs FPKM-UQ expression")
  dev.off()

