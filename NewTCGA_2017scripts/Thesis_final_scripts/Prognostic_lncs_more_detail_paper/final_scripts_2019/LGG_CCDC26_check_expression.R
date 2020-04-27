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

z = which(colnames(all) %in% c("ENSG00000229140", "patient", "Cancer", "PFI", "OS", "PFI.time", "OS.time"))
lnc = all[,z]
lnc=subset(lnc, Cancer %in% c("Brain Lower Grade Glioma", "Glioblastoma multiforme"))
lnc$CCDC26_log1p = log1p(lnc$ENSG00000229140)

maf = fread("/u/kisaev/LGG_maf_idh_tp53.txt")

#For the survival analysis, can you please sub-stratify the 
#IDH mutant cases in p53wt and p53 mutant cases (which gives you 
#LGG type 1 versus type 2) and then look at low and high CCDC6 
#expression – basically you will have 4 groups:

#IDHmt; p53 wt CCDC6 low versus IDHmt; p53 wt CCDC6 high
#and
#IDHmt; p53 mt CCDC6 low versus IDHmt; p53 mt) CCDC6 high

#boxplot

pdf("lgg_CCDC26_expresion.pdf", width=10, height=10)

ggboxplot(lnc, x="Cancer", y="ENSG00000229140", fill="Cancer") + stat_compare_means() + stat_n_text()
ggboxplot(lnc, x="Cancer", y="CCDC26_log1p", fill="Cancer") + stat_compare_means() + stat_n_text()

#--------ADD CLINICAL VARIABLES----------------------------------------

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

lgg = clin[[1]]
lgg = merge(lgg, lnc, by="patient")

lgg = lgg[,which(colnames(lgg) %in% c("patient", "ENSG00000229140", "OS.y", "OS.time.y", "PFI", "PFI.time", "IDH.status", "Chr.19.20.co.gain", "X1p.19q.codeletion"))]
lgg = subset(lgg, !(is.na(IDH.status)))
colnames(lgg)[5:6] = c("OS", "OS.time")

  ggboxplot(lgg, x="IDH.status", y="ENSG00000229140", fill="IDH.status") + stat_compare_means() + stat_n_text()
  ggboxplot(lgg, x="IDH.status", y="CCDC26_log1p", fill="IDH.status") + stat_compare_means() + stat_n_text()
   
  lgg = subset(lgg, IDH.status == "Mutant")

  med = median(lgg$ENSG00000229140)

  #median = 0
  if(med==0){
  lgg$med = ""
  lgg$med[which(lgg$ENSG00000229140 > 0)] = "High"
  lgg$med[which(lgg$ENSG00000229140 == 0)] = "Low"
  }
  if(!(med==0)){
  lgg$med = ""
  lgg$med[which(lgg$ENSG00000229140 >= med)] = "High"
  lgg$med[which(lgg$ENSG00000229140 < med)] = "Low"
  }

  lgg$PFI = as.numeric(lgg$PFI)
  lgg$PFI.time = as.numeric(lgg$PFI.time)
  lgg$med = factor(lgg$med, levels = c("Low","High"))

  lgg$PFI.time = lgg$PFI.time/365
  #KM plot
  fit <- survfit(Surv(PFI.time, PFI) ~ med, data = lgg)

            s <- ggsurvplot(
            title = paste("LGG IDH mutants", "CCDC26 expression"),
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

  ggboxplot(lgg, x="med", y="ENSG00000229140", fill="med", title="IDH mut LGGs FPKM-UQ expression")
  dev.off()

#merge with maf file
lgg=as.data.table(lgg)
lgg$OS.time = lgg$OS.time/365
  
colnames(maf)[16] = "patient"
maf$idh_mut = ""
maf$idh_mut[maf$Hugo_Symbol %in% c("IDH1", "IDH2")] = "yes" 

z1 = which(lgg$patient %in% maf$patient[maf$idh_mut=="yes"])
z2 = which(lgg$patient %in% maf$patient[maf$Hugo_Symbol == "TP53"])

lgg$idh = ""
lgg$tp53 = ""
lgg$idh[z1] = "IDH"
lgg$tp53[z2] = "TP53mut"
lgg$tp53[-z2] = "TP53WT"

lgg = as.data.table(filter(lgg, idh=="IDH"))
lgg_53 = as.data.table(filter(lgg, tp53=="TP53"))
lgg_no53 = as.data.table(filter(lgg, !(tp53=="TP53")))

pdf("/u/kisaev/LGG_CCDC26_expresion_TP53_boxplot.pdf")
ggboxplot(lgg, x="tp53", y="ENSG00000229140", add="jitter", fill="tp53", title="IDH mut LGGs FPKM-UQ expression vs TP53")+stat_n_text()+
stat_compare_means() + theme_bw()
dev.off()



med = median(lgg_53$ENSG00000229140)

  #median = 0
  if(med==0){
  lgg_53$med = ""
  lgg_53$med[which(lgg_53$ENSG00000229140 > 0)] = "High"
  lgg_53$med[which(lgg_53$ENSG00000229140 == 0)] = "Low"
  }
  if(!(med==0)){
  lgg_53$med = ""
  lgg_53$med[which(lgg_53$ENSG00000229140 >= med)] = "High"
  lgg_53$med[which(lgg_53$ENSG00000229140 < med)] = "Low"
  }


med = median(lgg_no53$ENSG00000229140)

  #median = 0
  if(med==0){
  lgg_no53$med = ""
  lgg_no53$med[which(lgg_no53$ENSG00000229140 > 0)] = "High"
  lgg_no53$med[which(lgg_no53$ENSG00000229140 == 0)] = "Low"
  }
  if(!(med==0)){
  lgg_no53$med = ""
  lgg_no53$med[which(lgg_no53$ENSG00000229140 >= med)] = "High"
  lgg_no53$med[which(lgg_no53$ENSG00000229140 < med)] = "Low"
  }



pdf("/u/kisaev/LGG_CCDC26_expresion_TP53.pdf")

#IDHmt; p53 wt CCDC6 low versus IDHmt; p53 wt CCDC6 high
fit <- survfit(Surv(PFI.time, PFI) ~ med, data = lgg_no53)
  s <- ggsurvplot(
          title = paste("LGG IDH mutants TP53 WT", "CCDC26 expression"),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          data = lgg_no53,      # data used to fit survival curves. 
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
#and
#IDHmt; p53 mt CCDC6 low versus IDHmt; p53 mt) CCDC6 high
fit <- survfit(Surv(PFI.time, PFI) ~ med, data = lgg_53)
  s <- ggsurvplot(
          title = paste("LGG IDH mutants TP53 MT", "CCDC26 expression"),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          data = lgg_53,      # data used to fit survival curves. 
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

dev.off()







