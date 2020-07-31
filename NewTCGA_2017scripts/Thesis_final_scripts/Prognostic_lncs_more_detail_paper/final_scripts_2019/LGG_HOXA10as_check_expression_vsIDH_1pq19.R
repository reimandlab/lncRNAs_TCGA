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

#library(TCGAbiolinks)

maf = fread("/u/kisaev/LGG_maf_idh_tp53.txt")
colnames(maf)[16] = "patient"
maf$idh_mut = ""
maf$idh_mut[maf$Hugo_Symbol %in% c("IDH1", "IDH2")] = "yes" 

#--------This script ------------------------------------------------

z = which(colnames(all) %in% c("ENSG00000253187", "patient", "Cancer", "PFI", "OS", "PFI.time", "OS.time"))
lnc = all[,z]
lnc=subset(lnc, Cancer %in% c("Brain Lower Grade Glioma"))
lnc$exp_log1p = log1p(lnc$ENSG00000253187)

#boxplot

#--------ADD CLINICAL VARIABLES----------------------------------------

#take IDH mutants and then stratify by 1pq19 codel
#then stratify by HOXA10AS
#IDH mutants with no 1pq19 stratify by HOXA10as high low 
#stratify by ATRX and TERT

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

lgg = clin[[1]]
lgg = merge(lgg, lnc, by="patient")

lgg = lgg[,which(colnames(lgg) %in% c("patient", "ENSG00000253187.y", "OS.y", "OS.time.y", "PFI", "PFI.time", 
  "IDH.status", "Chr.19.20.co.gain", "X1p.19q.codeletion", "TERT.promoter.status", 
  "ATRX.status"))]
lgg = subset(lgg, !(is.na(IDH.status))) #lost 3 patients 

colnames(lgg)[c(7:8, 11)] = c("OS", "OS.time", "ENSG00000253187")
lgg$exp_log1p = log1p(lgg$ENSG00000253187)

z1 = which(lgg$patient %in% maf$patient[maf$idh_mut=="yes"])
z2 = which(lgg$patient %in% maf$patient[maf$Hugo_Symbol == "TP53"])

lgg$tp53 = ""
lgg$tp53[z2] = "TP53mut"
lgg$tp53[-z2] = "TP53WT"

pdf("/u/kisaev/LGG_HOXA10AS_with_IDH_1pq19codel.pdf", width=15)

ggboxplot(lgg, x="IDH.status", y="exp_log1p", fill="IDH.status") + stat_compare_means() + stat_n_text()
ggboxplot(lgg, x="tp53", y="exp_log1p", fill="tp53") + stat_compare_means() + stat_n_text()
ggboxplot(lgg, x="X1p.19q.codeletion", y="exp_log1p", fill="X1p.19q.codeletion") + stat_compare_means() + stat_n_text()
ggboxplot(lgg, x="TERT.promoter.status", y="exp_log1p", fill="TERT.promoter.status") + stat_compare_means() + stat_n_text()

lgg = subset(lgg, IDH.status == "Mutant")

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
  lgg$med = factor(lgg$med, levels = c("High", "Low"))
  lgg$X1p.19q.codeletion = factor(lgg$X1p.19q.codeletion, levels = c("non-codel", "codel"))

  lgg$OS.time = lgg$OS.time/365
  
  #KM plot
  fit <- survfit(Surv(OS.time, OS) ~ med, data = lgg)

  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = lgg)
  pval = round(tidy(res.cox)$p.value, digits=6)

  #HOXA10-as in LGG IDH mutants 
  s <- ggsurvplot(
            title = paste("LGG IDH mutants", "HOXA10-AS expression\n",
              "pval=", pval),
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

  ggboxplot(lgg, x="med", y="ENSG00000253187", fill="med", title="IDH Mutants LGGs FPKM-UQ expression")+
  stat_compare_means() + stat_n_text()
  
  #IDH mutants with or without 1pq19 codel 
  fit <- survfit(Surv(OS.time, OS) ~ X1p.19q.codeletion , data = lgg)
  res.cox <- coxph(Surv(OS.time, OS) ~ X1p.19q.codeletion, data = lgg)
  pval = round(tidy(res.cox)$p.value, digits=6)

  s <- ggsurvplot(
            title = paste("LGG IDH mutants", "X1p.19q.codeletion\n", "pval=", pval),
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


  #1p19q codel + HOXA10as 
  fit <- survfit(Surv(OS.time, OS) ~ med + X1p.19q.codeletion , data = lgg)
  res.cox <- coxph(Surv(OS.time, OS) ~ med + X1p.19q.codeletion, data = lgg)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants", "HOXA10AS+X1p.19q.codeletion\n", "pval=", pval),
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


  #take only IDH mutants + 1pq19 codel stratify by hoxa10as 
  lgg_idh_m_1pq19 = as.data.table(filter(lgg, X1p.19q.codeletion == "codel"))
  fit <- survfit(Surv(OS.time, OS) ~ med , data = lgg_idh_m_1pq19)
  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = lgg_idh_m_1pq19)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants with 1p19q codel", "HOXA10-AS expression\n",
              "pval=", pval),
            fit, 
            xlab = "Time (Years)", 
            #surv.median.line = "hv",
            font.main = c(14, "bold", "black"),
            font.x = c(12, "plain", "black"),
            font.y = c(12, "plain", "black"),
            font.tickslab = c(11, "plain", "black"),
            font.legend = 10,
            risk.table.fontsize = 5, 
            data = lgg_idh_m_1pq19,      # data used to fit survival curves. 
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


  #take only IDH mutants + 1pq19 codel stratify by hoxa10as 
  lgg_idh_m_n1pq19 = as.data.table(filter(lgg, X1p.19q.codeletion == "non-codel"))
  fit <- survfit(Surv(OS.time, OS) ~ med , data = lgg_idh_m_n1pq19)
  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = lgg_idh_m_n1pq19)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants with no 1pq19 codel", "HOXA10-AS expression\n",
              "pval=", pval),
            fit, 
            xlab = "Time (Years)", 
            #surv.median.line = "hv",
            font.main = c(14, "bold", "black"),
            font.x = c(12, "plain", "black"),
            font.y = c(12, "plain", "black"),
            font.tickslab = c(11, "plain", "black"),
            font.legend = 10,
            risk.table.fontsize = 5, 
            data = lgg_idh_m_n1pq19,      # data used to fit survival curves. 
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

  #take LGG IDH mutants and look at TP53+HOXA10A
  lgg_idh_m_tp53_mut = as.data.table(filter(lgg, tp53 == "TP53mut"))
  fit <- survfit(Surv(OS.time, OS) ~ med , data = lgg_idh_m_tp53_mut)
  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = lgg_idh_m_tp53_mut)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants with TP53 mutation", "HOXA10-AS expression\n",
              "pval=", pval),
            fit, 
            xlab = "Time (Years)", 
            #surv.median.line = "hv",
            font.main = c(14, "bold", "black"),
            font.x = c(12, "plain", "black"),
            font.y = c(12, "plain", "black"),
            font.tickslab = c(11, "plain", "black"),
            font.legend = 10,
            risk.table.fontsize = 5, 
            data = lgg_idh_m_tp53_mut,      # data used to fit survival curves. 
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


  #take LGG IDH mutants and look at TP53+HOXA10A
  lgg_idh_m_tp53_wt = as.data.table(filter(lgg, tp53 == "TP53WT"))
  fit <- survfit(Surv(OS.time, OS) ~ med , data = lgg_idh_m_tp53_wt)
  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = lgg_idh_m_tp53_wt)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants with TP53 WT", "HOXA10-AS expression\n",
              "pval=", pval),
            fit, 
            xlab = "Time (Years)", 
            #surv.median.line = "hv",
            font.main = c(14, "bold", "black"),
            font.x = c(12, "plain", "black"),
            font.y = c(12, "plain", "black"),
            font.tickslab = c(11, "plain", "black"),
            font.legend = 10,
            risk.table.fontsize = 5, 
            data = lgg_idh_m_tp53_wt,      # data used to fit survival curves. 
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


  #take LGG IDH mutants and look at TP53+HOXA10A
  fit <- survfit(Surv(OS.time, OS) ~ med + tp53 , data = lgg)
  res.cox <- coxph(Surv(OS.time, OS) ~ med + tp53, data = lgg)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants", "HOXA10-AS vs TP53\n",
              "pval=", pval),
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


  #take LGG IDH mutants and look at just hOXA10-AS high and check TP53 mutant vs WT 
  lnc_high_tp53 = as.data.table(filter(lgg, med == "High"))
  fit <- survfit(Surv(OS.time, OS) ~ tp53 , data = lnc_high_tp53)
  res.cox <- coxph(Surv(OS.time, OS) ~ tp53, data = lnc_high_tp53)
  pval = round(tidy(res.cox)$p.value, digits=6)[1]

  s <- ggsurvplot(
            title = paste("LGG IDH mutants", "HOXA10-AS high only vs TP53\n",
              "pval=", pval),
            fit, 
            xlab = "Time (Years)", 
            #surv.median.line = "hv",
            font.main = c(14, "bold", "black"),
            font.x = c(12, "plain", "black"),
            font.y = c(12, "plain", "black"),
            font.tickslab = c(11, "plain", "black"),
            font.legend = 10,
            risk.table.fontsize = 5, 
            data = lnc_high_tp53,      # data used to fit survival curves. 
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

