source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

get_median_risk_group = function(PI_centered, PI_lower_thres, epsilon = 1e-16) {

  if (!any(PI_centered > PI_lower_thres)) {
    PI_lower_thres = PI_lower_thres - epsilon
  } else if (all(PI_centered > PI_lower_thres)) {
    PI_lower_thres = PI_lower_thres + epsilon
  }

  risk_group = 1 + (PI_centered > PI_lower_thres)
  risk_group = c("low_risk", "high_risk")[risk_group]
  risk_group = factor(risk_group, levels = c("low_risk", "high_risk"))
  risk_group
}

#data read in

lgg=readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")
lgg_wsurv = unique(lgg[,c("patient", "ENSG00000253187",
"Transcriptome.Subtype", "treatment_outcome_first_course",
"IDH.status", "OS", "OS.time")])

lgg_rna = as.data.table(filter(rna, type =="LGG"))

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

lgg_wsurv$HOXA10AS = ""
lgg_wsurv$HOXA10AS[lgg_wsurv$ENSG00000253187 == 0 ] = "Low HOXA10-AS"
lgg_wsurv$HOXA10AS[lgg_wsurv$ENSG00000253187 > 0 ] = "High HOXA10-AS"

dat_save_figure = lgg_wsurv
dat_save_figure$OS.time = dat_save_figure$OS.time/365

#add censoring to match other panels in figure 5
dat_save_figure$OS[dat_save_figure$OS.time >8] = 0
dat_save_figure$OS.time[dat_save_figure$OS.time >8] = 8

dat_save_figure$TS = dat_save_figure$Transcriptome.Subtype
dat_save_figure$TOFC = dat_save_figure$treatment_outcome_first_course
dat_save_figure$TOFC[dat_save_figure$TOFC == "Complete Remission/Response"] = "CR"
dat_save_figure$TOFC[dat_save_figure$TOFC == "Partial Remission/Response"] = "PR"
dat_save_figure$TOFC[dat_save_figure$TOFC == "Progressive Disease"] = "PD"
dat_save_figure$TOFC[dat_save_figure$TOFC == "Stable Disease"] = "SD"
dat_save_figure$TOFC[dat_save_figure$TOFC %in% c("[Not Available]", "[Unknown]",
"[Discrepancy]", "[Not Applicable]")] = "NA"
dat_save_figure = dat_save_figure[,c("patient", "ENSG00000253187",
"HOXA10AS", "TOFC", "TS",
"OS", "OS.time")]
write.csv(dat_save_figure, "/u/kisaev/data_for_figure5_TCGA_TOFC_TS.csv", quote=F, row.names=F)

#transcriptome subtype
tum_subtype = lgg_wsurv
z= which(is.na(tum_subtype$Transcriptome.Subtype))
tum_subtype = tum_subtype[-z,]
tum_subtype$OS.time = tum_subtype$OS.time/365
tum_subtype$TS = tum_subtype$Transcriptome.Subtype
#add censoring to match other panels in figure 5
tum_subtype$OS[tum_subtype$OS.time >8] = 0
tum_subtype$OS.time[tum_subtype$OS.time >8] = 8

pdf("/u/kisaev/HOXA10AS_transcriptome_subtype_KM_plot.pdf")
fit <- survfit(Surv(OS.time, OS) ~ TS, data = tum_subtype)
s1 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = tum_subtype,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        risk.table = TRUE,      # present narrower X axis, but not affect
        break.time.by = 1)#,     # break X axis in time intervals by 500.
        #palette = c("purple", "orange"))
print(s1)
tum_subtype$HOXA10AS = factor(tum_subtype$HOXA10AS, levels=c("Low HOXA10-AS", "High HOXA10-AS"))
fit <- survfit(Surv(OS.time, OS) ~ HOXA10AS, data = tum_subtype)
s2 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = tum_subtype,      # data used to fit survival curves.
        facet.by = "TS",
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        palette =c("#4DBBD5FF", "#E64B35FF"),
        break.time.by = 1)#,     # break X axis in time intervals by 500.
print(s2)
dev.off()

summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="CL")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="ME")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="NE")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="PN")))$coefficients

clin_only = coxph(Surv(OS.time, OS) ~ TS, data =tum_subtype)
clin_plus_lnc = coxph(Surv(OS.time, OS) ~ TS + HOXA10AS, data =tum_subtype)
anova(clin_only, clin_plus_lnc) #pvalue = 5.58e-06

xx = table(tum_subtype$HOXA10AS,tum_subtype$TS)
chisq.test(xx) #pvalue = 8.006e-14

#get risk based models clinical only
cox_lnc = coxph(Surv(OS.time, OS) ~ TS, data = tum_subtype)
relRisk <- predict(cox_lnc, tum_subtype, type="risk")   # relative risk
tum_subtype$rel_risk_clin_only = relRisk

# split into two risk groups based on median
PI_lower_thres = median(tum_subtype$rel_risk_clin_only)
PI_max_threshold = summary(tum_subtype$rel_risk_clin_only)[5]
tum_subtype$risk_group_clin_only = get_median_risk_group(
  tum_subtype$rel_risk_clin_only,
  PI_lower_thres)
tum_subtype$risk_group_clin_only = factor(tum_subtype$risk_group_clin_only, levels=c("low_risk", "high_risk"))

#get risk based models clinical plus lncRNA only
cox_lnc = coxph(Surv(OS.time, OS) ~ TS + HOXA10AS, data = tum_subtype)
relRisk <- predict(cox_lnc, tum_subtype, type="risk")   # relative risk
tum_subtype$rel_risk_clin_plus_lncRNA = relRisk

# split into two risk groups based on median
PI_lower_thres = median(tum_subtype$rel_risk_clin_plus_lncRNA)
PI_max_threshold = summary(tum_subtype$rel_risk_clin_plus_lncRNA)[5]
tum_subtype$risk_group_clin_plus_lncRNA = get_median_risk_group(
  tum_subtype$rel_risk_clin_plus_lncRNA,
  PI_lower_thres)
tum_subtype$risk_group_clin_plus_lncRNA = factor(tum_subtype$risk_group_clin_plus_lncRNA, levels=c("low_risk", "high_risk"))

pdf("/u/kisaev/HOXA10AS_transcriptome_subtype_risk_based_KM_plot.pdf", width=6, height=5)
fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_only, data = tum_subtype)
cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_only, data = tum_subtype)
hr=round(summary(cox_lnc)$coefficients[2], digits=3)
pval=round(summary(cox_lnc)$coefficients[5], digits=5)
lowrisksamps = table(tum_subtype$risk_group_clin_only)[1]
highrisksamps = table(tum_subtype$risk_group_clin_only)[2]

s1 <- ggsurvplot(fit ,
        title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
        xlab = "Time (Years)",
        data = tum_subtype,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        risk.table = FALSE,      # present narrower X axis, but not affect
        break.time.by = 1,     # break X axis in time intervals by 500.
        palette =c("#4DBBD5FF", "#E64B35FF"))
print(s1)
fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = tum_subtype)

cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = tum_subtype)
hr=round(summary(cox_lnc)$coefficients[2], digits=3)
pval=round(summary(cox_lnc)$coefficients[5], digits=14)
lowrisksamps = table(tum_subtype$risk_group_clin_plus_lncRNA)[1]
highrisksamps = table(tum_subtype$risk_group_clin_plus_lncRNA)[2]

s2 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
        data = tum_subtype,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        palette =c("#4DBBD5FF", "#E64B35FF"),
        break.time.by = 1)#,     # break X axis in time intervals by 500.
print(s2)
dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##treatment outcome
tum_outcome = lgg_wsurv
z= which(tum_outcome$treatment_outcome_first_course %in% c("Complete Remission/Response", "Partial Remission/Response",
"Progressive Disease", "Stable Disease"))
tum_outcome = tum_outcome[z,]
tum_outcome$OS.time = tum_outcome$OS.time/365
#add censoring to match other panels in figure 5
tum_outcome$OS[tum_outcome$OS.time >8] = 0
tum_outcome$OS.time[tum_outcome$OS.time >8] = 8
tum_outcome$TOFC = tum_outcome$treatment_outcome_first_course
tum_outcome$TOFC[tum_outcome$TOFC == "Complete Remission/Response"] = "CR"
tum_outcome$TOFC[tum_outcome$TOFC == "Partial Remission/Response"] = "PR"
tum_outcome$TOFC[tum_outcome$TOFC == "Progressive Disease"] = "PD"
tum_outcome$TOFC[tum_outcome$TOFC == "Stable Disease"] = "SD"
tum_outcome$HOXA10AS = factor(tum_outcome$HOXA10AS, levels=c("Low HOXA10-AS", "High HOXA10-AS"))

pdf("/u/kisaev/HOXA10AS_treatment_outcome_KM_plot.pdf")

fit <- survfit(Surv(OS.time, OS) ~ TOFC, data = tum_outcome)
s1 <- ggsurvplot(fit ,
      font.main = c(7, "bold", "black"),
        xlab = "Time (Years)",
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        risk.table = TRUE,      # present narrower X axis, but not affect
        break.time.by = 1)#,     # break X axis in time intervals by 500.
        #palette = c("purple", "orange"))
print(s1)

fit <- survfit(Surv(OS.time, OS) ~ HOXA10AS, data = tum_outcome)
s2 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        font.main = c(7, "bold", "black"),
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        palette =c("#4DBBD5FF", "#E64B35FF"),
        facet.by = "TOFC",
        break.time.by = 1)#,     # break X axis in time intervals by 500.
        #palette = c("purple", "orange"))
print(s2)

dev.off()

#get hazard ratios and p-values
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_outcome, TOFC=="CR")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_outcome, TOFC=="PR")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_outcome, TOFC=="PD")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_outcome, TOFC=="SD")))$coefficients

clin_only = coxph(Surv(OS.time, OS) ~ TOFC, data =tum_outcome)
clin_plus_lnc = coxph(Surv(OS.time, OS) ~ TOFC + HOXA10AS, data =tum_outcome)
anova(clin_only, clin_plus_lnc) #pvalue = 1.567e-05

xx = table(tum_outcome$TOFC, tum_outcome$HOXA10AS)
chisq.test(xx) #pvalue = 0.0001013

#get risk based models clinical only
cox_lnc = coxph(Surv(OS.time, OS) ~ TOFC, data = tum_outcome)
relRisk <- predict(cox_lnc, tum_outcome, type="risk")   # relative risk
tum_outcome$rel_risk_clin_only = relRisk

# split into two risk groups based on median
PI_lower_thres = median(tum_outcome$rel_risk_clin_only)
PI_max_threshold = summary(tum_outcome$rel_risk_clin_only)[5]
tum_outcome$risk_group_clin_only = get_median_risk_group(
  tum_outcome$rel_risk_clin_only,
  PI_lower_thres)
tum_outcome$risk_group_clin_only = factor(tum_outcome$risk_group_clin_only, levels=c("low_risk", "high_risk"))

#get risk based models clinical plus lncRNA only
cox_lnc = coxph(Surv(OS.time, OS) ~ TOFC + HOXA10AS, data = tum_outcome)
relRisk <- predict(cox_lnc, tum_outcome, type="risk")   # relative risk
tum_outcome$rel_risk_clin_plus_lncRNA = relRisk

# split into two risk groups based on median
PI_lower_thres = median(tum_outcome$rel_risk_clin_plus_lncRNA)
PI_max_threshold = summary(tum_outcome$rel_risk_clin_plus_lncRNA)[5]
tum_outcome$risk_group_clin_plus_lncRNA = get_median_risk_group(
  tum_outcome$rel_risk_clin_plus_lncRNA,
  PI_lower_thres)
tum_outcome$risk_group_clin_plus_lncRNA = factor(tum_outcome$risk_group_clin_plus_lncRNA, levels=c("low_risk", "high_risk"))

pdf("/u/kisaev/HOXA10AS_TOFC_risk_based_KM_plot.pdf", width=6, height=5)
fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_only, data = tum_outcome)
cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_only, data = tum_outcome)
hr=round(summary(cox_lnc)$coefficients[2], digits=3)
pval=round(summary(cox_lnc)$coefficients[5], digits=15)
lowrisksamps = table(tum_outcome$risk_group_clin_only)[1]
highrisksamps = table(tum_outcome$risk_group_clin_only)[2]

s1 <- ggsurvplot(fit ,
        title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
        xlab = "Time (Years)",
        font.main = c(7, "bold", "black"),
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        risk.table = FALSE,      # present narrower X axis, but not affect
        break.time.by = 1,     # break X axis in time intervals by 500.
        palette =c("#4DBBD5FF", "#E64B35FF"))
print(s1)
fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = tum_outcome)

cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = tum_outcome)
hr=round(summary(cox_lnc)$coefficients[2], digits=3)
pval=round(summary(cox_lnc)$coefficients[5], digits=15)
lowrisksamps = table(tum_outcome$risk_group_clin_plus_lncRNA)[1]
highrisksamps = table(tum_outcome$risk_group_clin_plus_lncRNA)[2]

s2 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        font.main = c(7, "bold", "black"),
        title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        #xlim = c(0,8),
        palette =c("#4DBBD5FF", "#E64B35FF"),
        break.time.by = 1)#,     # break X axis in time intervals by 500.
print(s2)
dev.off()
