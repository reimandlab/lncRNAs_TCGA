source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#library(M3C)
#library(umap)

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

#transcriptome subtype
tum_subtype = lgg_wsurv
z= which(is.na(tum_subtype$Transcriptome.Subtype))
tum_subtype = tum_subtype[-z,]
tum_subtype$OS.time = tum_subtype$OS.time/365
tum_subtype$TS = tum_subtype$Transcriptome.Subtype

pdf("/u/kisaev/HOXA10AS_transcriptome_subtype_KM_plot.pdf")
fit <- survfit(Surv(OS.time, OS) ~ TS, data = tum_subtype)
s1 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = tum_subtype,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        xlim = c(0,10),
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
        xlim = c(0,10),
        palette =c("#4DBBD5FF", "#E64B35FF"),
        break.time.by = 1)#,     # break X axis in time intervals by 500.
print(s2)
dev.off()

summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="CL")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="ME")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="NE")))$coefficients
summary(coxph(Surv(OS.time, OS) ~ HOXA10AS, data = filter(tum_subtype, TS=="PN")))$coefficients


##treatment outcome
tum_outcome = lgg_wsurv
z= which(tum_outcome$treatment_outcome_first_course %in% c("Complete Remission/Response", "Partial Remission/Response",
"Progressive Disease", "Stable Disease"))
tum_outcome = tum_outcome[z,]
tum_outcome$OS.time = tum_outcome$OS.time/365
tum_outcome$TOFC = tum_outcome$treatment_outcome_first_course
tum_outcome$TOFC[tum_outcome$TOFC == "Complete Remission/Response"] = "CR"
tum_outcome$TOFC[tum_outcome$TOFC == "Partial Remission/Response"] = "PR"
tum_outcome$TOFC[tum_outcome$TOFC == "Progressive Disease"] = "PD"
tum_outcome$TOFC[tum_outcome$TOFC == "Stable Disease"] = "SD"
tum_outcome$HOXA10AS = factor(tum_outcome$HOXA10AS, levels=c("Low HOXA10-AS", "High HOXA10-AS"))


pdf("/u/kisaev/HOXA10AS_treatment_outcome_KM_plot.pdf")

fit <- survfit(Surv(OS.time, OS) ~ TOFC, data = tum_outcome)
s1 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        xlim = c(0,10),
        risk.table = TRUE,      # present narrower X axis, but not affect
        break.time.by = 1)#,     # break X axis in time intervals by 500.
        #palette = c("purple", "orange"))
print(s1)

fit <- survfit(Surv(OS.time, OS) ~ HOXA10AS, data = tum_outcome)
s2 <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = tum_outcome,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        xlim = c(0,10),
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
