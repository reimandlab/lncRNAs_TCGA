source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#library(M3C)
#library(umap)

lgg=readRDS("TCGA_gbm_wsubtype_info_biolinks.rds")

gbm_exp = filter(rna, type=="GBM")
z = which(colnames(gbm_exp) %in% c("patient", "ENSG00000253187"))
gbm_exp =gbm_exp[,..z]

lgg_wsurv = unique(lgg[,c("patient",
"Transcriptome.Subtype", "treatment_outcome_first_course",
"IDH.status", "OS", "OS.time")])
dim(lgg_wsurv)
lgg_wsurv = merge(lgg_wsurv, gbm_exp, by="patient")

lgg_rna = as.data.table(filter(rna, type =="GBM"))
dim(lgg_rna)

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

pdf("/u/kisaev/HOXA10AS_transcriptome_subtype_KM_plot_GBM.pdf")
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#using median

lgg_wsurv$HOXA10AS = ""
med = median(lgg_wsurv$ENSG00000253187)
lgg_wsurv$HOXA10AS[lgg_wsurv$ENSG00000253187 < med ] = "Low HOXA10-AS"
lgg_wsurv$HOXA10AS[lgg_wsurv$ENSG00000253187 >= med ] = "High HOXA10-AS"

#transcriptome subtype
tum_subtype = lgg_wsurv
z= which(is.na(tum_subtype$Transcriptome.Subtype))
tum_subtype = tum_subtype[-z,]
tum_subtype$OS.time = tum_subtype$OS.time/365
tum_subtype$TS = tum_subtype$Transcriptome.Subtype

pdf("/u/kisaev/HOXA10AS_transcriptome_subtype_KM_plot_GBM_via_median_dichom.pdf")
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

#Overall GBM survival curve
dat = lgg_rna
dat$OS = as.numeric(dat$OS)
dat$OS.time = as.numeric(dat$OS.time)
dat$OS.time = dat$OS.time/365
dat$gene = dat$ENSG00000253187
dat$gene_zero = ""
dat$gene_zero[dat$ENSG00000253187 == 0 ] = "Low HOXA10-AS"
dat$gene_zero[dat$ENSG00000253187 > 0 ] = "High HOXA10-AS"

dat$gene_med = ""
med = median(dat$ENSG00000253187)
dat$gene_med[dat$ENSG00000253187 < med ] = "Low HOXA10-AS"
dat$gene_med[dat$ENSG00000253187 >= med ] = "High HOXA10-AS"

dat$gene_zero = factor(dat$gene_zero, levels=c("Low HOXA10-AS","High HOXA10-AS"))
dat$gene_med = factor(dat$gene_med, levels=c("Low HOXA10-AS","High HOXA10-AS"))

cox_mod = coxph(Surv(OS.time, OS) ~ gene_zero, data = dat)
hr = summary(cox_mod)$coefficients[2]
dtt = dat

fit <- survfit(Surv(OS.time, OS) ~ gene_zero, data = dtt)

s1 <- ggsurvplot(title = paste("HOXA10-AS in GBM zero dicho", "HR =", round(hr, digits=2)),
          fit,
          ylab = "Survival Probability" ,
          xlab = "Time (Years)",
          data = dtt,      # data used to fit survival curves.
          risk.table = FALSE,       # show risk table.
          legend = "bottom",
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for
          xlim = c(0,10),        # present narrower X axis, but not affect
          break.time.by = 1,     # break X axis in time intervals by 500.
          palette = c("#4DBBD5FF", "#E64B35FF"))

cox_mod = coxph(Surv(OS.time, OS) ~ gene_med, data = dat)
hr = summary(cox_mod)$coefficients[2]
dtt = dat
fit <- survfit(Surv(OS.time, OS) ~ gene_med, data = dtt)

s2 <- ggsurvplot(title = paste("HOXA10-AS in GBM median dicho", "HR =", round(hr, digits=2)),
                    fit,
                    ylab = "Survival Probability" ,
                    xlab = "Time (Years)",
                    #surv.median.line = "hv",
                    #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
                    data = dtt,      # data used to fit survival curves.
                    risk.table = FALSE,       # show risk table.
                    legend = "bottom",
                    pval = TRUE,             # show p-value of log-rank test.
                    conf.int = FALSE,        # show confidence intervals for
                    xlim = c(0,10),        # present narrower X axis, but not affect
                    break.time.by = 1,     # break X axis in time intervals by 500.
                    palette = c("#4DBBD5FF", "#E64B35FF"))

pdf("/u/kisaev/HOXA10AS_OS_KM_plot_GBM.pdf")
print(s1)
print(s2)
dev.off()
