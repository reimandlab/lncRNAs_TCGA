set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files

#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs

res = allCands
res_lncs = subset(res, cancer == "Liver hepatocellular carcinoma" )$gene

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

#ENSG00000264026

#DATA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===

#TCGA+
tcga =readRDS("/u/kisaev/LIHC_TCGA_data_2021.rds")
tcga=as.data.frame(tcga)
#z = which(!(colnames(tcga) %in% c(res_lncs, "patient", "OS", "OS.time")))
#tcga = tcga[,-z]
head(tcga)

lncs = which(str_detect(colnames(tcga), "ENSG"))

tcga_full = tcga

#label lncs as high versus low
for(i in lncs){print(i)
  print(i)
  print(colnames(tcga)[i])
  med = median(unlist(tcga[,i]))
  if(med == 0){
    z1 = which(unlist(tcga[,i]) > 0)
    z2 = which(unlist(tcga[,i]) == 0)
    tcga[z1,i] = 1
    tcga[z2,i] = 0
  }

  if(med > 0){
    z1 = which(unlist(tcga[,i]) >= med)
    z2 = which(unlist(tcga[,i]) < med)
    tcga[z1,i] = 1
    tcga[z2,i] = 0
  }
  table(tcga[,i])
  tcga[,i] = factor(tcga[,i], levels=c(0,1))
}

#PCAWG+
pcawg =readRDS("/u/kisaev/LIHC_PCAWG_data.rds")
pcawg=as.data.frame(pcawg)
z = which(!(colnames(pcawg) %in% c(res_lncs, "patient", "status", "time")))
pcawg = pcawg[,-z]
head(pcawg)

lncs = which(str_detect(colnames(pcawg), "ENSG"))
pcawg_full = pcawg

#label lncs as high versus low
del = c()

for(i in lncs){print(i)
  print(i)
  print(colnames(pcawg)[i])
  med = median(unlist(pcawg[,i]))
  if(med == 0){
    z1 = which(unlist(pcawg[,i]) > 0)
    z2 = which(unlist(pcawg[,i]) == 0)
    check = median(unlist(pcawg[z1,i]))
    if((check >= 0.001) & (length(z1) >= 5)){
      pcawg[z1,i] = 1
      pcawg[z2,i] = 0
    }
    if((!(check >= 0.001)) | is.na(check)){
      del = c(del, colnames(pcawg)[i])
    }
    if(!((check >= 0.001) & (length(z1) >= 5))){
      del = c(del, colnames(pcawg)[i])
    }
  }

  if(med > 0){
    if(med >= 0.001){
    z1 = which(unlist(pcawg[,i]) >= med)
    z2 = which(unlist(pcawg[,i]) < med)
    pcawg[z1,i] = 1
    pcawg[z2,i] = 0
    }

    if(med < 0.001){
      del = c(del, colnames(pcawg)[i])
    }
    }
  #table(pcawg[,i])
  pcawg[,i] = factor(pcawg[,i], levels=c(0,1))
}

z = which(colnames(pcawg) %in% del)
if(!(length(z)==0)){
  pcawg = pcawg[,-z]}

lncs_keep = c("ENSG00000178977", "ENSG00000264026")

#make labels survival the same in both datasets
colnames(pcawg)[which(colnames(pcawg)=="status")] = "OS"
colnames(pcawg)[which(colnames(pcawg)=="time")] = "OS.time"
pcawg$OS[pcawg$OS == "deceased"]= 1
pcawg$OS[pcawg$OS == "alive"]= 0
pcawg$OS.time = as.numeric(pcawg$OS.time)/365
tcga$OS.time = as.numeric(tcga$OS.time)/365

rownames(pcawg) = pcawg$patient
pcawg$patient = NULL
pcawg$canc = NULL
pcawg$histo = NULL
pcawg$sex = NULL
pcawg$donor_age_at_diagnosis = NULL

rownames(tcga) = tcga$patient
tcga$patient = NULL

#make column orders same two datasets
tcga = tcga[,c("ENSG00000178977", "ENSG00000264026", "OS", "OS.time")]
pcawg = pcawg[,colnames(tcga)]

#MODEL+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===

#cols we want to use to build model

# train a multivariate model to build prognostic index (PI)
# keep training coefficients because these will be used for validation data
cox_lnc = coxph(Surv(OS.time, OS) ~ ., data = tcga)
coef_training = coef(cox_lnc)

mean_lnc1  <- mean(as.numeric(tcga$ENSG00000178977)-1)  # average of financial aid dummy
mean_lnc2  <- mean(as.numeric(tcga$ENSG00000264026)-1)  # average of financial aid dummy

rMean <- exp(coef_training[[1]] * mean_lnc1 + coef_training[[2]]  * mean_lnc2)

r1234 <- exp(coef_training[[1]] * (as.numeric(tcga[1:6, "ENSG00000178977"])-1) +
               coef_training[[2]] * (as.numeric(tcga[1:6, "ENSG00000264026"])-1))

#relative risk
r1234 / rMean
relRisk <- predict(cox_lnc, tcga, type="risk")   # relative risk
relRisk[1:6]
tcga$rel_risk = relRisk

# split into two risk groups based on median
PI_lower_thres = median(tcga$rel_risk)
PI_max_threshold = summary(tcga$rel_risk)[5]
tcga$risk_group = get_median_risk_group(
  tcga$rel_risk,
  PI_lower_thres)
tcga$risk_group = factor(tcga$risk_group, levels=c("low_risk", "high_risk"))

#tcga$risk_group=""
#tcga$rel_risk = as.factor(tcga$rel_risk)
#tcga$risk_group[tcga$rel_risk==2.41070925571496] = "high_risk"
#tcga$risk_group[tcga$rel_risk==1.16903890886235] = "int_risk"
#tcga$risk_group[tcga$rel_risk==1.11588407304556 ] = "int_risk"
#tcga$risk_group[tcga$rel_risk==0.541131990959718] = "low_risk"
#table(tcga$risk_group)
#tcga$risk_group = factor(tcga$risk_group, levels=c("high_risk", "int_risk", "low_risk"))

# train a model on training data, based on risk groups
h1_train = coxph(Surv(OS.time, OS) ~ risk_group, data = tcga)

#apply to validation data
pcawg_vals <- exp(#coef_training[[1]] * (as.numeric(pcawg[1:42, "ENSG00000231918"])-1)
             #+ coef_training[[2]] * (as.numeric(pcawg[1:42, "ENSG00000236760"])-1)
            # + coef_training[[1]] * (as.numeric(pcawg[1:42, "ENSG00000236437"])-1)
              coef_training[[1]] * (as.numeric(pcawg[1:42, "ENSG00000178977"])-1)
             + coef_training[[2]] * (as.numeric(pcawg[1:42, "ENSG00000264026"])-1))
           #  + coef_training[[2]] * (as.numeric(pcawg[1:42, "ENSG00000230490"])-1))
             #+ coef_training[[7]] * (as.numeric(pcawg[1:42, "ENSG00000237027"])-1)
           #  + coef_training[[4]] * (as.numeric(pcawg[1:42, "ENSG00000258731"])-1))

relRisk_pcawg = pcawg_vals / rMean
pcawg$rel_risk = relRisk_pcawg

pcawg$risk_group = get_median_risk_group(
  relRisk_pcawg,
  PI_lower_thres)

#pcawg$rel_risk = relRisk_pcawg
#pcawg$risk_group=""
#pcawg$rel_risk = as.factor(pcawg$rel_risk)
#pcawg$risk_group[pcawg$rel_risk==2.41070925571496] = "high_risk"
#pcawg$risk_group[pcawg$rel_risk==1.16903890886235] = "int_risk"
#pcawg$risk_group[pcawg$rel_risk==1.11588407304556] = "int_risk"
#pcawg$risk_group[pcawg$rel_risk==0.541131990959718] = "low_risk"
#table(pcawg$risk_group)

#pcawg$risk_group = factor(pcawg$risk_group, levels=c("high_risk", "int_risk", "low_risk"))
pcawg$risk_group = factor(pcawg$risk_group, levels=c("low_risk", "high_risk"))

pcawg$OS = as.numeric(pcawg$OS)

# train a model on validation data, based on risk groups
h1_valid = coxph(Surv(OS.time, OS) ~ risk_group, data = pcawg)
h1_valid

#make KM plots for TCGA and Shanghai
tcga$risk_group = factor(tcga$risk_group, levels=c("high_risk", "low_risk"))
pcawg$risk_group = factor(pcawg$risk_group, levels=c("high_risk", "low_risk"))

tcga_hr=round(summary(h1_train)$coefficients[2], digits=2)
tcga_pval = round(summary(h1_train)$coefficients[5], digits=8)

pcawg_hr=round(summary(h1_valid)$coefficients[2], digits=2)
pcawg_pval = round(summary(h1_valid)$coefficients[5], digits=4)

fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = tcga)
t1 <- ggsurvplot(
  title=paste("lncRNAs risk-score (TCGA)", "HR=", tcga_hr, "Wald p=", tcga_pval),
    fit,
  xlab = "Time (Years)",
#  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = tcga,      # data used to fit survival curves.
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
#print(t1)


fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = pcawg)
v1 <- ggsurvplot(
  title=paste("lncRNAs risk-score (PCAWG)", "HR=", pcawg_hr, "Wald p=", pcawg_pval),
  fit,
  xlab = "Time (Years)",
#  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = pcawg,      # data used to fit survival curves.
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
#print(v1)

pdf("LIHC_case_study_lncRNA_candidates_combined_risk_model_two_lncRNAs_2021.pdf")
print(t1)
print(v1)
dev.off()

#plots for just individual lncsRNAs

#plots for just individual lncsRNAs
tcga$ENSG00000178977 = factor(tcga$ENSG00000178977, levels=c(0,1))
fit <- survfit(Surv(OS.time, OS) ~ ENSG00000178977, data = tcga)
mo = coxph(Surv(OS.time, OS) ~ ENSG00000178977, data = tcga)
tcga_hr=round(summary(mo)$coefficients[2], digits=2)
tcga_pval = round(summary(mo)$coefficients[5], digits=6)
t1 <- ggsurvplot(
  title=paste("LINC00324 (TCGA)", "HR=", tcga_hr, "Wald p=", tcga_pval),
  fit,
  xlab = "Time (Years)",
  #  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = tcga,      # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  legend = "right",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for
  # point estimaes of survival curves.
  xlim = c(0,10),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #palette = colorRampPalette(mypal)(14),
  palette = c("#4DBBD5FF", "#E64B35FF"),
  #palette = "npg",
  #ggtheme = theme_bw(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

pcawg$ENSG00000178977 = factor(pcawg$ENSG00000178977, levels=c(0,1))
fit <- survfit(Surv(OS.time, OS) ~ ENSG00000178977, data = pcawg)
mo = coxph(Surv(OS.time, OS) ~ ENSG00000178977, data = pcawg)
pcawg_hr=round(summary(mo)$coefficients[2], digits=2)
pcawg_pval = round(summary(mo)$coefficients[5], digits=2)
v1 <- ggsurvplot(
  title=paste("LINC00324 (PCAWG)", "HR=", pcawg_hr, "Wald p=", pcawg_pval),
  fit,
  xlab = "Time (Years)",
  #  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = pcawg,      # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  legend = "right",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for
  # point estimaes of survival curves.
  xlim = c(0,10),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #palette = colorRampPalette(mypal)(14),
  palette = c("#4DBBD5FF", "#E64B35FF"),
  #palette = "npg",
  #ggtheme = theme_bw(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)
pdf("LIHC_TCGA_PCAWG_LINC00324.pdf")
print(t1)
print(v1)
dev.off()

#ENSG00000264026
tcga$ENSG00000264026 = factor(tcga$ENSG00000264026, levels=c(0,1))
fit <- survfit(Surv(OS.time, OS) ~ ENSG00000264026, data = tcga)
mo = coxph(Surv(OS.time, OS) ~ ENSG00000264026, data = tcga)
tcga_hr=round(summary(mo)$coefficients[2], digits=2)
tcga_pval = round(summary(mo)$coefficients[5], digits=6)
t1 <- ggsurvplot(
  title=paste("LINC02003 (TCGA)", "HR=", tcga_hr, "Wald p=", tcga_pval),
  fit,
  xlab = "Time (Years)",
  #  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = tcga,      # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  legend = "right",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for
  # point estimaes of survival curves.
  xlim = c(0,10),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #palette = colorRampPalette(mypal)(14),
  palette = c("#4DBBD5FF", "#E64B35FF"),
  #palette = "npg",
  #ggtheme = theme_bw(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

pcawg$ENSG00000264026 = factor(pcawg$ENSG00000264026, levels=c(0,1))
fit <- survfit(Surv(OS.time, OS) ~ ENSG00000264026, data = pcawg)
mo = coxph(Surv(OS.time, OS) ~ ENSG00000264026, data = pcawg)
pcawg_hr=round(summary(mo)$coefficients[2], digits=2)
pcawg_pval = round(summary(mo)$coefficients[5], digits=2)
v1 <- ggsurvplot(
  title=paste("LINC02003 (PCAWG)", "HR=", pcawg_hr, "Wald p=", pcawg_pval),
  fit,
  xlab = "Time (Years)",
  #  surv.median.line = "hv",
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  font.legend = 10,
  risk.table.fontsize = 5,
  data = pcawg,      # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  legend = "right",
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for
  # point estimaes of survival curves.
  xlim = c(0,10),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 1,     # break X axis in time intervals by 500.
  #palette = colorRampPalette(mypal)(14),
  palette = c("#4DBBD5FF", "#E64B35FF"),
  #palette = "npg",
  #ggtheme = theme_bw(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)
pdf("LIHC_TCGA_PCAWG_LINC02003.pdf")
print(t1)
print(v1)
dev.off()

#combine data to save as supplemetnary tables
#get original expression of genes
colnames(tcga)[1:2] = paste(colnames(tcga)[1:2], "tag", sep="_")
colnames(pcawg)[1:2] = paste(colnames(pcawg)[1:2], "tag", sep="_")
tcga$patient = rownames(tcga)
pcawg$patient = rownames(pcawg)

tcga_full = tcga_full[,c("patient", "ENSG00000178977", "ENSG00000264026")]
tcga_full$ENSG00000178977_median_cutoff = median(tcga_full$ENSG00000178977)
tcga_full$ENSG00000264026_median_cutoff = median(tcga_full$ENSG00000264026)

pcawg_full = pcawg_full[,c("patient", "ENSG00000178977", "ENSG00000264026")]
pcawg_full$ENSG00000178977_median_cutoff = median(pcawg_full$ENSG00000178977)
pcawg_full$ENSG00000264026_median_cutoff = median(pcawg_full$ENSG00000264026)

tcga = merge(tcga, tcga_full, by="patient")
pcawg = merge(pcawg, pcawg_full, by="patient")

tcga$data = "TCGA"
pcawg$data = "PCAWG"

all_data = rbind(tcga, pcawg)
write.csv(all_data, file="/u/kisaev/LIHC_case_study_dataset.csv", quote=F, row.names=F)
