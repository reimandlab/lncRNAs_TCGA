library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)

#------FEATURES-----------------------------------------------------

all_cands = readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")

#-------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% all_cands$Cancer)
cancer_data = canc_datas[z]

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(all_cands$Name[all_cands$Cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

add_tags = function(dtt){
  rownames(dtt) = dtt$patient
  dtt$patient = NULL

  #log1p 
  z = which(str_detect(colnames(dtt), "ENSG"))
  if(length(z)>1){
  medians = apply(dtt[,z], 2, median)}
  if(length(z)==1){
    medians = median(dtt[,z])
  }
  #add high low tag
  for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dtt[,k] > 0)
    l2 = which(dtt[,k] ==0)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(dtt[,k] >= med)
    l2 = which(dtt[,k] < med)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }
  }  
  return(dtt)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")

get_survival_models = function(dtt){
  results_cox1 <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95")

  dat = dtt
  dat$Cancer = NULL
  dat$new_tumor_event_type = NULL
  dat$treatment_outcome_first_course = NULL

  z = which(str_detect(colnames(dtt), "ENSG"))
  check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
  }
  keep = unlist(apply(dat, 2, check_contrasts))
  
  dat = dat[,c(which(colnames(dat) %in% names(keep)))]

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)
  dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)

  for(i in 1:length(z)){

  gene = z[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]
  
  lncs = coxph(Surv(OS.time, OS)  ~ ., data = newdat)
  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)])
    
  names(row) <- names(results_cox1) 
  results_cox1 = rbind(results_cox1, row)
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$OS.time = newdat$OS.time/365
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, dtt$Cancer[1], "HR =", round(as.numeric(row[3]), digits=4)),
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = newdat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}

results_cox1 = results_cox1[-1,]
return(results_cox1)

}

pdf("TCGA_candidates_survival_plots_final_cands_May3rd.pdf")
llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()











#------PCAWG DATA---------------------------------------------------
pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March20.rds")
pcawg_data = subset(pcawg_data, canc == "Ovary Serous cystadenocarcinoma")
z = which(colnames(pcawg_data) %in% c(colnames(canc_data), "time", "status"))
pcawg_data = pcawg_data[,z]
#add high low tag
medians = apply(pcawg_data[,1:(ncol(pcawg_data)-2)], 2, median)
for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(pcawg_data[,k] > 0)
    l2 = which(pcawg_data[,k] ==0)
    pcawg_data[l1,k] = 1
    pcawg_data[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(pcawg_data[,k] >= med)
    l2 = which(pcawg_data[,k] < med)
    pcawg_data[l1,k] = 1
    pcawg_data[l2, k] = 0
    }
}  

pcawg_data$time = as.numeric(pcawg_data$time)
pcawg_data$status = as.numeric(pcawg_data$status)

#------individual lncs-----------------------------------------

pdf("topCands_OV_individual_survivaplots_TCGA.pdf")
results_cox1 <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
for(i in 1:(ncol(canc_data)-2)){
	dat = canc_data[,c(i, ncol(canc_data), (ncol(canc_data)-1))]
	justlncs_pcawg = coxph(Surv(time, status)  ~ ., data = dat)
	row <- c(colnames(canc_data)[i], summary(justlncs_pcawg)$coefficients[1,c(1,2,5)],  summary(justlncs_pcawg)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox1)	
  	results_cox1 = rbind(results_cox1, row)
  	gene = colnames(canc_data)[i]
  	colnames(dat)[1] = "gene"
  	dat$time = dat$time/365
  	fit <- survfit(Surv(time, status) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene, "HR =", round(as.numeric(row[3]), digits=4)),
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}
dev.off()
results_cox1 = results_cox1[-1,]

#------individual lncs-----------------------------------------

pdf("topCands_OV_individual_survivaplots_PCAWG.pdf")
results_cox2 <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox2) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
for(i in 1:(ncol(pcawg_data)-2)){
	dat = pcawg_data[,c(i, ncol(pcawg_data), (ncol(pcawg_data)-1))]
	justlncs_pcawg = coxph(Surv(time, status)  ~ ., data = dat)
	row <- c(colnames(pcawg_data)[i], summary(justlncs_pcawg)$coefficients[1,c(1,2,5)],  summary(justlncs_pcawg)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox2)	
  	results_cox2 = rbind(results_cox2, row)
  	gene = colnames(pcawg_data)[i]
  	colnames(dat)[1] = "gene"
  	dat$time = dat$time/365
  	fit <- survfit(Surv(time, status) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene, "HR =", round(as.numeric(row[3]), digits=4)),
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}
dev.off()
results_cox2 = results_cox2[-1,]












####Test on the PCAWG data
canc_data$time = as.numeric(canc_data$time)
canc_data$status[canc_data$status=="Alive"] <- 0
canc_data$status[canc_data$status=="Dead"] <- 1
canc_data$status = as.numeric(canc_data$status)

#####Train model using all TCGA data and the chosen predictor lncRNAs 
canc_data = canc_data[,which(colnames(canc_data) %in% c(ov_features$V1, "time", "status"))]
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)
#keep = names(which(summary(justlncs)$coefficients[,5] <=0.05))
#if(!(length(keep)==0)){
#canc_data = canc_data[,c(which(colnames(canc_data) %in% c(keep,"time", "status")))]
#}
#updated model
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)


#using all 8 candidate lncRNAs 
pdf("timedependentAUC_externalPCAWG_OV.pdf")
lpnew <- predict(justlncs, newdata=pcawg_data)
Surv.rsp <- Surv(canc_data$time, canc_data$status)
Surv.rsp.new <- Surv(pcawg_data$time, pcawg_data$status)
times <- seq(10, 1000, 10)
AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
names(AUC_Uno)
AUC_Uno$iauc
plot(AUC_Uno)
abline(h = 0.5)
BrierScore <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times,
type = "brier", int.type = "weighted")
plot(BrierScore)
abline(h = 0.25)
dev.off()


