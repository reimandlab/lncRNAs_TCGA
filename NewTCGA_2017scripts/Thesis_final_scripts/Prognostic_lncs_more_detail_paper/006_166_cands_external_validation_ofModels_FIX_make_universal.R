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
library(caret)  
library(Rtsne)

#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands


#--------This script ------------------------------------------------

#build models using all of TCGA for each lncRNA candidate
#conduct bootstrapping on PCAWG sample 
#and evluate perforamnce of TCGA candidate model 

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% all_cands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(all_cands$gene[all_cands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

add_tags = function(dtt){
  print(dtt$Cancer[1])
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
  results_cox1 <- as.data.frame(matrix(ncol=7)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer")

  dat = dtt
  dat$Cancer = NULL
  dat$new_tumor_event_type = NULL
  dat$treatment_outcome_first_course = NULL

  num_genes = which(str_detect(colnames(dtt), "ENSG"))
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

  for(i in 1:length(num_genes)){
  gene = num_genes[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]
  
  lncs = coxph(Surv(OS.time, OS)  ~ ., data = newdat)
  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$Cancer[1])
    
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

   #generate boxplot 
   z = which(cancers == dtt$Cancer[1])
   exp_data = canc_datas[[z]]
   exp_data = exp_data[,which(colnames(exp_data) %in% c(gene, "patient"))]
   newdat$patient = rownames(newdat)
   exp_data = merge(exp_data, newdat, by="patient")
   colnames(exp_data)[2] = "geneexp"
   p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", dtt$Cancer[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_minimal())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test")
  print(p)
}

results_cox1 = results_cox1[-1,]
return(results_cox1)

}

pdf("TCGA_candidates_survival_plots_final_cands_May3rd.pdf")
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$pval = as.numeric(tcga_results1$pval)
tcga_results1 = filter(tcga_results1, pval <=0.05)

#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_May2.rds")
pcawg_data$canc[pcawg_data$histo == "Adenocarcinoma, endometrioid"] = "Uterine Corpus Endometrial Carcinoma"
pcawg_data$canc[pcawg_data$histo == "Adenocarcinoma, clear cell type"] = "Kidney renal clear cell carcinoma"
pcawg_data$canc[pcawg_data$histo == "Pancreatic ductal carcinoma"] = "Pancreatic adenocarcinoma"
pcawg_data$canc[pcawg_data$histo == "Hepatocellular carcinoma"] = "Liver hepatocellular carcinoma"
pcawg_data$canc[pcawg_data$histo == "Squamous cell carcinoma"] = "Lung squamous cell carcinoma"
pcawg_data$canc[pcawg_data$histo == "Adenocarcinoma, invasive"] = "Lung adenocarcinoma"
pcawg_data$canc[pcawg_data$histo == "Serous cystadenocarcinoma"] = "Ovarian serous cystadenocarcinoma"
pcawg_data$canc[pcawg_data$histo == "Infiltrating duct carcinoma"] = "Breast invasive carcinoma"

cancers_tests = as.list(unique(tcga_results1$cancer[which(tcga_results1$cancer %in% pcawg_data$canc)]))

get_matched_data = function(cancer){
    dtt = subset(pcawg_data, canc == cancer)
    z = which(colnames(dtt) %in% c(as.character(all_cands$Name[all_cands$Cancer == dtt$canc[1]]), "canc", 
    "histo", "time", "status", "sex"))
    dtt = dtt[,z]
    return(dtt)
}

filtered_data = llply(cancers_tests, get_matched_data)

add_tags = function(dtt){

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
  results_cox1 <- as.data.frame(matrix(ncol=7)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer")

  dat = dtt
  dat$canc = NULL
  dat$histo = NULL
  dat$sex = NULL
  dat$status[dat$status=="alive"] =0
  dat$status[dat$status=="deceased"] =1

  z = which(str_detect(colnames(dat), "ENSG"))
  check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
  }
  keep = unlist(apply(dat, 2, check_contrasts))
  
  dat = dat[,c(which(colnames(dat) %in% names(keep)))]

  dat$status = as.numeric(dat$status)
  dat$time = as.numeric(dat$time)
  num_genes = which(str_detect(colnames(dat), "ENSG"))

  for(i in 1:length(num_genes)){
  print(i)
  gene = num_genes[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]
  
  lncs = coxph(Surv(time, status)  ~ ., data = newdat)
  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$canc[1])
    
  names(row) <- names(results_cox1) 
  results_cox1 = rbind(results_cox1, row)
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$time = newdat$time/365
  fit <- survfit(Surv(time, status) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, dtt$canc[1], "HR =", round(as.numeric(row[3]), digits=4)),
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
    #generate boxplot 
      z = which(cancers_tests == dtt$canc[1])
      exp_data = filtered_data[[z]]
      exp_data$patient = rownames(exp_data)
      exp_data = exp_data[,which(colnames(exp_data) %in% c(gene, "patient"))]
      newdat$patient = rownames(newdat)
      exp_data = merge(exp_data, newdat, by="patient")
      colnames(exp_data)[2] = "geneexp"
      p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", dtt$canc[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_minimal())
        # Change method
       p = p + stat_compare_means(method = "wilcox.test")
       print(p)
}

results_cox1 = results_cox1[-1,]
return(results_cox1)
}

pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_final_cands_May3rd.pdf")
pcawg_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
pcawg_results1 = ldply(pcawg_results, data.frame)

#combine results from TCGA and PCAWG
tcga_results1$data = "TCGA"
pcawg_results1$data = "PCAWG"
all_results = as.data.table(rbind(tcga_results1, pcawg_results1))
all_results = all_results[order(gene, cancer)]

#add more info on lncRNAs
lnc_info = read.csv("fantom_genebased_evidence_supp_table_17.csv")
lnc_info = lnc_info[which(str_detect(lnc_info$CAT_geneID, "ENSG")),]

#shorten the gene names 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
lnc_info[,1] <- apply(lnc_info[,1:2], 1, extract3) #5049 lncRNAs 

#how many conserved 
#1207 have conserved exons 
z1 = which(!(lnc_info$exon_RS_score == "__na"))

#924 have conserved TIRs 
z2 = which(!(lnc_info$TIR_RS_score == "__na"))

conserved = lnc_info[unique(c(z1, z2)),]
colnames(lnc_info)[1] = "gene"


all_results = merge(all_results, lnc_info, by="gene")
colnames(all_cands)[1] = "gene"
all_results = merge(all_results, all_cands, by ="gene")

all_results = all_results[order(data, pval)]


saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")














