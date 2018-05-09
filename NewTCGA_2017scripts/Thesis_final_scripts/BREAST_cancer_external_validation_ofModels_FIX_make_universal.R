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

#------BREAST CANCER CLINICAL---------------------------------------

breast_cancer = fread("nationwidechildrens.org_clinical_patient_brca.txt")
cols = c("bcr_patient_barcode", "tumor_status", "vital_status", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc", 
  "histological_type", "clinical_stage")
breast_cancer = as.data.frame(breast_cancer)
breast_cancer = breast_cancer[,which(colnames(breast_cancer) %in% cols)]
breast_cancer = breast_cancer[-c(1:2),]

#-------------------------------------------------------------------


#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% all_cands$Cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

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

#Breast cancer 
dtt = filtered_data_tagged[[13]]

#breast cancer protein coding gene data
gene = "ENSG00000237870"

#PCG data
z = which(pcg$patient %in% rownames(dtt))
breastpcg = pcg[z,]

#get CAV2 exp
z = which(colnames(pcg) %in% c("patient", "ENSG00000105971"))
pcg = pcg[,z]

#merge together
dtt$patient = rownames(dtt)
dtt =merge(dtt, pcg, by = "patient")

#look at how it's associated with survival 
z = which(rna$patient %in% dtt$patient)
breastrna = rna[z,]

z = which(colnames(breastrna) %in% c("patient", "ENSG00000237870"))
breastrna = as.data.frame(breastrna)
breastrna = breastrna[,z]

genexp = merge(breastrna, pcg, by="patient")
genexp[,2:3] = log1p(genexp[,2:3])
colnames(genexp)[2:3] = c("AC073130", "CAV2")

dtt = dtt[,which(colnames(dtt) %in% c("patient", "ENSG00000237870", "ENSG00000105971", "OS", "OS.time", "Cancer"))]
genexp = merge(dtt, genexp, by="patient")
genexp$OS = as.numeric(genexp$OS)
genexp$OS.time = as.numeric(genexp$OS.time)
cav2_med = median(as.numeric(genexp$CAV2))

for(y in 1:nrow(genexp)){
  exp = genexp$CAV2[y]
  if(exp >=cav2_med){
    genexp$ENSG00000105971[y] = 1
  }
  if(exp < cav2_med){
    genexp$ENSG00000105971[y] = 0
  }
}

genexp$OS.time = genexp$OS.time/365

coxph(Surv(OS.time, OS)  ~ ENSG00000237870 + ENSG00000105971, data = genexp)
pdf("CAV2_survival_plot_breastcancer_May8th.pdf", width=9)
fit <- survfit(Surv(OS.time, OS) ~ ENSG00000105971, data = genexp)
          s <- ggsurvplot(
          title = "CAV2 expression in Breast Cancer",
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = genexp,      # data used to fit survival curves. 
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
          #ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s) 
dev.off()      



# Basic plot
# +++++++++++++++++++++++++++
pdf("breast_cancer_CAV2_lncRNA_correlation_May8.pdf")
ggscatter(genexp, x = "AC073130", y = "CAV2",
   color = "ENSG00000237870", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   #conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   )
dev.off()



pdf("breast_cancer_candidates_survival_plots_final_cands_May7th.pdf")

  results_cox1 <- as.data.frame(matrix(ncol=7)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer")

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
  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$Cancer[1])
    
  names(row) <- names(results_cox1) 
  results_cox1 = rbind(results_cox1, row)
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$OS.time = newdat$OS.time/365
  #add what mutation status they are 
  newdat$er_mut = ""
  newdat$pr_mut = ""
  newdat$her2_mut = ""

  for(l in 1:nrow(newdat)){
    z = which(breast_cancer$bcr_patient_barcode == rownames(newdat)[l])
    if(!(length(z)) == 0){
      newdat$er_mut[l] = breast_cancer$er_status_by_ihc[z]
      newdat$pr_mut[l] = breast_cancer$pr_status_by_ihc[z]
      newdat$her2_mut[l] = breast_cancer$her2_status_by_ihc[z]
    }
  }  

  newdat$er_mut[newdat$er_mut == "Negative"] =0
  newdat$er_mut[newdat$er_mut == "Positive"] = 1
  newdat$pr_mut[newdat$pr_mut == "Negative"] = 0
  newdat$pr_mut[newdat$pr_mut == "Positive"]= 1
  newdat$her2_mut[newdat$her2_mut == "Negative"] = 0
  newdat$her2_mut[newdat$her2_mut == "Positive"] = 1

  newdat = subset(newdat, er_mut %in% c(0,1))
  newdat = subset(newdat, pr_mut %in% c(0,1))
  newdat = subset(newdat, her2_mut %in% c(0,1))

  all_muts = as.data.table(table(newdat$gene, newdat$er_mut, newdat$pr_mut, newdat$her2_mut))
  all_muts = filter(all_muts, N >0)
  colnames(all_muts) = c("lncRNAhighlow", "er_mut", "pr_mut", "her2_mut", "count")
  all_muts = as.data.table(all_muts)
  all_muts = all_muts[order(lncRNAhighlow)]

  #combine mutation presence into one variable
  newdat$combo = ""
  z = which((newdat$her2_mut) == 0 & (newdat$pr_mut == 0))
  newdat$combo[z] = "HER2_PR_neg"
  z = which((newdat$her2_mut) == 0 & (newdat$er_mut == 0))
  newdat$combo[z] = "HER2_ER_neg"
  z = which((newdat$er_mut) == 0 & (newdat$pr_mut == 0))
  newdat$combo[z] = "ER_PR_neg"
  z = which((newdat$er_mut) == 0 & (newdat$pr_mut == 0) & (newdat$her2_mut == 0))
  newdat$combo[z] = "TripleNegative"
  
  z = which((newdat$her2_mut) == 1 & (newdat$pr_mut == 1))
  newdat$combo[z] = "HER2_PR_pos"
  z = which((newdat$her2_mut) == 1 & (newdat$er_mut == 1))
  newdat$combo[z] = "HER2_ER_pos"
  z = which((newdat$er_mut) == 1 & (newdat$pr_mut == 1))
  newdat$combo[z] = "ER_PR_pos"
  z = which((newdat$er_mut) == 1 & (newdat$pr_mut == 1) & (newdat$her2_mut == 1))
  newdat$combo[z] = "TriplePos"

  lncs = coxph(Surv(OS.time, OS)  ~ gene + combo, data = newdat)


  fit <- survfit(Surv(OS.time, OS) ~ gene + er_mut + pr_mut + her2_mut, data = newdat)
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
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
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
          #palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

          s = s$plot + theme_bw() + theme (legend.position = "right")+
                  facet_grid(er_mut ~ gene)
          print(s) 
        
}

results_cox1 = results_cox1[-1,]

dev.off()







































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
  z = which(str_detect(colnames(dat), "ENSG"))

  for(i in 1:length(z)){

  gene = z[i]  
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

saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")














