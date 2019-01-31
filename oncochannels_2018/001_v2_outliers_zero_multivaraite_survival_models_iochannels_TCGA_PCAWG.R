
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#author: Karina Isaev, karin.isaev@gmail.com 
#date updated: Sept 21, 2018
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#-------------------------------------------------------------------
#this script uses data from TCGA (gene expression and clinical) to evaluate 
#the prognositc value of ion channels across different cancer types 
#working directory (source files, data files)
#/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ

#------Load libraries and scripts-----------------------------------

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#this script below prepares the RNA and clinical files for analysis 
source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

source("check_lnc_exp_cancers.R")

#------DATA---------------------------------------------------------

#get full dataset of GBM patients 
ext = readRDS("all_genes_external_tcga_all_cancers_March13_wclinical_data.rds")

#check if cands are significant using data from ext 
pats = as.data.table(table(ext$type))
pats = as.data.table(filter(pats, N >= 15))
colnames(pats)[1] = "type"
pats = merge(pats, canc_conv, by="type")

#get gbm
gbm = subset(ext, type=="GBM")

z = which(colnames(all) %in% colnames(gbm))
cols = colnames(all)[z]
all = all[,z]

z = which(colnames(gbm) %in% colnames(all))
cols = colnames(gbm)[z]
gbm = gbm[,z]

r = rbind(all, gbm)

all = r

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
colnames(ucsc)[8] = "HGNC.symbol"
cands = merge(cands, ucsc, by = "HGNC.symbol")
z = which(duplicated(all$patient))
if(!(length(z)==0)){
  dups = all$patient[z]
  k = which(all$patient %in% dups)
  all = all[-k,]
}

#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set

#--------------------------------------------------------------------

#1. Get cancer data (gene expression and clinical for each cancer type)
cancers = unique(all$type)
get_canc_data = function(canc){
  sub = subset(all, type == canc)
  return(sub)
} 
cancer_data = llply(cancers, get_canc_data)

#2. Subset data to ion channels and clinical variables 
get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(cands$hg19.ensGene.name2, "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "type"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

#3. Add tags to each ion channel to indicate High or Low expression 
#if median zero --> zero based split
#if not, outlier based split 

add_tags = function(dtt){
  print(dtt$type[1])
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
  
  rm = c()

  for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dtt[,names(medians)[k]] > 0)
    l2 = which(dtt[,names(medians)[k]] ==0)
    dtt[l1,names(medians)[k]] = 1
    dtt[l2, names(medians)[k]] = 0
    }

    if(!(med ==0)){
    #use outlier to split them
    iqr = 1.5 * (summary(dtt[,which(colnames(dtt) == names(medians[k]))])[5] - summary(dtt[,which(colnames(dtt) == names(medians[k]))])[2])
    q3 = summary(dtt[,which(colnames(dtt) == names(medians[k]))])[5] + iqr

    new_med = q3
    l1 = which(dtt[,names(medians)[k]] >= new_med)
    l2 = which(dtt[,names(medians)[k]] < new_med)
    
    if((length(l1) >= 10) & (length(l2) >=10)){
      dtt[l1,names(medians)[k]] = 1
      dtt[l2, names(medians)[k]] = 0
    }
  
   if(!((length(l1) >= 10) & (length(l2) >=10))){
      rm = c(rm, names(medians)[k])
    }

  }
  }

  z = which(colnames(dtt) %in% rm)
  if(!(length(z)==0)){
  dtt = dtt[,-z]}

  return(dtt)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")
saveRDS(filtered_data_tagged, file="28_cancers_tagged_by_ion_channels_zeromed_outliers.rds")

#3. Fit survival models and plot Kaplan Meier plot 

get_survival_models = function(dtt){
  results_cox1 <- as.data.frame(matrix(ncol=11)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "ic_test_ph", 'global_test_ph', "num_risk", "perc_risk")

  dat = dtt
  dat$Cancer = NULL
  dat$new_tumor_event_type = NULL
  dat$treatment_outcome_first_course = NULL
  dat$type = NULL

  #only evaluate genes that have patients with both high and low
  #expression and that clinical variables have a balanced split of
  #variables

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
  #dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)
  num_genes = which(str_detect(colnames(dat), "ENSG"))

  #save KM plots for each lncRNA for each cancer type sepereatley 
  file = paste("ION_CHANNELS_TCGA_SEPT2018/", dtt$type[1], ".pdf", sep="_")
  pdf(file)
  
  for(i in 1:length(num_genes)){
  gene = num_genes[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]
  c1 = table(newdat[,1])[1] > 10
  c2 = table(newdat[,1])[2] > 10

  if(c1 & c2){
  ionchannel_model = coxph(Surv(OS.time, OS)  ~ ., data = newdat)
  test.ph <- cox.zph(ionchannel_model)
  ic_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  hr = summary(ionchannel_model)$coefficients[1,c(1,2,5)][2]
  print(hr)
  if(hr >1){
    risk_num = length(which(newdat[,1] == 1))
    perc = risk_num/nrow(newdat)
  }

  if(hr < 1){
    risk_num = length(which(newdat[,1] == 0))
    perc = risk_num/nrow(newdat)
  }

  row <- c(colnames(newdat)[1], summary(ionchannel_model)$coefficients[1,c(1,2,5)],  summary(ionchannel_model)$conf.int[1,c(3,4)], dtt$type[1], 
    ic_test_ph, global, risk_num, perc)

  names(row) <- names(results_cox1) 
  results_cox1 = rbind(results_cox1, row)
  gene = colnames(newdat)[1]
  name_ic = cands$HGNC.symbol[cands$hg19.ensGene.name2==gene]
  colnames(newdat)[1] = "gene"
  newdat$OS.time = newdat$OS.time/365
  print(gene)
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, dtt$type[1], name_ic, "\nHR =", round(as.numeric(row[3]), digits=4)),
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
          xlim = c(0,6),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          #print(s)

   #generate boxplot 
   z = which(cancers == dtt$type[1])
   exp_data = cancer_data[[z]]
   exp_data = exp_data[,which(colnames(exp_data) %in% c(gene, "patient"))]
   newdat$patient = rownames(newdat)
   exp_data = merge(exp_data, newdat, by="patient")
   colnames(exp_data)[2] = "geneexp"
   exp_data$geneexp = log1p(exp_data$geneexp)
   p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", dtt$Cancer[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw()) + stat_n_text()
        # Change method
  p = p + stat_compare_means(method = "wilcox.test") 
  print(p)
  print(dtt$type[1])
}
}
dev.off()

results_cox1 = results_cox1[-1,]
#fdr on p-values 
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")

return(results_cox1)

}

library(parallel)

#tcga_results = mclapply(filtered_data_tagged, get_survival_models, mc.cores = 4) 
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$ic_test_ph = as.numeric(tcga_results1$ic_test_ph)
tcga_results1$global_test_ph = as.numeric(tcga_results1$global_test_ph)
tcga_results1$fdr_pval = as.numeric(tcga_results1$fdr_pval)
tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval)]

saveRDS(tcga_results1, file="TCGA_ION_CHANNEL_results_Sept21.rds")

#check which models violate the PH assumption
#to those models add age * survival time interaction 
which(tcga_results1$global_test_ph <= 0.05)
tcga_results1$num_risk = as.numeric(tcga_results1$num_risk)
tcga_results1[tcga_results1$num_risk <15,]

#plot sig KM plots
sig_fdr = as.data.table(filter(tcga_results1, fdr_pval <= 0.05))
source("check_lnc_exp_cancers.R")

#convert back to cancer codes
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "cancer"
sig_fdr = merge(sig_fdr, canc_conv, by="cancer")
sig_fdr = sig_fdr[order(fdr_pval)]
write.csv(sig_fdr, file="231_sig_IonChannels_survival_fdr_corrected_Sept21.csv", quote=F, row.names=F)

genes = sig_fdr$gene
cancs = sig_fdr$type

get_km_plot = function(gene, cancer){
  all_g = all
  all_g = as.data.frame(all_g)
  dat = all[,c(which(colnames(all) %in% c("type", gene, "OS", "OS.time")))]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients 
  med = median(dat$gene)
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)
  dat$OS.time = dat$OS.time/365
  dat$gene = factor(dat$gene, levels = c(0,1))
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(get_name_pcg(gene), dat$type[1]),
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
          palette = c("blue", "red"),
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}
} 

pdf("fdr_sig_ion_channels_km_plots_sept21.pdf")
mapply(get_km_plot, genes, cancs)
dev.off()

#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

tcga_results1 = readRDS("TCGA_results_multivariate_results_June22.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")

head(cands)
cands$combo = paste(cands$Geneid, cands$Cancer, sep="_")

pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_May2.rds")
pcawg_data$combo = paste(pcawg_data$canc, pcawg_data$histo)

pcawg_data$canc[pcawg_data$combo == "Uterus Adenocarcinoma, endometrioid"] = "Uterine Corpus Endometrial Carcinoma" 
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, clear cell type"] = "Kidney renal clear cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Breast Infiltrating duct carcinoma"] = "Breast invasive carcinoma"
pcawg_data$canc[pcawg_data$combo == "Ovary Serous cystadenocarcinoma"] = "Ovarian serous cystadenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "Pancreas Pancreatic ductal carcinoma"] = "Pancreatic adenocarcinoma" 
pcawg_data$canc[pcawg_data$combo == "Liver Hepatocellular carcinoma"] = "Liver hepatocellular carcinoma"
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, papillary type"] = "Kidney renal papillary cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Squamous cell carcinoma"] = "Lung squamous cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Adenocarcinoma, invasive"] = "Lung adenocarcinoma"

cancers_tests = as.list(unique(tcga_results1$cancer[which(tcga_results1$cancer %in% pcawg_data$canc)]))

get_matched_data = function(cancer){
    dtt = subset(pcawg_data, canc == cancer)
    z = which(colnames(dtt) %in% c(as.character(cands$Geneid[cands$Cancer == dtt$canc[1]]), "canc", 
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
  results_cox1 <- as.data.frame(matrix(ncol=11)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "num_risk", "perc_risk")

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
  test.ph <- cox.zph(lncs)
  lnc_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  hr = summary(lncs)$coefficients[1,c(1,2,5)][2]
  if(hr >1){
    risk_num = length(which(newdat[,1] == 1))
    perc = risk_num/nrow(newdat)
  }

  if(hr < 1){
    risk_num = length(which(newdat[,1] == 0))
    perc = risk_num/nrow(newdat)
  }

  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2)],  glance(lncs)[4],  summary(lncs)$conf.int[1,c(3,4)], dtt$canc[1], 
    lnc_test_ph, global, risk_num, perc)
    
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
          #xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
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
      exp_data$geneexp = log1p(exp_data$geneexp)
      p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", dtt$canc[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw())
        # Change method
       p = p + stat_compare_means(method = "wilcox.test")
       print(p)
}

results_cox1 = results_cox1[-1,]
return(results_cox1)
}

pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_final_cands_May3rd_full_lifespan.pdf")
pcawg_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
pcawg_results1 = ldply(pcawg_results, data.frame)
pcawg_results1$pval = as.numeric(pcawg_results1$pval)
pcawg_results1$fdr_pval = p.adjust(pcawg_results1$pval, method="fdr")

#combine results from TCGA and PCAWG
pcawg_results1$data = "PCAWG"
pcawg_results1$num_risk = as.numeric(pcawg_results1$num_risk)
pcawg_results1 = filter(pcawg_results1, num_risk >=10)

z = which(pcawg_results1$upper95 == "Inf")
pcawg_results1 = pcawg_results1[-z,]

pcawg_results1$combo = paste(pcawg_results1$gene, pcawg_results1$cancer, sep="_")

#all-results
#tcga_results1$lnc_test_ph =NULL
#tcga_results1$global_test_ph = NULL
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

colnames(lnc_info)[1] = "gene"

all_results = merge(all_results, lnc_info, by="gene")
all_results = all_results[order(data, pval)]
all_results$best_cell_trait_fdr = NULL
all_results$CAT_browser_link = NULL

#this file was created without first filtering to only include the robust lncRNAs
#so contains results for all 173 lncRNAs 
saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#this file was created by FIRST filtering to inlclude ONLY robust lncRNAs
#saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_Aug8.rds")

#-----check which actually match---------------------------------------------------------------------------

all_results_orig = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#all_results_orig = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_Aug8.rds")

all_results = all_results_orig[!duplicated(all_results_orig), ]
all_results = filter(all_results, pval <= 0.05)

lncs = as.list(as.character(unique(all_results$gene[all_results$data == "PCAWG"])))

check_match = function(lnc){
  z = which(all_results$gene == lnc)
  res = as.data.table(all_results[z,])
  
  if(dim(res)[1] > 2){
  
  canc = res$cancer[which(duplicated(res$cancer))]
  res = filter(res, cancer == canc)
  }

  test = which(as.numeric(res$HR) >= 1)
  test = length(test)
  if(test ==2){
    match = "match"
  }
  if(test ==0){
    match = "match"
  }
  if(test ==1){
    match = "nomatch"
  }
  canc = unique(res$canc)
  return(c(lnc, match, canc))
}

matches = llply(lncs, check_match)
matches = do.call(rbind.data.frame, matches)
matches = as.data.table(matches)
colnames(matches) = c("lnc", "match", "cancer")
matches = filter(matches, match == "match")
colnames(matches)[1] = "gene"
matches = merge(matches, all_results, by=c("gene", "cancer"))
matches$combo = paste(matches$gene, matches$cancer, sep="_")

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]
colnames(ucsc)[6] = "gene"
matches = merge(matches, ucsc, by=c("gene"))

#remove weird KIRC one (already did above by filtering)
#z = which(matches$gene == "ENSG00000250360")
#matches = matches[-z,]

write.table(matches, file="6_unique_lncNRAs_validate_PCAWG.txt", quote=F, row.names=F, sep=";")

write.table(matches, file="4_unique_lncNRAs_validate_PCAWG.txt", quote=F, row.names=F, sep=";")
matches = as.data.frame(matches)
matches = matches[,1:13]
pdf("6_unique_lncNRAs_validate_PCAWG.pdf", width=24)
p<-tableGrob(matches)
grid.arrange(p)
dev.off()
matches$combo = paste(matches$gene, matches$cancer, sep="_")

#Clean up and save as spreadsheet
all_results_orig$combo = paste(all_results_orig$gene, all_results_orig$cancer, sep="_")
z = which(all_results_orig$combo %in% matches$combo)
all_results_orig$top_pcawg_val = ""
all_results_orig$top_pcawg_val[z] = "YES"
all_results_orig$combo = NULL
all_results_orig$best_eQTL_mRNA_fdr = NULL

all_results_orig$coef = as.numeric(all_results_orig$coef)
all_results_orig$coef = round(all_results_orig$coef, digits=4)

all_results_orig$HR = as.numeric(all_results_orig$HR)
all_results_orig$HR = round(all_results_orig$HR, digits=4)

all_results_orig$pval = as.numeric(all_results_orig$pval)
all_results_orig$pval = round(all_results_orig$pval, digits=6)

all_results_orig$low95 = as.numeric(all_results_orig$low95)
all_results_orig$low95 = round(all_results_orig$low95, digits=4)

all_results_orig$upper95 = as.numeric(all_results_orig$upper95)
all_results_orig$upper95 = round(all_results_orig$upper95, digits=4)

all_results_orig$lnc_test_ph = as.numeric(all_results_orig$lnc_test_ph)
all_results_orig$lnc_test_ph = round(all_results_orig$lnc_test_ph, digits=4)

all_results_orig$global_test_ph = as.numeric(all_results_orig$global_test_ph)
all_results_orig$global_test_ph = round(all_results_orig$global_test_ph, digits=4)

all_results_orig$perc_risk = as.numeric(all_results_orig$perc_risk)
all_results_orig$perc_risk = round(all_results_orig$perc_risk, digits=4)

all_results_orig$fdr_pval = as.numeric(all_results_orig$fdr_pval)
all_results_orig$fdr_pval = round(all_results_orig$fdr_pval, digits=4)

write.csv(all_results_orig, file="175_lncRNA_cancers_combos_23_cancer_types_july5.csv", quote=F, row.names=F)
#write.csv(all_results_orig, file="112_lncRNA_cancers_combos_22_cancer_types_aug8.csv", quote=F, row.names=F)





