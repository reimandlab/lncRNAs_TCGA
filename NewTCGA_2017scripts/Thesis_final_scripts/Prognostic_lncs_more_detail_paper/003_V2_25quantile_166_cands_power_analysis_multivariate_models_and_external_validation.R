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
require("powerSurvEpi")


#------FEATURES-----------------------------------------------------

#cands -- should be this file
cands = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")

#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% cands$Cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(cands$Geneid[cands$Cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
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
  medians = apply(dtt[,z], 2, function(x){quantile(x)[2]})}
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
  results_cox1 <- as.data.frame(matrix(ncol=18)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "num_risk", "perc_risk", "power", "median_nonzero", "sd_nonzero", "min_nonzero", "max_nonzero", "multi_model_concordance", "lnc_only_concordance")

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
  test.ph <- cox.zph(lncs)
  lnc_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  #mutlvariate concordance 
  cmulti = as.numeric(glance(lncs)[11])
  lnc_only = coxph(Surv(OS.time, OS)  ~ newdat[,1], data = newdat)
  clnconly = as.numeric(glance(lnc_only)[11])

  #calculate power
  #k is ratio of participants in group E (experimental group) compared to group C (controlgroup).
  kval = as.numeric(table(newdat[,1])[1]/table(newdat[,1])[2])
  #m = expected total number of events over both groups.
  mval = length(which(newdat$OS == 1))

  power = powerCT.default0(k=kval,m=mval, RR=2, alpha=0.05)

  hr = summary(lncs)$coefficients[1,c(1,2,5)][2]
  if(hr >1){
    risk_num = length(which(newdat[,1] == 1))
    perc = risk_num/nrow(newdat)
  }

  if(hr < 1){
    risk_num = length(which(newdat[,1] == 0))
    perc = risk_num/nrow(newdat)
  }

  gene_name = colnames(newdat)[1]
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$OS.time = newdat$OS.time/365
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, dtt$Cancer[1], "HR =", hr, digits=4),
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
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
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

   #get median non-zero expression and SD 
   if(perc < 0.45){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
   } else if (perc > 0.55){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
   } else {
      median_nonzero = "notavail"
      sd_nonzero = "notavail"
      min_nonzero = "notavail"
      max_nonzero = "notavail"
   }

   row <- c(gene_name, summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$Cancer[1], 
    lnc_test_ph, global, risk_num, perc, power, median_nonzero,
    sd_nonzero,
    min_nonzero,
    max_nonzero, cmulti, clnconly)

   names(row) <- names(results_cox1) 
   results_cox1 = rbind(results_cox1, row)

   #make density plot using FPKM-UQ values and logged values
   #visualize bimodality 
   exp_data$median=""
   exp_data$median[exp_data$gene==0] = "Low"
   exp_data$median[exp_data$gene==1] = "High"
   exp_data$median = factor(exp_data$median, levels=c("High", "Low"))
   
   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" "))
   print(gg)

   exp_data$geneexp = log1p(exp_data$geneexp)
  
   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" "))
   print(gg)

   p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", dtt$Cancer[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test")
  print(p)
}

results_cox1 = results_cox1[-1,]
#fdr on p-values 
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")

return(results_cox1)

}

#pdf("TCGA_candidates_survival_plots_final_cands_FULL_lifespan_May3rd.pdf")
pdf("TCGA_25_percentile_candidates_survival_plots_final_cands_FULL_5year_surv_oct3.pdf")
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$lnc_test_ph = as.numeric(tcga_results1$lnc_test_ph)
tcga_results1$global_test_ph = as.numeric(tcga_results1$global_test_ph)
tcga_results1$fdr_pval = as.numeric(tcga_results1$fdr_pval)
tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval)]

saveRDS(tcga_results1, file="TCGA_25_percentile_results_multivariate_results_Oct3.rds")

#check which models violate the PH assumption
#to those models add age * survival time interaction 
which(tcga_results1$global_test_ph <= 0.05)
tcga_results1$num_risk = as.numeric(tcga_results1$num_risk)
tcga_results1[tcga_results1$num_risk <15,]

#plot distribution, cut number of risk patients 
#tcga_results1$groupy = cut(tcga_results1$num_risk, breaks =c(1, 20, 40, 60, 80, 100, 200, 300, 
#  400, 500, 600, 700, 800, 900, 1000))

#SUMMARIZE
# Change line color and fill color
#riskplot = gghistogram(tcga_results1, x = "groupy", y = "..count..", stat="count", 
#   fill = "white", color="black",  palette = c("#00AFBB", "#E7B800")) + xlab("Number of Risk patients")+
#ylab("Frequency")
#pdf("Dist_num_risk_patients_per_lncRNA.pdf", width=10)
#riskplot
#dev.off()

tcga_results1$perc_risk = as.numeric(tcga_results1$perc_risk)
riskplot = gghistogram(tcga_results1, x = "perc_risk", 
   fill = "white", color="black",  palette = c("#00AFBB", "#E7B800")) + xlab("Percentage of Risk patients")+
ylab("Frequency")
pdf("Dist_perc_risk_patients_per_lncRNA.pdf", width=10)
riskplot
dev.off()


#tcga_results1 = filter(tcga_results1, fdr_pval <=0.05)

#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

tcga_results1 = readRDS("TCGA_25_percentile_results_multivariate_results_Oct3.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")
#z = which(tcga_results1$combo %in% robust$combo)
#tcga_results1 = tcga_results1[z,]

#only the robust ones
#robust = readRDS(file="112_combos_robust_internal_validation_survival_lncRNAs_aug8.rds")
#robust = readRDS("148_combos_robust_5perc_increase_internal_validation_survival_lncRNAs_aug9.rds")
head(cands)
cands$combo = paste(cands$Geneid, cands$Cancer, sep="_")
#z = which(cands$combo %in% robust$combo)
#cands = cands[z,]

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
  medians = apply(dtt[,z], 2, function(x){quantile(x)[2]})}

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
  results_cox1 <- as.data.frame(matrix(ncol=18)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "num_risk", "perc_risk", "power", "median_nonzero", "sd_nonzero", "min_nonzero", "max_nonzero",
    "multi_model_concordance", "lnc_only_concordance")

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


  #mutlvariate concordance 
  cmulti = "not_avail"
  clnconly = as.numeric(glance(lncs)[11])

  #calculate power
  #k is ratio of participants in group E (experimental group) compared to group C (controlgroup).
  kval = as.numeric(table(newdat[,1])[1]/table(newdat[,1])[2])
  #m = expected total number of events over both groups.
  mval = length(which(newdat$status == 1))

  power = powerCT.default0(k=kval,m=mval, RR=2, alpha=0.05)

  gene_name = colnames(newdat)[1]

  hr = summary(lncs)$coefficients[1,c(1,2,5)][2]
  if(hr >1){
    risk_num = length(which(newdat[,1] == 1))
    perc = risk_num/nrow(newdat)
  }

  if(hr < 1){
    risk_num = length(which(newdat[,1] == 0))
    perc = risk_num/nrow(newdat)
  }
    
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$time = newdat$time/365
  fit <- survfit(Surv(time, status) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, dtt$canc[1], "HR =", hr, digits=4),
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
      
      #get median non-zero expression and SD 
      if(perc < 0.45){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
      } else if (perc > 0.55){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
      } else {
      median_nonzero = "notavail"
      sd_nonzero = "notavail"
      min_nonzero = "notavail"
      max_nonzero = "notavail"
      }

    row <- c(gene_name, summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$canc[1], 
    lnc_test_ph, global, risk_num, perc, power, median_nonzero,
    sd_nonzero,
    min_nonzero,
    max_nonzero, cmulti, clnconly)

    names(row) <- names(results_cox1) 
    results_cox1 = rbind(results_cox1, row)


    #make density plot using FPKM-UQ values and logged values
   #visualize bimodality 
   exp_data$median=""
   exp_data$median[exp_data$gene==0] = "Low"
   exp_data$median[exp_data$gene==1] = "High"
   exp_data$median = factor(exp_data$median, levels=c("High", "Low"))
   
   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$canc[1] , sep=" "))
   print(gg)

   exp_data$geneexp = log1p(exp_data$geneexp)
  
   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$canc[1] , sep=" "))
   print(gg)      

      p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene_name, "Expression", dtt$canc[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw())
        # Change method
       p = p + stat_compare_means(method = "wilcox.test")
       print(p)
}

results_cox1 = results_cox1[-1,]
return(results_cox1)
}

#pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_final_cands_May3rd_full_lifespan.pdf")
pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_Oct3_five_year_survival.pdf")
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

saveRDS(pcawg_results1, file="PCAWG_25_percentile_external_validation_lncCands_oct3.rds")

#plot power versus Hazard Ratio 
pcawg_results1$pval = as.numeric(pcawg_results1$pval)
pcawg_results1$power = as.numeric(pcawg_results1$power)

ggscatter(pcawg_results1, x = "power", y = "pval",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )
dev.off()


#all-results
#tcga_results1$lnc_test_ph =NULL
#tcga_results1$global_test_ph = NULL
#tcga_results1$groupy = NULL
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
saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15_25_percentile.rds")

#this file was created by FIRST filtering to inlclude ONLY robust lncRNAs
#saveRDS(all_results, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_Aug8.rds")

#-----check which actually match---------------------------------------------------------------------------

all_results_orig = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15_25_percentile.rds")
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

write.table(matches, file="5_unique_lncNRAs_validate_PCAWG_from_25_percentile_splitting.txt", quote=F, row.names=F, sep=";")

#write.table(matches, file="4_unique_lncNRAs_validate_PCAWG.txt", quote=F, row.names=F, sep=";")
matches = as.data.frame(matches)
matches = matches[,1:13]
pdf("5_unique_lncNRAs_validate_PCAWG.pdf", width=24)
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

all_results_orig$combo = paste(all_results_orig$gene, all_results_orig$cancer, sep="_")

write.csv(all_results_orig, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15_25_percentile.csv", quote=F, row.names=F)
#write.csv(all_results_orig, file="112_lncRNA_cancers_combos_22_cancer_types_aug8.csv", quote=F, row.names=F)
saveRDS(all_results_orig, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15_25_percentile.rds")





