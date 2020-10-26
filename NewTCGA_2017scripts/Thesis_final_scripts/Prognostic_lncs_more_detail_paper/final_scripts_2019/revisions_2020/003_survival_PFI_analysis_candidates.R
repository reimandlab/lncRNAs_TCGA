set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FEATURES-----------------------------------------------------

#cands -- should be this file
#cands = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")

cands = readRDS("lncRNAs_selected_by_EN_april14.rds") #1000 runs of cross-validations using new updated dataset (GBM=124, OV and LUAD)

#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set
#fantom 

#--------------------------------------------------------------------

cancers = unique(cands$canc)
get_canc = function(canc){
  rna_dat = subset(rna, Cancer == canc)
  return(rna_dat)
}

cancer_data = llply(cancers, get_canc) #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(cands$gene[cands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "PFI.time", "PFI", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "type"))
  print(dtt$type[1])
  #print(cands$cancer[1])
  dtt = dtt[,..z]

  t=filter(as.data.table(table(dtt$histological_grade)), N <=5)$V1
  print(t)
  if((!(length(t)==0))){
    z=which(dtt$histological_grade %in% t)
    dtt = dtt[-z,]
  }

  t=filter(as.data.table(table(dtt$clinical_stage)), N <=5)$V1
  print(t)
  if((!(length(t)==0))){
    z=which(dtt$clinical_stage %in% t)
    dtt = dtt[-z,]
  }

  t=filter(as.data.table(table(dtt$race)), N <=5)$V1
  print(t)
  if((!(length(t)==0))){
    z=which(dtt$race %in% t)
    dtt = dtt[-z,]
  }

  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)
all_tcga_cancers = as.data.table(ldply(filtered_data))
all_tcga_cancers = unique(all_tcga_cancers[,c("patient", "type", "age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "new_tumor_event_type", "treatment_outcome_first_course", "PFI", "PFI.time", "Cancer")])
write.csv(all_tcga_cancers, file="TCGA_clinical_data_patients_used_external_survival_analysis_PFI.csv", quote=F, row.names=F)

add_tags = function(dtt){
  dtt = as.data.frame(dtt)
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
saveRDS(filtered_data_tagged, file="22_cancer_types_with_lncRNA_candidates_labelled_high_low.rds") #21 actually but just keep file names the same for now

get_survival_models = function(dtt){

  print(dtt$Cancer[1])

  results_cox1 <- as.data.frame(matrix(ncol=21)) ; colnames(results_cox1) <- c("gene", "coef", "pval", "HR", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "num_risk", "perc_risk", "median_nonzero", "sd_nonzero", "min_nonzero", "max_nonzero", "multi_model_concordance", 
    "lnc_only_concordance", "clinical_only_concordance", "num_events", "perc_wevents", "anova_pval")

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

  dat$PFI = as.numeric(dat$PFI)
  dat$PFI.time = as.numeric(dat$PFI.time)
  dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)
  #num events
  num_events = length(which(dat$PFI == 1))
  perc_events = num_events/nrow(dtt)

  for(i in 1:length(num_genes)){
  gene = num_genes[i]  
  print(gene)
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]

  lncs = coxph(Surv(PFI.time, PFI)  ~ ., data = newdat)
  test.ph <- cox.zph(lncs)
  lnc_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  #mutlvariate concordance 
  cmulti = as.numeric(glance(lncs)[13])
  lnc_only_model = coxph(Surv(PFI.time, PFI)  ~ newdat[,1], data = newdat)
  lnc_only = as.numeric(glance(lnc_only_model)[13])
  #hr = summary(lnc_only_model)$coefficients[2]

  z = which(str_detect(colnames(newdat), "ENSG"))
  clin_only = newdat[,-z]
  clinical_only_model = coxph(Surv(PFI.time, PFI)  ~ ., data = clin_only)
  clinical_only = as.numeric(glance(clinical_only_model)[13])
  lr = anova(clinical_only_model, lncs)
  lr_pval = lr[2,4]

  #calculate power
  #k is ratio of participants in group E (experimental group) compared to group C (controlgroup).
  kval = as.numeric(table(newdat[,1])[1]/table(newdat[,1])[2])
  #m = expected total number of events over both groups.
  mval = length(which(newdat$PFI == 1))

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
  newdat$PFI.time = newdat$PFI.time/365
  newdat$gene = factor(newdat$gene, levels=c(1,0))
  fit <- survfit(Surv(PFI.time, PFI) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(get_name(gene_name), canc_conv$type[canc_conv$Cancer == dtt$Cancer[1]][1], "HR =", round(hr, digits=2)),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = newdat,      # data used to fit survival curves. 
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
          print(s)

   #generate boxplot 
   z = which(rna$Cancer == dtt$Cancer[1])
   exp_data = rna[z,]
   exp_data = as.data.frame(exp_data)
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

   row <- c(gene_name, summary(lncs)$coefficients[1,c(1,5)], hr,  summary(lncs)$conf.int[1,c(3,4)], dtt$Cancer[1], 
    lnc_test_ph, global, risk_num, perc, median_nonzero,
    sd_nonzero,
    min_nonzero,
    max_nonzero, cmulti, lnc_only, clinical_only, num_events, perc_events, lr_pval)

   names(row) <- names(results_cox1) 
   results_cox1 = rbind(results_cox1, row)

   #make density plot using FPKM-UQ values and logged values
   #visualize bimodality 
   exp_data$median=""
   exp_data$median[exp_data$gene==0] = "Low"
   exp_data$median[exp_data$gene==1] = "High"
   exp_data$median = factor(exp_data$median, levels=c("High", "Low"))
   
   #gg <- ggplot(exp_data)
   #gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   #gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" "))
   #print(gg)

   exp_data$geneexp = log1p(exp_data$geneexp)
  
   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" ")) + labs(y="log1p(FPKM-UQ)")
   print(gg)

   p <- ggboxplot(exp_data, x = "median", y = "geneexp",
          color = "median",
         palette = mypal[c(4,1)], title = paste(get_name(gene_name), "Expression", canc_conv$type[canc_conv$Cancer == dtt$Cancer[1]][1], sep=" "), 
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_classic())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test") + stat_n_text() + scale_color_npg() 
  print(p)
}

results_cox1 = results_cox1[-1,]
#fdr on p-values 
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")

return(results_cox1)

}

pdf("/u/kisaev/TCGA_candidates_survival_plots_final_cands_FULL_10year_PFI_2019.pdf", width=6, height=5)
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$lnc_test_ph = as.numeric(tcga_results1$lnc_test_ph)
tcga_results1$global_test_ph = as.numeric(tcga_results1$global_test_ph)

tcga_results1$fdr_pval = p.adjust(as.numeric(tcga_results1$pval), method="fdr")
tcga_results1$fdr_anova_lr = p.adjust(as.numeric(tcga_results1$anova_pval), method="fdr")

tcga_results1 = as.data.table(tcga_results1)
#tcga_results1 = as.data.table(filter(tcga_results1, fdr_pval < 0.05))

tcga_results1$perc_wevents = as.numeric(tcga_results1$perc_wevents)
tcga_results1$num_events = as.numeric(tcga_results1$num_events)
tcga_results1$lnc_better = ""
z= which(tcga_results1$lnc_only_concordance >= tcga_results1$clinical_only_concordance)
tcga_results1$lnc_better[z] = "yes"

tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval)]

#check which models violate the PH assumption
#to those models add age * survival time interaction 
which(tcga_results1$global_test_ph <= 0.05)
tcga_results1$num_risk = as.numeric(tcga_results1$num_risk)

#plot distribution, cut number of risk patients 
tcga_results1$groupy = cut(tcga_results1$num_risk, breaks =c(1, 20, 40, 60, 80, 100, 200, 300, 
  400, 500, 600, 700, 800, 900, 1000))

tcga_results1$perc_risk = as.numeric(tcga_results1$perc_risk)
#tcga_results1 = filter(tcga_results1, fdr_pval <=0.05)
tcga_results1$gene_name = sapply(tcga_results1$gene, get_name)
saveRDS(tcga_results1, file="TCGA_results_multivariate_results_Oct3_PFI.rds")

colnames(fantom)[1] = "gene"
tcga_results1 = merge(fantom, tcga_results1, by="gene")
write.table(tcga_results1, file="SuppTable4_PFI.txt", quote=F, row.names=F, sep=";")



