set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FEATURES-----------------------------------------------------

#cands -- should be this file
#cands = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")

cands = readRDS("lncRNAs_selected_by_EN_april14.rds") #1000 runs of cross-validations using new updated dataset (GBM=150, OV and LUAD)

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
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course",
    "new_tumor_event_type", "Cancer", "type"))
  print(dtt$type[1])
  dtt = as.data.table(dtt)
  dtt = dtt[,..z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)
all_tcga_cancers = as.data.table(ldply(filtered_data))
all_tcga_cancers = unique(all_tcga_cancers[,c("patient", "type", "age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "new_tumor_event_type", "treatment_outcome_first_course", "OS", "OS.time", "Cancer")])
write.csv(all_tcga_cancers, file="TCGA_clinical_data_patients_used_external_survival_analysis.csv", quote=F, row.names=F)

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

  results_cox1 <- as.data.frame(matrix(ncol=28)) ; colnames(results_cox1) <- c("gene", "coef", "pval", "HR", "low95", "upper95", "cancer",
    "lnc_test_ph", "num_risk", "perc_risk", "median_nonzero", "sd_nonzero", "min_nonzero", "max_nonzero", "multi_model_concordance",
    "lnc_only_concordance", "clinical_only_concordance",
    "num_events", "perc_wevents", "anova_pval", "gene_symbol", "canc_type",
  "hr_adjusted", "pval_adjusted", "HR_adj_low95", "HR_adj_high95", "perc_zeroes", "median_exp")

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
  #num events
  num_events = length(which(dat$OS == 1))
  perc_events = num_events/nrow(dtt)

  for(i in 1:length(num_genes)){
    gene = num_genes[i]
    print(gene)
    k = which(!(str_detect(colnames(dat), "ENSG")))
    newdat = dat[,c(gene,k)]

    #get number of people with zero expression of this gene
    zeros_matrix = filter(rna, type == dtt$type[1])
    z = which(colnames(zeros_matrix) == colnames(newdat)[1])
    exp_of_gene = zeros_matrix[,..z]
    num_zeroes = length(which(exp_of_gene == 0))/dim(newdat)[1]

    #  newdat$gene=newdat[,gene]
    lncs = coxph(Surv(OS.time, OS)  ~ ., data = newdat)
    clin_vars_used = paste(colnames(newdat), collapse=",")
    print(clin_vars_used)

    #  global = test.ph$table[nrow(test.ph$table),3]

    #mutlvariate concordance
    cmulti = as.numeric(glance(lncs)[13])
    lnc_only_model = coxph(Surv(OS.time, OS)  ~ newdat[,1], data = newdat)
    lnc_only = as.numeric(glance(lnc_only_model)[13])
    #  hr = summary(lnc_only_model)$coefficients[2]
    test.ph <- cox.zph(lnc_only_model)
    lnc_test_ph = test.ph$table[1,3]

    z = which(str_detect(colnames(newdat), "ENSG"))
    clin_only = newdat[,-z]
    clinical_only_model = coxph(Surv(OS.time, OS)  ~ ., data = clin_only)
    clinical_only = as.numeric(glance(clinical_only_model)[13])
    lr = anova(clinical_only_model, lncs)
    lr_pval = lr[2,4]

    #calculate power
    #k is ratio of participants in group E (experimental group) compared to group C (controlgroup).
    kval = as.numeric(table(newdat[,1])[1]/table(newdat[,1])[2])
    #m = expected total number of events over both groups.
    mval = length(which(newdat$OS == 1))

    hr = summary(lnc_only_model)$coefficients[1,c(1,2,5)][2]
    pval = summary(lnc_only_model)$coefficients[1,c(1,2,5)][3]

    hr_adjusted = summary(lncs)$coefficients[1,c(1,2,5)][2]
    hr_adjust_low95=summary(lncs)$conf.int[1,3]
    hr_adjust_high95=summary(lncs)$conf.int[1,4]
    pval_adjusted = summary(lncs)$coefficients[1,c(1,2,5)][3]

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
    newdat$gene = factor(newdat$gene, levels=c(1,0))
    fit <- survfit(Surv(OS.time, OS) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(hg38$symbol[hg38$ensgene==gene_name], canc_conv$type[canc_conv$Cancer == dtt$Cancer[1]][1], "HR =", round(hr, digits=2)),
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

  gene_symbol=hg38$symbol[hg38$ensgene==gene_name]
  if((perc > 0.1) & (perc < 0.9)){
  print(s)}

   #generate boxplot
   z = which(rna$Cancer == dtt$Cancer[1])
   exp_data = rna[z,]
   exp_data = as.data.frame(exp_data)
   exp_data = exp_data[,which(colnames(exp_data) %in% c(gene, "patient"))]

   newdat$patient = rownames(newdat)
   exp_data = merge(exp_data, newdat, by="patient")
   colnames(exp_data)[2] = "geneexp"

   median_exp = median(exp_data$geneexp)

   #get median non-zero expression and SD
   if(median_exp == 0){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
   }

   if(!(median_exp == 0)){
     median_nonzero = "notavail"
     sd_nonzero = "notavail"
     min_nonzero = "notavail"
     max_nonzero = "notavail"
   }

   row <- c(gene_name, summary(lnc_only_model)$coefficients[1,c(1,5)], hr,
   summary(lnc_only_model)$conf.int[1,c(3,4)], dtt$Cancer[1],
    lnc_test_ph, risk_num, perc, median_nonzero,
    sd_nonzero,
    min_nonzero,
    max_nonzero, cmulti, lnc_only, clinical_only, num_events,
    perc_events, lr_pval, gene_symbol, canc_type=dtt$type[1], hr_adjusted, pval_adjusted,
    hr_adjust_low95, hr_adjust_high95, num_zeroes, median_exp)

   names(row) <- names(results_cox1)
   results_cox1 = rbind(results_cox1, row)

   #make density plot using FPKM-UQ values and logged values
   #visualize bimodality
   exp_data$median=""
   exp_data$median[exp_data$gene==0] = "Low"
   exp_data$median[exp_data$gene==1] = "High"
   exp_data$median = factor(exp_data$median, levels=c("High", "Low"))
   exp_data$geneexp = log1p(exp_data$geneexp)

   gg <- ggplot(exp_data)
   gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
   gg <- gg + theme_bw() + ggtitle(paste(gene_symbol, "Expression", dtt$Cancer[1] , sep=" ")) + labs(y="log1p(FPKM-UQ)")

#   print(gg)

   p <- ggboxplot(exp_data, x = "median", y = "geneexp",
          color = "median",
         palette = mypal[c(4,1)], title = paste(gene_symbol, "Expression", canc_conv$type[canc_conv$Cancer == dtt$Cancer[1]][1], sep=" "),
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_classic())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test") + stat_n_text() + scale_color_npg()
  if((perc > 0.1) & (perc < 0.9)){
  print(p)}
}

results_cox1 = results_cox1[-1,]
#fdr on p-values
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")
results_cox1$fdr_pval_adjusted = p.adjust(results_cox1$pval_adjusted, method="fdr")
results_cox1$clin_vars = clin_vars_used

return(results_cox1)

}

pdf("/u/kisaev/Jan2021/TCGA_candidates_survival_plots_final_cands_FULL_10year_OS.pdf", width=6, height=5)
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$lnc_test_ph = as.numeric(tcga_results1$lnc_test_ph)
tcga_results1$fdr_pval = p.adjust(as.numeric(tcga_results1$pval), method="fdr")
tcga_results1$fdr_pval_adjusted = p.adjust(as.numeric(tcga_results1$pval_adjusted), method="fdr")

tcga_results1$fdr_anova_lr = p.adjust(as.numeric(tcga_results1$anova_pval), method="fdr")

tcga_results1 = as.data.table(tcga_results1)

tcga_results1$perc_wevents = as.numeric(tcga_results1$perc_wevents)
tcga_results1$num_events = as.numeric(tcga_results1$num_events)
tcga_results1$num_risk = as.numeric(tcga_results1$num_risk)
tcga_results1$perc_risk = as.numeric(tcga_results1$perc_risk)
tcga_results1 = as.data.table(filter(tcga_results1, fdr_pval < 0.05, fdr_pval_adjusted < 0.05))
tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval, fdr_pval_adjusted)]

riskplot = gghistogram(tcga_results1, x = "perc_risk",
   fill = "white", color="black",  palette = c("#00AFBB", "#E7B800")) + xlab("Percentage of Risk patients")+
ylab("Frequency")

#pdf("/u/kisaev/Jan2021/Distribution_percent_risk_patients_per_lncRNA.pdf", width=10)
#riskplot
#dev.off()

tcga_results1$lnc_test_ph = as.numeric(tcga_results1$lnc_test_ph)
saveRDS(tcga_results1, file="TCGA_results_multivariate_results_Oct3.rds")

colnames(fantom)[1] = "gene"
tcga_results1 = merge(fantom, tcga_results1, by="gene")
write.table(tcga_results1, file="/u/kisaev/Jan2021/SuppTable4.txt", quote=F, row.names=F, sep="}")

tcga_results1 = as.data.table(tcga_results1)

dim(filter(tcga_results1, HR >1))
tcga_results1$HR = as.numeric(tcga_results1$HR)
median(filter(tcga_results1, HR >1)$HR)

dim(filter(tcga_results1, HR <1))
tcga_results1$HR = as.numeric(tcga_results1$HR)
median(filter(tcga_results1, HR <1)$HR)

tcga_results1 = readRDS("TCGA_results_multivariate_results_Oct3.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")
tcga_results1$clinical_only_concordance = NULL
tcga_results1$num_events = NULL
tcga_results1$perc_wevents = NULL
tcga_results1 = as.data.table(tcga_results1)
saveRDS(tcga_results1, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
