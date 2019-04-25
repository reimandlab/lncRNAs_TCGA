
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#author: Karina Isaev, karin.isaev@gmail.com 
#date updated: Sept 21, 2018
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#-------------------------------------------------------------------
#this script uses data from TCGA (gene expression and clinical) to evaluate 
#the prognositc value of ion channels across different cancer types 
#working directory (source files, data files)
#/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#this script below prepares the RNA and clinical files for analysis 
source("universal_LASSO_survival_script_oncochannels.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(EnvStats)

#source("check_lnc_exp_cancers.R")

#------DATA---------------------------------------------------------

dim(rna)
print(table(rna$type))

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
colnames(ucsc)[8] = "HGNC.symbol"
cands = merge(cands, ucsc, by = "HGNC.symbol")

#--------This script ------------------------------------------------

#first assign outlier based binary labels to everyone
#then conduct outlier based survival analysis 

#--------------------------------------------------------------------

#1. Get cancer data (gene expression and clinical for each cancer type)
all = rna

t = as.data.table(table(all$type))
t = filter(t, N >=50)

cancers = unique(t$V1) #look at only cancers with minimum 100 patients 

#just lgg 
cancers = c("LGG", "GBM")

get_canc_data = function(canc){
  sub = subset(all, type == canc)
  return(sub)
} 
cancer_data = llply(cancers, get_canc_data)

#2. Subset data to ion channels and clinical variables 

z = which(cands$hg19.ensGene.name2 %in% c("ENSG00000165474", "ENSG00000169432", "ENSG00000175294", "ENSG00000103569"))
cands = cands[z,]

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

add_tags = function(dtt){
  print(dtt$type[1])
  rownames(dtt) = dtt$patient
  dtt$patient = NULL

  #log1p 
  z = which(str_detect(colnames(dtt), "ENSG"))
  dtt[,z] = log1p(dtt[,z])

  if(length(z)>1){
  medians = apply(dtt[,z], 2, median)}
  if(length(z)==1){
    medians = median(dtt[,z])
  }
  #add high low tag but first make boxplot
  rm = c()
  for(k in 1:length(medians)){
    med = medians[k]
    if(med == 0){rm = c(rm, names(medians[k]))}
    if(!(med ==0)){
    
    l1 = which(dtt[,names(medians)[k]] >= med)
    l2 = which(dtt[,names(medians)[k]] < med)
    dtt$med_tag = ""
    dtt$med_tag[l1] = "High"
    dtt$med_tag[l2] = "Low"

    z = which(colnames(dtt) == names(medians)[k])
    dtt$gene = dtt[,z]
    dtt$med_tag = factor(dtt$med_tag, levels = c("Low", "High"))
    
    file = paste(names(medians)[k], dtt$type[1], "median split.png", sep="_")
    png(file, units="in", width=5, height=5, res=300)
    g = ggboxplot(dtt, "med_tag", "gene",
     fill = "med_tag", palette = c("#00AFBB", "#FC4E07"), title=paste(names(medians)[k], dtt$type[1], "median split"))+
    xlab("Median tag") + ylab("log1p(FPKM-UQ)") + stat_n_text()
    print(g)
    dev.off()

    dtt[l1,names(medians)[k]] = 1
    dtt[l2, names(medians)[k]] = 0
    }
  }
  z = which(colnames(dtt) %in% rm)
  if(!(length(z)==0)){
  dtt = dtt[,-z]}  
  return(dtt)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")

#saveRDS(filtered_data_tagged, file="29cancers_median_tagged_by_ion_channels.rds")

#3. Fit survival models and plot Kaplan Meier plot 

get_survival_models = function(dtt){
  dtt$gene = NULL
  dtt$med_tag = NULL

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
  z = which(colnames(dat) == "age_at_initial_pathologic_diagnosis")
  if(length(z)==1){  
    dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)}

  num_genes = which(str_detect(colnames(dat), "ENSG"))

  #save KM plots for each lncRNA for each cancer type sepereatley 
  
  for(i in 1:length(num_genes)){
  gene = num_genes[i]  
  file = paste(colnames(dat)[num_genes[i]], dtt$type[1], "median_split_km_plot.png", sep="_")
  png(file, units="in", width=5, height=5, res=300)

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
          #surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 9,
          risk.table.fontsize = 4, 
          legend.labs = c("Low Exp", "High Exp"),             # survfit object with calculated statistics.
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
          ggtheme = theme_classic(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
          dev.off()
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
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test")
  #print(p)
  print(dtt$type[1])
}
}
#dev.off()

results_cox1 = results_cox1[-1,]
#fdr on p-values 
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")

return(results_cox1)

}

tcga_results = llply(filtered_data_tagged, get_survival_models) 
#saveRDS(tcga_results, "median_based_IC_survival_analysis.R")

#for now just need GBM 
#the other cancer types shouldn't be affected 
#before i just didnt have the right number of patients 
#gbm_tagged = filtered_data_tagged[[1]]
#tcga_results = get_survival_models(gbm_tagged)

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results)
tcga_results1$ic_test_ph = as.numeric(tcga_results1$ic_test_ph)
tcga_results1$global_test_ph = as.numeric(tcga_results1$global_test_ph)
tcga_results1$fdr_pval = as.numeric(tcga_results1$fdr_pval)
tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval)]

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

get_name_pcg = function(pcg){
  z = which(ucsc$hg19.ensGene.name2 == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensemblToGeneName.value[z])
}

tcga_results1$name = sapply(tcga_results1$gene, get_name_pcg)

saveRDS(tcga_results1, file="TCGA_ION_CHANNEL_results_april11_median_based.rds")
#saveRDS(tcga_results1, file="GBM_median_splits_IonCHannels.rds")




