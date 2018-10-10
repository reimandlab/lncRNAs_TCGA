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

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

cands = fread("FMRE_difexp_target_genes.txt")

#--------This script ------------------------------------------------

#Candidate gene survival analysis, Vol2 
#These are differentially expressed genes in mutated samples in PCAWG
#here we summarize the survival results for each of these candidates 
#within each cancer type 
#Author = Karina Isaev
#Date = July 3rd, 2018

#Can you please perform median dichotomized survival analysis for 
#each of these genes across all cancer types in TCGA?

#Two expected outputs:

###----1. a geom_tile plot of gene X cancer type, HR to color the tile, asterisks to show significance of 
#HR from CoxPH (. 0.1 * 0.05 ** 0.01 *** 0.001 **** 0.0001). Genes ranked by significance from table. 
#We should somehow color or annotate gene names in plot to show those with positive and negative coefficients 
#shown in the table (pos coefficient - gene up regulated with mutations). 

###----2. all gene X cancer type KM plots, one file per gene, p-value ranked. 


#------PROCESS CANDIDATES--------------------------------------------

#1. get genes ids 
clean_gene = function(ge){
  newge = unlist(strsplit(ge, "::"))[2]
  #get ensembl id
  z = which(ucsc$hg19.ensemblToGeneName.value == newge)
  newge = ucsc$hg19.ensGene.name2[z]
  return(newge)
}

cands$gene = unlist(llply(cands$gene, clean_gene))

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cancer_data = canc_datas_pcgs #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(cands$gene, "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = as.data.frame(dtt)
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
  results_cox1 <- as.data.frame(matrix(ncol=11)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "sample_size", "num_events")

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
  num_genes = which(str_detect(colnames(dat), "ENSG"))

  for(i in 1:length(num_genes)){
  gene = num_genes[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  #for now don't need to have all clinical variables 
  #newdat = dat[,c(gene,k)]
  k = which(colnames(dat) %in% c("OS", "OS.time"))
  newdat = dat[,c(gene, k)]  

  lncs = coxph(Surv(OS.time, OS)  ~ ., data = newdat)
  test.ph <- cox.zph(lncs)
  lnc_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  row <- c(colnames(newdat)[1], summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], dtt$Cancer[1], lnc_test_ph, global, 
    dim(dat)[1], length(which(dat$OS==1)))
    
  names(row) <- names(results_cox1) 
  results_cox1 = rbind(results_cox1, row)
  #print(ggcoxzph(test.ph))

}

results_cox1 = results_cox1[-1,]
#fdr on p-values 
results_cox1$pval = as.numeric(results_cox1$pval)
results_cox1$fdr_pval = p.adjust(results_cox1$pval, method="fdr")

return(results_cox1)

}

pdf("TCGA_candidates_survival_plots_final_cands_Oct10.pdf")
tcga_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
tcga_results1 = ldply(tcga_results, data.frame)
tcga_results1$lnc_test_ph = as.numeric(tcga_results1$lnc_test_ph)
tcga_results1$global_test_ph = as.numeric(tcga_results1$global_test_ph)
tcga_results1$fdr_pval = as.numeric(tcga_results1$fdr_pval)
tcga_results1 = as.data.table(tcga_results1)
tcga_results1 = tcga_results1[order(fdr_pval)]

#------rename gene IDs-----------
change_gene = function(ge){
  newge = ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 == ge)]
  return(newge)
}
tcga_results1$name = ""
tcga_results1$name = unlist(llply(tcga_results1$gene, change_gene))
saveRDS(tcga_results1, file="TCGA_FMRE_diff_exp_PCGS_survival_results_Oct10.rds")

#------plot survival curves-------

#-------------------------------------------------------
#DO NOT RUN---------------------------------------------
#-------------------------------------------------------

pdf("TCGA_FMRE_10_diff_Exp_cands_surv_plots_Oct10.pdf")

for(i in 1:nrow(tcga_results1)){
  canc = tcga_results1$cancer[i]
  gene = tcga_results1$gene[i]
  name = tcga_results1$name[i]
  hr = round(as.numeric(tcga_results1$HR[i]), digits=3)

  newdat = subset(pcg, Cancer %in% canc)
  newdat = as.data.frame(newdat) 
  z = which(colnames(newdat) %in% c(gene, "Cancer", "type", "OS", "OS.time"))
  newdat = newdat[,z]
  colnames(newdat)[1] = name
  gene = colnames(newdat)[1]
  colnames(newdat)[1] = "gene"
  newdat$OS = as.numeric(newdat$OS)
  newdat$OS.time = as.numeric(newdat$OS.time)
  newdat$OS.time = newdat$OS.time/365
  newdat$geneexp = newdat[,1]

  #add tags
  medians = median(newdat[,1])
  #add high low tag
  med = medians
    k = 1
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(newdat[,k] > 0)
    l2 = which(newdat[,k] ==0)
    newdat[l1,k] = 1
    newdat[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(newdat[,k] >= med)
    l2 = which(newdat[,k] < med)
    newdat[l1,k] = 1
    newdat[l2, k] = 0
    } 
  
 order = c(0, 1)
 newdat$gene = factor(newdat$gene, levels = order) 
 fit <- survfit(Surv(OS.time, OS) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene, newdat$type[1], "HR =", hr),
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
   newdat$geneexp = log1p(newdat$geneexp)
   p <- ggboxplot(newdat, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene, "Expression", newdat$type[1] , sep=" "), 
          add = "jitter", ylab = "FPKM",  ggtheme = theme_bw())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test")
  print(p)

}#end tcga_results1

dev.off()






