
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

##look into subtype rather than IDH effect 

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
rna = rna[-which(rna$Cancer %in% canc_rm),]

#Combined into one dataframe because need to get ranks 
com = colnames(pcg)[which(colnames(pcg) %in% colnames(rna))]
all <- merge(rna, pcg, by = com)

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
cands = merge(cands, ucsc, by = "HGNC.symbol")

#outlier gene exp 
load("expr_discr_cc.rsav")
gbm = expr_discr_cc[[1]]
lgg = expr_discr_cc[[27]]

#change colnames 
colnames(gbm) = sapply(colnames(gbm), function(x){gsub("\\.", "-", x)})
colnames(lgg) = sapply(colnames(lgg), function(x){gsub("\\.", "-", x)})

gbm = t(gbm)
lgg = t(lgg)

genes_test = c("CATSPER1", "SCN9A", "AQP9", "GJB2")
z = which(colnames(gbm) %in% genes_test)
gbm = gbm[,z]
z = which(colnames(lgg) %in% genes_test)
lgg = lgg[,z]

lgg = as.data.frame(lgg)
gbm = as.data.frame(gbm)

gbm$patient = rownames(gbm)
lgg$patient = rownames(lgg)

#add surv info 

#add IDH mutation status 
lgg_clin = readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")
gbm_clin = readRDS("TCGA_gbm_wsubtype_info_biolinks.rds")

lgg = merge(lgg_clin, lgg, by = "patient")
gbm = merge(gbm_clin, gbm, by = "patient")

z = which(colnames(gbm) %in% c(genes_test, "Original.Subtype", "ESTIMATE.immune.score", "ATRX.status", "IDH.status", "patient", "OS", "OS.time"))
gbm = gbm[,z]

lgg$Transcriptome.Subtype = as.character(lgg$Transcriptome.Subtype)
z = which(colnames(lgg) %in% c(genes_test, "Original.Subtype", "ESTIMATE.immune.score", "ATRX.status", "IDH.status", "patient", "OS", "OS.time", "Transcriptome.Subtype"))
lgg = lgg[,z]

gbm$type = "GBM"
lgg$type ="LGG"
colnames(gbm)[colnames(gbm) == "Original.Subtype"] = "Transcriptome.Subtype"

#-----Survival models----------------------------------------------------

canc_dats = list(lgg, gbm)

#these survival models using Juri's high/low labels 

get_surv = function(dat){

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)

  genes = as.list(genes_test)

  get_ic_km = function(gene){

  print(gene)
  #get data for gene 
  gene_dat = dat
  z = which(colnames(gene_dat) == gene)
  colnames(gene_dat)[z] = "IC" 
    
  #Surv model just Ion channel 
  ic_model = coxph(Surv(OS.time, OS) ~ IC, data = gene_dat)
  ic_conc = glance(ic_model)$concordance
  ic_hr = summary(ic_model)$coefficients[2]
  ic_pval = summary(ic_model)$coefficients[5]
  print(ic_model)

  fit <- survfit(Surv(OS.time, OS) ~ IC, data = gene_dat)
          s <- ggsurvplot(
          title = paste(gene, gene_dat$type[1], "\nConcordance=", round(ic_conc, digits=3), "\nHR=", round(ic_hr, digits=3)),
          fit, 
          xlab = "Time (days)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = gene_dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          #xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          #break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

  #Full model with covariates
  ic_idh_model = coxph(Surv(OS.time, OS) ~ IC + Transcriptome.Subtype, data = gene_dat)
  ic_full_conc = glance(ic_idh_model)$concordance
  ic_full_hr = summary(ic_idh_model)$coefficients[1,2]
  ic_full_pval = summary(ic_idh_model)$coefficients[1,5]
  print(ic_idh_model)

  #KM plot with summary 
  fit <- survfit(Surv(OS.time, OS) ~ IC + Transcriptome.Subtype, data = gene_dat)
          s <- ggsurvplot(
          title = paste(gene, gene_dat$type[1], "\nConcordance=", round(ic_full_conc, digits=3), "\nHR=", round(ic_full_hr, digits=3)),
          fit, 
          xlab = "Time (days)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = gene_dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          #xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          #break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

    #compare full model to model without IDH 
    idh_only_model = coxph(Surv(OS.time, OS) ~ Transcriptome.Subtype, data = gene_dat)

    anov_ps = round(anova(idh_only_model, ic_idh_model)[2,4], digits=4)
    res = c(gene, gene_dat$type[1], ic_conc, ic_hr, ic_pval, ic_full_conc, ic_full_hr, ic_full_pval , anov_ps)
    names(res) = c("gene", "cancer", "ic_conc", "ic_hr", "ic_pval", "ic_full_conc", "ic_full_hr", "ic_full_pval" , "anov_ps")
    print(anov_ps)
    print(res)

    return(res)
  
  }
    
  gene_results = llply(genes, get_ic_km)
  gene_results = ldply(gene_results)
  return(gene_results)

}


pdf("lgg_gbm_ics_wSubtypes_km_plots.pdf", width=15, height=10)
results = llply(canc_dats, get_surv, .progress="text")
dev.off()


results = as.data.table(ldply(results))

results$ic_conc = as.numeric(results$ic_conc)
results$ic_hr = as.numeric(results$ic_hr)
results$ic_pval = as.numeric(results$ic_pval)
results$ic_full_conc = as.numeric(results$ic_full_conc)
results$ic_full_hr = as.numeric(results$ic_full_hr)
results$ic_full_pval = as.numeric(results$ic_full_pval)

results$ic_conc = round(results$ic_conc, digits=4)
results$ic_hr = round(results$ic_hr, digits=4)
results$ic_pval = round(results$ic_pval, digits=4)
results$ic_full_conc = round(results$ic_full_conc, digits=4)
results$ic_full_hr = round(results$ic_full_hr, digits=4)
results$ic_full_pval = round(results$ic_full_pval, digits=4)


write.csv(results, file="outlier_survival_subtypes_4exp_ion_channels_results_310119.csv", quote=F, row.names=F)




