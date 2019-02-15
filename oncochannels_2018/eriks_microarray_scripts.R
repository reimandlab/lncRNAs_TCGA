###now have methylation files for each cancer 
###with probes that overlap any of the lncRNA candidates (might not be the cancer specific candidate)

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
gbm = as.data.table(filter(all, type == "GBM"))
dup = gbm$patient[which(duplicated(gbm$patient))]
z = which(gbm$patient %in% dup)
gbm = gbm[-z,]

#-------------------------------------------------------------------------

microarray = read.csv("IC_ONLY_TCGA_Affymetrix_RMA.csv")
surv = microarray[1:3,]
surv = t(surv)
colnames(surv) = surv[1,]
surv=surv[-1,]
surv = as.data.table(surv)
colnames(surv)[1:2] = c("OS", "OS.time")
gbm_surv = gbm[,c("OS", "OS.time", "patient")]
gbm_surv = gbm_surv[which(gbm_surv$patient %in% surv$patient),]
rm = which(surv$patient %in% gbm_surv$patient)
surv = surv[-rm, ]
surv = rbind(surv, gbm_surv)

microarray = microarray[3:nrow(microarray),]
colnames(microarray) = as.character(unlist(microarray[1,]))
microarray = microarray[-1,]
colnames(microarray)[1] = "gene"
microarray$gene = as.character(microarray$gene)

#1. survival function 

###EASY WAY TO MAKE KM PLOT
get_km_plot_microarray = function(gene){
  
  print(gene)
  #merge expression with survival data
  z = which(microarray[,1] == gene)
  if(!(length(z)==0)){
  exp = (t(as.data.frame(microarray[z,])))
  exp = as.data.frame(exp)
  exp$patient = ""
  exp$patient = rownames(exp)
  colnames(exp)[1] = gene
  exp = exp[-1,]
  exp = merge(exp, surv, by="patient")
  
  print("part 1")

  #split patients 
  exp[,2] = as.numeric(exp[,2])
  med = median(exp[,2])
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(exp[,2] > 0)
    l2 = which(exp[,2] ==0)
    exp[l1,2] = 1
    exp[l2, 2] = 0
    }

    if(!(med ==0)){
    l1 = which(exp[,2] >= med)
    l2 = which(exp[,2] < med)
    exp[l1,2] = 1
    exp[l2, 2] = 0
    }

  print("part 2")

  exp$OS = as.numeric(exp$OS)
  exp$OS.time = as.numeric(exp$OS.time)
  exp$OS.time = exp$OS.time/365
  colnames(exp)[2] = "gene"
  exp$gene = factor(exp$gene, levels = c(0,1))
  
  #balance check
  bal_check = (table(exp[,2])[1] >= 5) & (table(exp[,2])[2] >= 5)
  if(bal_check){
  cox_mod = coxph(Surv(OS.time, OS) ~ gene, data = exp)
  print(glance(cox_mod)$concordance)
  conc = round(glance(cox_mod)$concordance, digits=2)
  print("part 3")
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = exp)
          s <- ggsurvplot(
          title = paste(gene, "GBM Microarray", "\nConcordance=", conc),
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
          data = exp,      # data used to fit survival curves. 
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
          #print(s)
    print("part 4")
    row = c(gene, summary(cox_mod)$coefficients[2], glance(cox_mod)[8], conc)
    row = unlist(row)
    names(row) = c("gene", "HR", "pval", "conc")
    print(row)
    return(row)
}
}	
}

#make plots for just the 4 candidates 
#cands = as.list(c("AQP9", "CATSPER1", "GJB2", "SCN9A"))

#all ics
genes = microarray$gene

#pdf("gbm_microarray_ionchannels_4_cands_exp_nov8.pdf")
res = llply(genes, get_km_plot_microarray, .progress="text")
#dev.off()
res = ldply(res)
res$HR = as.numeric(res$HR)
res$pval = as.numeric(res$pval)
res$conc = as.numeric(res$conc)
res$fdr = p.adjust(res$pval, method="fdr")
res = as.data.table(res)
res = res[order(fdr)]
colnames(res)[2:ncol(res)] = paste(colnames(res)[2:ncol(res)], "microarray", sep="_")
saveRDS(res, file="GBM_ics_microarray_results_feb9.rds")








#order all LGG ion channels 
r = readRDS("TCGA_ION_CHANNEL_results_Sept21.rds")
r$name = unlist(llply(r$gene, get_name_pcg))
r$HR = as.numeric(r$HR)
#write.csv(r, file="272_ion_channels_w_survival_information_all_cancer_types_KI.csv", quote=F, row.names=F)

r = read.csv("272_ion_channels_w_survival_information_all_cancer_types_KI.csv")
r=as.data.table(r)
r$HR = as.numeric(r$HR)

#272 ion channels 
lgg = as.data.table(filter(r, cancer == "Brain Lower Grade Glioma"))
lgg$name = unlist(llply(lgg$gene, get_name_pcg))
lgg$HR = as.numeric(lgg$HR)
lgg$HR = log2(lgg$HR)

lgg = lgg[order(pval, -abs(HR))]
filter(lgg, name %in% cands)

#make facet for each tumour type 
#gene name on x-axis 
#hazard ratio on y-axis
z = which(r$name %in% cands)
r$exp_val = ""
#r$exp_val[z] = "V"
r$exp_val[z] = r$name[z]

r$fdr_pval = as.numeric(r$fdr_pval)
r = as.data.table(filter(r, fdr_pval < 0.05))
r$HR = log2(r$HR)

#cols = sample(mypal5, 13)
cols = c("#96897F" ,"#6F46E6" ,"#49B98D", "#EAD286", "#70A2A4" ,"#C86ED7" ,"#B2B47A",
 "#7C9BE1", "#CFA0E0" ,"#BBE6DF" ,"#786DDA" ,"#DEAEC7", "#64709B")

canc_conv = unique(rna[,c("type", "Cancer")])
colnames(canc_conv)[2] = "cancer"
r = merge(r, canc_conv, by="cancer")

r = r[order(cancer)]
r$cancer = factor(r$cancer, levels = unique(r$cancer))
z = which(duplicated(r$name))
r$name[z] = paste(r$name[z], "*", sep="")
z = which(duplicated(r$name))
r$name[z] = paste(r$name[z], "**", sep="")
z = which(duplicated(r$name))
r$name[z] = paste(r$name[z], "***", sep="")
z = which(duplicated(r$name))
r$name[z] = paste(r$name[z], "****", sep="")
z = which(duplicated(r$name))
r$name[z] = paste(r$name[z], "*****", sep="")
t = as.data.table(table(r$cancer))
t=t[order(-N)]

r[, cancFac := factor(cancer, levels=t$V1)]
setorder(r, cancFac)

r =r[order(type, -(HR))]

r$name = factor(r$name, levels = unique(r$name))
r$type = factor(r$type, levels = unique(r$type))

#cols2 = sample(mypal5, length(unique(r$type)))
cols2= c("#67E9D0", "#D94753", "#70A2A4", "#786DDA", "#7BEE95", "#C7CBE7", "#B2B47A",
"#BBE6DF" ,"#65B9E0", "#C0EC3E", "#96897F", "#AAE6B0", "#74DAE3", "#EAD286",
"#7C9BE1", "#61EA4F", "#D1EB7B", "#D97B8F", "#6F46E6", "#49B98D", "#E5C0A6",
"#DF4FA6", "#DEAEC7", "#CFA0E0")

#for each cancer type - make barplot of genes versus hazard ratios
pdf("summary_231_ion_channels_hazard_ratios_fdr0.05_faceted.pdf", width=12)
#pdf("summary_231_ion_channels_hazard_ratios_pval0.01_faceted.pdf", width=10)
ggplot(r, aes(name, HR, label = exp_val)) + geom_bar(stat = "identity",
     aes(fill = type)) + scale_fill_manual(values=cols2)+
     theme(legend.position="bottom") + labs(x="Tumour type", y="log2(Hazard-Ratio)")+
     theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
     geom_text(aes(label=exp_val), position=position_dodge(width=0.9), size=4)+facet_grid(~type, scales="free", space="free")+
      theme(strip.text.x = element_text(size=3))
dev.off()


pdf("summary_231_ion_channels_hazard_ratios_fdr0.05_nonfaceted.pdf", width=10)
#pdf("summary_231_ion_channels_hazard_ratios_pval0.01_nonfaceted.pdf", width=10)
ggplot(r, aes(name, HR, label = exp_val)) + geom_bar(stat = "identity",
     aes(fill = type)) + scale_fill_manual(values=cols2)+
     theme(legend.position="bottom") + labs(x="Tumour type", y="log2(Hazard-Ratio)")+
     theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
     geom_text(aes(label=exp_val), position=position_dodge(width=0.9), size=4)
dev.off()














