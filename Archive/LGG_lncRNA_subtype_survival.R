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

#------DATA---------------------------------------------------------

#get full dataset of GBM patients 
#ext = readRDS("all_genes_external_tcga_all_cancers_March13_wclinical_data.rds")

#check if cands are significant using data from ext 
pats = as.data.table(table(ext$type))
pats = as.data.table(filter(pats, N >= 15))
colnames(pats)[1] = "type"

canc_conv = unique(rna[,c("type", "Cancer")])
pats = merge(pats, canc_conv, by="type")

#get gbm
gbm = subset(ext, type=="GBM")

z = which(colnames(pcg) %in% colnames(rna))
cols = colnames(pcg)[z]
all = merge(rna, pcg, by=cols)

z = which(colnames(all) %in% colnames(gbm))
cols = colnames(all)[z]
all = all[,z]

z = which(colnames(gbm) %in% colnames(all))
cols = colnames(gbm)[z]
gbm = gbm[,z]

r = rbind(all, gbm)

all = r

lgg_sub = readRDS("lgg_subtype_info_molecular.rds")
lgg_gene_exp = as.data.table(filter(all, type == "LGG"))
#hoxa10-as/ and other GBM one 
cols = c("ENSG00000240889", "patient", "OS", "OS.time", "type")
lgg_gene_exp = lgg_gene_exp[,..cols]

#label 
#gene = first column
gene_exp = as.numeric(unlist(lgg_gene_exp[,1]))
med = median(gene_exp)
if(med ==0){
  z1 = which(lgg_gene_exp[,1] ==0)
  z2 = which(lgg_gene_exp[,1] > 0)
  lgg_gene_exp[z1,1] = 0
  lgg_gene_exp[z2,1] = 1
}
if(!(med ==0)){
  z1 = which(lgg_gene_exp[,1] >= med)
  z2 = which(lgg_gene_exp[,1] < med)
  lgg_gene_exp[z1,1] = 1
  lgg_gene_exp[z2,1] = 0
}

lgg_gene_exp = merge(lgg_gene_exp, lgg_sub, by="patient")
z = which(is.na(lgg_gene_exp$IDH.status))
lgg_gene_exp = lgg_gene_exp[-z,]

z = which(is.na(lgg_gene_exp$TERT.promoter.status))
lgg_gene_exp = lgg_gene_exp[-z,]

z = which(is.na(lgg_gene_exp$TERT.expression.status))
lgg_gene_exp = lgg_gene_exp[-z,]

lgg_gene_exp$subtype = paste(lgg_gene_exp$MGMT.promoter.status, lgg_gene_exp$IDH.status, lgg_gene_exp$X1p.19q.codeletion, lgg_gene_exp$TERT.promoter.status, lgg_gene_exp$TERT.expression.status, 
  lgg_gene_exp$ATRX.status, sep="_")

#--------This script ------------------------------------------------

#HOXA10-AS

#--------------------------------------------------------------------


#combine clinical variables 
#first IDH splits people 

lgg_gene_exp$subtype[lgg_gene_exp$IDH.status == "WT"] = "IDH_WT"
lgg_gene_exp$subtype[lgg_gene_exp$IDH.status == "Mutant"] = "IDH_mut"

lgg_gene_exp$subtype[lgg_gene_exp$IDH.status == "Mutant"] = "IDH_mut"

#add codeletion
z1 = which((lgg_gene_exp$subtype == "IDH_mut") & (lgg_gene_exp$X1p.19q.codeletion == "non-codel"))
z2 = which((lgg_gene_exp$subtype == "IDH_mut") & (lgg_gene_exp$X1p.19q.codeletion == "codel"))
lgg_gene_exp$subtype[z1] = "idh_mut_nocodel"
lgg_gene_exp$subtype[z2] = "idh_mut_codel"

colnames(lgg_gene_exp)[2] = "lncRNA"

  #dat$OS = as.numeric(dat$OS)
  #dat$OS.time = as.numeric(dat$OS.time)
  #dat$OS.time = dat$OS.time/365

#check anova model wtih lncRNA versus one without
lnc = coxph(Surv(OS.time, OS) ~ subtype + lncRNA, data = lgg_gene_exp)
nolnc = coxph(Surv(OS.time, OS) ~ subtype, data = lgg_gene_exp)

anova(lnc, nolnc)

pdf("NDUFB2-AS1_lncRNA_lgg_subtypes_KM_curve.pdf", width=11)
    fit <- survfit(Surv(OS.time, OS) ~ subtype + lncRNA, data = lgg_gene_exp)
            s <- ggsurvplot(
            title = paste("LGG subtypes", lgg_gene_exp$type[1]),
            fit, 
            xlab = "Time (Years)", 
            surv.median.line = "hv",
            font.main = c(16, "bold", "black"),
            font.x = c(14, "plain", "black"),
            font.y = c(14, "plain", "black"),
            font.tickslab = c(14, "plain", "black"),
            font.legend = 7,
            risk.table.fontsize = 5, 
            #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
            data = lgg_gene_exp,      # data used to fit survival curves. 
            risk.table = TRUE,       # show risk table.
            legend = "right", 
            pval = TRUE,             # show p-value of log-rank test.
            conf.int = FALSE,        # show confidence intervals for 
                              # point estimaes of survival curves.
            #xlim = c(0,5),        # present narrower X axis, but not affect
                              # survival estimates.
            break.time.by = 1,     # break X axis in time intervals by 500.
            #palette = colorRampPalette(mypal)(14), 
            #palette = c("blue", "red"),
            ggtheme = theme_bw(), # customize plot and risk table with a theme.
            risk.table.y.text.col = T, # colour risk table text annotations.
            risk.table.y.text = FALSE # show bars instead of names in text annotations
                              # in legend of risk table
            )
            print(s)
dev.off()


