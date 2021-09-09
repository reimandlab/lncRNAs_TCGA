library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands


#--------This script ------------------------------------------------

#look at number of lncs with uneven split and even split
#seperate by HRs

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

all_results = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
all_results = as.data.table(filter(all_results, data == "TCGA"))

#remove insignificant ones
all_results = as.data.table(filter(all_results, fdr_pval <= 0.05))

#-------------------PLOT---------------------------------------------

all_results$exp_type = ""
all_results$exp_type[all_results$HR > 1] = "Unfavourable"
all_results$exp_type[all_results$HR < 1] = "Favourable"

all_results$split_type = ""
all_results$perc_risk = round(as.numeric(all_results$perc_risk), digits=2)
all_results$split_type[all_results$perc_risk == 0.5] = "Even"
all_results$split_type[all_results$perc_risk < 0.5] = "Uneven_smaller_risk_group"
all_results$split_type[all_results$perc_risk > 0.5] = "Uneven_bigger_risk_group"

#summary
summ = as.data.table(table(all_results$exp_type, all_results$split_type))
colnames(summ) = c("Expression_type", "Split_type", "Num_lncRNAs")
summ$perc = summ$Num_lncRNAs/170

pdf("summary_lncRNAs_risk_ratio_split.pdf")
g = ggbarplot(summ, x="Split_type", y="perc", fill="Expression_type", palette=c("lightcyan4", "darksalmon"))
ggpar(g, font.xtickslab = c(6, "plain", "black"), legend = "right", legend.title="Expression Risk") + xlab("Split Type") + ylab("Percent of lncRNAs")

dev.off()

#by cancer type

summ = as.data.table(table(all_results$exp_type, all_results$split_type, all_results$cancer))
summ = as.data.table(filter(summ, N >0))
colnames(summ) = c("Expression_type", "Split_type", "Cancer", "Num_lncRNAs")
#summ$perc = summ$Num_lncRNAs/170

cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]

summ = merge(summ, cancers_conv, by="Cancer")

pdf("summary_by_cancer_type_lncRNAs_risk_ratio_split.pdf", width=10)
g = ggbarplot(summ, x="Split_type", y="Num_lncRNAs", fill="Expression_type", facet.by="type", palette=c("lightcyan4", "darksalmon"))
ggpar(g, font.xtickslab = c(6, "plain", "black"), legend = "right", legend.title="Expression Risk", x.text.angle = 90) + xlab("Split Type") + ylab("Percent of lncRNAs")

dev.off()
