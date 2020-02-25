library(survAUC)
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(EnvStats)
library(reshape2)

source("check_lnc_exp_cancers.R")

#------FEATURES-----------------------------------------------------

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
colnames(allCands)[7] = "Cancer"
allCands = merge(allCands, canc_conv, by="Cancer")
allCands$combo = paste(allCands$gene, allCands$type)

lgg_dataset = as.data.table(filter(all, type =="LGG"))

z = which(colnames(lgg_dataset) %in% c("ENSG00000253293", "ENSG00000253187"))
lgg_dataset = lgg_dataset[,..z]
lgg_dataset = log1p(lgg_dataset)

library("ggpubr")

pdf("HOXA10AS_HOXA10_correlations.pdf")
ggscatter(lgg_dataset, x = "ENSG00000253187", y = "ENSG00000253293", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "HOXA10-AS", ylab = "HOXA10")
dev.off()

