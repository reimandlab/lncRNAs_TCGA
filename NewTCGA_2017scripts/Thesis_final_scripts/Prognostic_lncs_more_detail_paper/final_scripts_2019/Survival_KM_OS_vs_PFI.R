setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

library(survAUC)
require(caTools)
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
library(SIBER)
library(EnvStats)
source("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/universal_LASSO_survival_script.R")

#------FEATURES-----------------------------------------------------

#cands -- should be this file
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

pfi = fread("PFI_candidates_survival_times.txt")
colnames(pfi)[7:28] = paste("PFI", colnames(pfi)[7:28], sep="_")
os = fread("SuppTable4.txt")
pfi_cancers = c("LGG", "BRCA", "READ")
merged = merge(os, pfi)
colnames(canc_conv)[2]="cancer"
merged=merge(merged, canc_conv, by="cancer")
merged$plot = ""
merged$plot[merged$type == "LGG"] = "LGG"
merged$plot[merged$type == "READ"] = "READ" 
merged$plot[merged$type == "BRCA"] = "BRCA" 
merged$plot = factor(merged$plot, levels=c("LGG", "", "BRCA", "READ"))

merged$HR  = log2(merged$HR)
merged$PFI_HR = log2(merged$PFI_HR)
z = which((merged$fdr_pval < 0.05) & (merged$PFI_fdr_pval < 0.05))
merged$significant = ""
merged$significant[z] = "sig"

#1. get correlation between hazard ratios OS vs PFI and highlight correlation for LGG, READ and BRCA 

pdf("/u/kisaev/figure4X_PFI_vs_OS.pdf", width=6, height=6)
ggscatter(merged, x = "HR", y = "PFI_HR", palette="npg", alpha = 0.7, 
          add = "reg.line",  color="plot",  shape="significant",                             # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
  stat_cor(method = "spearman") + theme_bw() + xlim(c(-3,3)) +  ylim(c(-3,3)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") + xlab("log2HR (OS)")+ylab("log2HR (PFI)") 
dev.off()

merged = merged[,c(1:13,27,30:34, 49, 52)]
write.csv(merged, "/u/kisaev/supp_table_input_XX_PFI_wOStimes.csv", quote=F, row.names=F)


#2. get numbers of lncs that remain significant in PFI for LGG, READ and BRCA

#save supplementary table 

