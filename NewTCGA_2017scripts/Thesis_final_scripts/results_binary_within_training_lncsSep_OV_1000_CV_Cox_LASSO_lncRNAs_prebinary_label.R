###source_file.R

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 
lnc_info = read.csv("fantom_genebased_evidence_supp_table_17.csv")
lnc_info = lnc_info[which(str_detect(lnc_info$CAT_geneID, "ENSG")),]
#lnc_info = subset(lnc_info, lnc_info$CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_divergent", "p_lncRNA_intergenic"))
#shorten the gene names 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
lnc_info[,1] <- apply(lnc_info[,1:2], 1, extract3) #5049 lncRNAs 
fantom = lnc_info
colnames(fantom)[1] = "gene"

#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(limma)
library(stringr)
library(scater)
require(maftools)

library(gridExtra)
library(grid)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#3. ov
ov_cinds_clin = readRDS(file="8020_OV_100CV_justclin_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
ov_cinds_justlncs = readRDS(file="8020_OV_100CV_justlncs_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
ov_cinds_combined = readRDS(file="8020_OV_100CV_combined_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
ov_genes_results = readRDS(file="8020_OV_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
ov_cinds_clin = as.data.frame(ov_cinds_clin)
colnames(ov_cinds_clin) = "cindex"
ov_cinds_clin$predictor = "clinical"
ov_cinds_combined = as.data.frame(ov_cinds_combined)
colnames(ov_cinds_combined) = "cindex"
ov_cinds_combined$predictor = "combined"
####
ov_cinds_justlncs = as.data.frame(ov_cinds_justlncs)
colnames(ov_cinds_justlncs) = "cindex"
ov_cinds_justlncs$predictor = "justlncs"
####
ov = rbind(ov_cinds_clin, ov_cinds_combined, ov_cinds_justlncs)
ov$canc = "ov"

#Make box plots for cindices#-----------------------------
#all_cancers = rbind(ov, ov, paad, ov)
pdf("8020_split_1000CVs_ov_march20th_pval0.05_pcawgdet_correlated_lncs.pdf", width=9, height=9)
ggboxplot(ov, x="predictor", y="cindex", fill="predictor", palette=mypal) +
stat_boxplot(geom = "errorbar", width = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group="clinical") +
geom_hline(yintercept = 0.5, linetype = 2, colour="red")
dev.off()

ov_features = as.data.table(table(unlist(ov_genes_results)))
ov_features = ov_features[order(N)]
ov_features = dplyr::filter(ov_features, N >=500)
ov_features$canc = "ov"
ov_features$name = ""
for(i in 1:nrow(ov_features)){
	z = which(fantom$gene == ov_features$V1[i])
	ov_features$name[i] = fantom$CAT_geneName[z]
}

