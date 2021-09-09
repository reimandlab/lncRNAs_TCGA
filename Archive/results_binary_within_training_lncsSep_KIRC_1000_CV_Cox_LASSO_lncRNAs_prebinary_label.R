###source_file.R

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

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

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

#3. KIRC
kirc_cinds_clin = readRDS(file="KIRC_100CV_justclin_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
kirc_cinds_justlncs = readRDS(file="KIRC_100CV_justlncs_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
kirc_cinds_combined = readRDS(file="KIRC_100CV_combined_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
kirc_genes_results = readRDS(file="KIRC_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
kirc_cinds_clin = as.data.frame(kirc_cinds_clin)
colnames(kirc_cinds_clin) = "cindex"
kirc_cinds_clin$predictor = "clinical"
kirc_cinds_combined = as.data.frame(kirc_cinds_combined)
colnames(kirc_cinds_combined) = "cindex"
kirc_cinds_combined$predictor = "combined"
####
kirc_cinds_justlncs = as.data.frame(kirc_cinds_justlncs)
colnames(kirc_cinds_justlncs) = "cindex"
kirc_cinds_justlncs$predictor = "justlncs"
####
kirc = rbind(kirc_cinds_clin, kirc_cinds_combined, kirc_cinds_justlncs)
kirc$canc = "kirc"

#Make box plots for cindices#-----------------------------
#all_cancers = rbind(lihc, kirc, paad, ov)
pdf("1000CVs_KIRC_march20th_fdr0.05_pcawgdet_correlated_lncs.pdf", width=9, height=9)
ggboxplot(kirc, x="predictor", y="cindex", fill="predictor", palette=mypal) +
stat_boxplot(geom = "errorbar", width = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group="clinical") +
geom_hline(yintercept = 0.5, linetype = 2, colour="red")
dev.off()

kirc_features = as.data.table(table(unlist(kirc_genes_results)))
kirc_features = kirc_features[order(N)]
kirc_features = dplyr::filter(kirc_features, N >=500)
kirc_features$canc = "kirc"
kirc_features$name = ""
for(i in 1:nrow(kirc_features)){
	z = which(fantom$gene == kirc_features$V1[i])
	kirc_features$name[i] = fantom$CAT_geneName[z]
}

