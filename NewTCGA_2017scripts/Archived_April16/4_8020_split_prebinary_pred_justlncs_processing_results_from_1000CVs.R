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

#Cancers#-----------------------------------------------

#1. LIHC
lihc_cinds_clin = readRDS(file="8020_LIHC_100CV_justclin_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
lihc_cinds_justlncs = readRDS(file="8020_LIHC_100CV_justlncs_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
lihc_cinds_combined = readRDS(file="8020_LIHC_100CV_combined_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
lihc_genes_results = readRDS(file="8020_LIHC_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
lihc_cinds_clin = as.data.frame(lihc_cinds_clin)
colnames(lihc_cinds_clin) = "cindex"
lihc_cinds_clin$predictor = "clinical"
lihc_cinds_combined = as.data.frame(lihc_cinds_combined)
colnames(lihc_cinds_combined) = "cindex"
lihc_cinds_combined$predictor = "combined"
####
lihc_cinds_justlncs = as.data.frame(lihc_cinds_justlncs)
colnames(lihc_cinds_justlncs) = "cindex"
lihc_cinds_justlncs$predictor = "justlncs"
####
lihc = rbind(lihc_cinds_clin, lihc_cinds_combined, lihc_cinds_justlncs)
lihc$canc = "lihc"

#2. OV
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

#4. PAAD
paad_cinds_clin = readRDS(file="PAAD_8020_100CV_justclin_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
paad_cinds_justlncs = readRDS(file="PAAD_8020_100CV_justlncs_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
paad_cinds_combined = readRDS(file="PAAD_8020_100CV_combined_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
paad_genes_results = readRDS(file="PAAD_8020_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
paad_cinds_clin = as.data.frame(paad_cinds_clin)
colnames(paad_cinds_clin) = "cindex"
paad_cinds_clin$predictor = "clinical"
paad_cinds_combined = as.data.frame(paad_cinds_combined)
colnames(paad_cinds_combined) = "cindex"
paad_cinds_combined$predictor = "combined"
####
paad_cinds_justlncs = as.data.frame(paad_cinds_justlncs)
colnames(paad_cinds_justlncs) = "cindex"
paad_cinds_justlncs$predictor = "justlncs"
####
paad = rbind(paad_cinds_clin, paad_cinds_combined, paad_cinds_justlncs)
paad$canc = "paad"

#Make box plots for cindices#-----------------------------
all_cancers = rbind(lihc, kirc, paad, ov)
pdf("8020_split_1000CVs_allcancers_mar22.pdf", width=9, height=9)
ggboxplot(all_cancers, x="predictor", y="cindex", facet.by="canc", fill="predictor", palette=mypal) +
stat_boxplot(geom = "errorbar", width = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group="clinical") +
geom_hline(yintercept = 0.5, linetype = 2, colour="red")
dev.off()

#Which features were selected#----------------------------
lihc_features = as.data.table(table(unlist(lihc_genes_results)))
lihc_features = lihc_features[order(N)]
lihc_features = dplyr::filter(lihc_features, N >=500)
lihc_features$canc = "lihc"
lihc_features$name = ""
for(i in 1:nrow(lihc_features)){
	z = which(fantom$gene == lihc_features$V1[i])
	lihc_features$name[i] = fantom$CAT_geneName[z]
}

ov_features = as.data.table(table(unlist(ov_genes_results)))
ov_features = ov_features[order(N)]
ov_features = dplyr::filter(ov_features, N >=500)
ov_features$canc = "ov"
ov_features$name = ""
for(i in 1:nrow(ov_features)){
	z = which(fantom$gene == ov_features$V1[i])
	ov_features$name[i] = fantom$CAT_geneName[z]
}

kirc_features = as.data.table(table(unlist(kirc_genes_results)))
kirc_features = kirc_features[order(N)]
kirc_features = dplyr::filter(kirc_features, N >=500)
kirc_features$canc = "kirc"
kirc_features$name = ""
for(i in 1:nrow(kirc_features)){
	z = which(fantom$gene == kirc_features$V1[i])
	kirc_features$name[i] = fantom$CAT_geneName[z]
}


paad_features = as.data.table(table(unlist(paad_genes_results)))
paad_features = paad_features[order(N)]
paad_features = dplyr::filter(paad_features, N >=500)
paad_features$canc = "paad"
paad_features$name = ""
for(i in 1:nrow(paad_features)){
	z = which(fantom$gene == paad_features$V1[i])
	paad_features$name[i] = fantom$CAT_geneName[z]
}

all_chosen_features = rbind(lihc_features, ov_features, kirc_features, paad_features)
#saveRDS(all_chosen_features, file="chosen_features_all_cancesr_Mar22_1000CVs_8020split.rds")

chosen_features = readRDS("chosen_features_all_cancesr_Mar22_1000CVs_8020split.rds")
colnames(chosen_features)[1:3] =c("gene", "numChosen", "Cancer")

chosen_features$TCGA_HR = as.numeric(chosen_features$TCGA_HR)
chosen_features$TCGA_HR = round(chosen_features$TCGA_HR, digits=3)

chosen_features$TCGA_pval = as.numeric(chosen_features$TCGA_pval)
chosen_features$TCGA_pval = round(chosen_features$TCGA_pval, digits=6)

chosen_features$PCAWG_HR = as.numeric(chosen_features$PCAWG_HR)
chosen_features$PCAWG_HR = round(chosen_features$PCAWG_HR, digits=3)

chosen_features$PCAWG_pval = as.numeric(chosen_features$PCAWG_pval)
chosen_features$PCAWG_pval = round(chosen_features$PCAWG_pval, digits=4)

chosen_features = merge(chosen_features, fantom, by="gene")
chosen_features = chosen_features[,c(1:8, 10:11, 13:17)]

chosen_features = as.data.table(chosen_features)
chosen_features = chosen_features[order(PCAWG_pval)]

saveRDS(chosen_features, file="chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")

pdf("chosen_justlncsaswell_features_1000CVs_wPCAWG_validation_Mar22.pdf", height=12, width=23)
p<-tableGrob(chosen_features)
grid.arrange(p)
dev.off()




















