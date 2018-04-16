###source_file.R

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

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

#Cancers#-----------------------------------------------

#1. LIHC
lihc_cinds_clin= readRDS(file="ALL_LIHC_pats_prebin_predictors_clin_cindeces_CV_100_March16_justlncsaswell.rds")
lihc_cinds_justlncs = readRDS(file="ALL_LIHC_prebin_predictors_justlncs_cindeces_CV_100_March16_justlncsaswell.rds")
lihc_cinds_combined=readRDS(file="ALL_LIHC_prebin_predictors_combined_cindeces_CV_100_March16_justlncsaswell.rds")
lihc_genes_results=readRDS(file="ALL_LIHC_pats_prebin_predictors_list_of_sig_genes_CV_100March16_justlncsaswell.rds")
####
lihc_cinds_clin = as.data.frame(lihc_cinds_clin)
colnames(lihc_cinds_clin) = "cindex"
lihc_cinds_clin$predictor = "clinical"
####
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
ov_cinds_clin = readRDS(file="ALL_OV_pats_prebin_predictors_clin_cindeces_CV_100_March16_justlncsaswell.rds")
ov_cinds_justlncs = readRDS(file="ALL_OV_prebin_predictors_justlncs_cindeces_CV_100_March16_justlncsaswell.rds")
ov_cinds_combined = readRDS(file="ALL_OV_prebin_predictors_combined_cindeces_CV_100_March16_justlncsaswell.rds")
ov_genes_results = readRDS(file="ALL_OV_pats_predictors_list_of_sig_genes_CV_100March16_justlncsaswell.rds")
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
kirc_cinds_clin = readRDS(file="ALL_KIRC_pats_prebin_predictors_clin_cindeces_CV_100_March16_justlncsaswell.rds")
kirc_cinds_justlncs = readRDS(file="ALL_KIRC_prebin_predictors_justlncs_cindeces_CV_100_March16_justlncsaswell.rds")
kirc_cinds_combined = readRDS(file="ALL_KIRC_prebin_predictors_combined_cindeces_CV_100_March16_justlncsaswell.rds")
kirc_genes_results = readRDS(file="ALL_KIRC_pats_prebin_predictors_list_of_sig_genes_CV_100March16_justlncsaswell.rds")
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
paad_cinds_clin = readRDS(file="ALL_PAAD_pats_prebin_predictors_clin_cindeces_CV_100_March16_justlncsaswell.rds")
paad_cinds_justlncs = readRDS(file="ALL_PAAD_prebin_predictors_justlncs_cindeces_CV_100_March16_justlncsaswell.rds")
paad_cinds_combined = readRDS(file="ALL_PAAD_prebin_predictors_combined_cindeces_CV_100_March16_justlncsaswell.rds")
paad_genes_results = readRDS(file="ALL_PAAD_pats_prebin_predictors_list_of_sig_genes_CV_100March16_justlncsaswell.rds")
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
pdf("100CVs_lncsASwell_prebinary_labelled_predictor_cindecides_fixed_dichotomization.pdf", width=9, height=9)
ggboxplot(all_cancers, x="predictor", y="cindex", facet.by="canc", fill="predictor", palette=mypal) +
stat_boxplot(geom = "errorbar", width = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group="clinical") +
geom_hline(yintercept = 0.5, linetype = 2, colour="red")
dev.off()

#Which features were selected#----------------------------
lihc_features = as.data.table(table(unlist(lihc_genes_results)))
lihc_features = lihc_features[order(N)]
lihc_features = dplyr::filter(lihc_features, N >=50)
lihc_features$canc = "lihc"

ov_features = as.data.table(table(unlist(ov_genes_results)))
ov_features = ov_features[order(N)]
ov_features = dplyr::filter(ov_features, N >=50)
ov_features$canc = "ov"

kirc_features = as.data.table(table(unlist(kirc_genes_results)))
kirc_features = kirc_features[order(N)]
kirc_features = dplyr::filter(kirc_features, N >=50)
kirc_features$canc = "kirc"

paad_features = as.data.table(table(unlist(paad_genes_results)))
paad_features = paad_features[order(N)]
paad_features = dplyr::filter(paad_features, N >=50)
paad_features$canc = "paad"

all_chosen_features = rbind(lihc_features, ov_features, kirc_features, paad_features)
saveRDS(all_chosen_features, file="chosen_features_justlncsaswell_100CVs_prebinary_predictors_fixed_dichotomization_March16.rds")
colnames(all_chosen_features) =c("gene", "numChosen", "Cancer")
pdf("chosen_justlncsaswell_features_100CVs_prebinary_predictors_fixed_dichomtization_just_clinical_vs_combined.pdf", height=12)
p<-tableGrob(all_chosen_features)
grid.arrange(p)
dev.off()




















