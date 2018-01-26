###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:5790)]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

library(glmnet)
library(survcomp)
library(caret)
library(stringr)

###---------------------------------------------------------------
###Results - using all ovarian cancer patients ~305
###---------------------------------------------------------------

#1. 1000 runs batch1
##cindeces
genes1 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch1.rds")
genes1 = unlist(genes1)
genes1 = as.data.table(table(genes1))
genes1 = genes1[order(N)]
genes1 = dplyr::filter(genes1, N>=400)
genes1$batch = 1
colnames(genes1)[1] = "gene"

genes2 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch2.rds")
genes2 = unlist(genes2)
genes2 = as.data.table(table(genes2))
genes2 = genes2[order(N)]
genes2 = dplyr::filter(genes2, N>=400)
genes2$batch = 2
colnames(genes2)[1] = "gene"


genes3 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch3.rds")
genes3 = unlist(genes3)
genes3 = as.data.table(table(genes3))
genes3 = genes3[order(N)]
genes3 = dplyr::filter(genes3, N>=400)
genes3$batch = 3
colnames(genes3)[1] = "gene"


genes4 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch4.rds")
genes4 = unlist(genes4)
genes4 = as.data.table(table(genes4))
genes4 = genes4[order(N)]
genes4 = dplyr::filter(genes4, N>=400)
genes4$batch = 4
colnames(genes4)[1] = "gene"


genes5 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch5.rds")
genes5 = unlist(genes5)
genes5 = as.data.table(table(genes5))
genes5 = genes5[order(N)]
genes5 = dplyr::filter(genes5, N>=400)
genes5$batch = 5
colnames(genes5)[1] = "gene"


genes6 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch6.rds")
genes6 = unlist(genes6)
genes6 = as.data.table(table(genes6))
genes6 = genes6[order(N)]
genes6 = dplyr::filter(genes6, N>=400)
genes6$batch = 6
colnames(genes6)[1] = "gene"


genes7 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch7.rds")
genes7 = unlist(genes7)
genes7 = as.data.table(table(genes7))
genes7 = genes7[order(N)]
genes7 = dplyr::filter(genes7, N>=400)
genes7$batch = 7
colnames(genes7)[1] = "gene"


genes8 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch8.rds")
genes8 = unlist(genes8)
genes8 = as.data.table(table(genes8))
genes8 = genes8[order(N)]
genes8 = dplyr::filter(genes8, N>=400)
genes8$batch = 8
colnames(genes8)[1] = "gene"


genes9 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch9.rds")
genes9 = unlist(genes9)
genes9 = as.data.table(table(genes9))
genes9 = genes9[order(N)]
genes9 = dplyr::filter(genes9, N>=400)
genes9$batch = 9
colnames(genes9)[1] = "gene"


genes10 = readRDS("ALL_OV_pats305_binary_predictors_list_of_sig_genes_CV1000Jan23_batch10.rds")
genes10 = unlist(genes10)
genes10 = as.data.table(table(genes10))
genes10 = genes10[order(N)]
genes10 = dplyr::filter(genes10, N>=400)
genes10$batch = 10
colnames(genes10)[1] = "gene"

genes = rbind(genes1, genes2, genes3, genes4, genes5, genes6, genes7, genes8, genes9, genes10)

saveRDS(genes, file="final_candidates_10batches_of1000CV_305ovarian_cancer_patientsJan23.RDS")


















