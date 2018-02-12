###methylation files 
###identify probes overlapping lncRNAs 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#############in command line###################################################################################################################################
head -n 1 OV.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt > OV_patients.txt

#in R
ov_pats = read.table("OV_patients.txt")
ov_pats = ov_pats[,3:ncol(ov_pats)]
ov_pats = as.character(unlist(ov_pats))
ov_meth = fread("OV.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt")
colnames(ov_meth) = c("probe", ov_pats)
ov_meth = ov_meth[-1,]
saveRDS(ov_meth, file="OV_miRNA_Expression_data_Feb12.rds")

head -n 1 PAAD.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt > PAAD_patients.txt

#in R
paad_pats = read.table("PAAD_patients.txt")
paad_pats = paad_pats[,3:ncol(paad_pats)]
paad_pats = as.character(unlist(paad_pats))
paad_meth = fread("PAAD.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt")
colnames(paad_meth) = c("probe", paad_pats)
paad_meth = paad_meth[-1,]
saveRDS(paad_meth, file="PAAD_miRNA_Expression_data_Feb12.rds")

head -n 1 KIRC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt > KIRC_patients.txt

#in R
kirc_pats = read.table("KIRC_patients.txt")
kirc_pats = kirc_pats[,3:ncol(kirc_pats)]
kirc_pats = as.character(unlist(kirc_pats))
kirc_meth = fread("KIRC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt")
colnames(kirc_meth) = c("probe", kirc_pats)
kirc_meth = kirc_meth[-1,]
saveRDS(kirc_meth, file="KIRC_miRNA_Expression_data_Feb12.rds")

head -n 1 LIHC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt > LIHC_patients.txt

lihc_pats = read.table("LIHC_patients.txt")
lihc_pats = lihc_pats[,3:ncol(lihc_pats)]
lihc_pats = as.character(unlist(lihc_pats))
lihc_meth = fread("LIHC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt")
colnames(lihc_meth) = c("probe", lihc_pats)
lihc_meth = lihc_meth[-1,]
saveRDS(lihc_meth, file="LIHC_miRNA_Expression_data_Feb12.rds")

