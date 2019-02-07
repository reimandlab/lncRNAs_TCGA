setwd("/Users/kisaev/remote3/IC_exp_cands")

library(data.table)
library(dplyr)
library(stringr)

####################################################################################################
#LOAD DATA
####################################################################################################

#----RNA-Seq----------------------------------------------------------------------------------------

#SUBID results
subid = readRDS("subid_results_KI_IC_feb7.rds")
subid$method = "SUBID_erik_data"
colnames(subid)[c(1,6)] = paste(colnames(subid)[c(1,6)], "Eriks_data", sep="_")

#MEDIAN DICHOTOMIZED 
gbm_med = readRDS("GBM_median_splits_IonCHannels.rds")
colnames(gbm_med)[2:ncol(gbm_med)] = paste(colnames(gbm_med)[2:ncol(gbm_med)], "median", sep="_") 
colnames(gbm_med)[7] = "Cancer"
gbm_med = gbm_med[,c("gene", "HR_median", "pval_median", "fdr_pval_median", "name_median")]
colnames(gbm_med)[c(1,5)] = c("ensg", "gene")

#OUTLIER DICHOTOMIZED 
outliers = readRDS("TCGA_ION_CHANNEL_results_Jan3119.rds")
gbm_out = as.data.table(filter(outliers, cancer == "GBM"))
colnames(gbm_out)[2:ncol(gbm_out)] = paste(colnames(gbm_out)[2:ncol(gbm_out)], "outlier", sep="_") 
colnames(gbm_out)[7] = "Cancer"
gbm_out = gbm_out[,c("gene", "HR_outlier", "pval_outlier", "num_risk_outlier", "perc_risk_outlier", "fdr_pval_outlier")]
gbm_out$perc_risk_outlier = as.numeric(gbm_out$perc_risk_outlier)
colnames(gbm_out)[1] = c("ensg")

#----MICROARRAY-data--------------------------------------------------------------------------------
microarray = readRDS("GBM_ics_microarray_results_feb9.rds")

####################################################################################################
#SUMMARIZE DATA
####################################################################################################

#[1] Combine RNA-Seq data 
kidata = merge(gbm_med, gbm_out, by="ensg") #126 ion channels with both median and outlier available analyiss 
kidata_sub = merge(kidata, subid, by="gene") #123 ion channels with both median, outlier and subid available 
kidata_sub$HR_median = as.numeric(kidata_sub$HR_median)
kidata_sub$HR_outlier = as.numeric(kidata_sub$HR_outlier)
kidata_sub$pval_outlier = as.numeric(kidata_sub$pval_outlier)

#[2] Combine microarray or add tag because most ICs may not have affy data
z = which(kidata_sub$gene %in% microarray$gene)
noaffy = kidata_sub[-z,]
noaffy$HR_microarray = ""
noaffy$pval_microarray = ""
noaffy$conc_microarray = ""
noaffy$fdr_microarray = ""
alldat = merge(kidata_sub, microarray, by="gene")
alldat = rbind(alldat, noaffy)

alldat$HR_microarray = as.numeric(alldat$HR_microarray)
alldat$pval_microarray = as.numeric(alldat$pval_microarray)
alldat$conc_microarray = as.numeric(alldat$conc_microarray)
alldat$fdr_microarray = as.numeric(alldat$fdr_microarray)

#round all digits
z = which(colnames(alldat) %in% c("HR_median", "pval_median", "fdr_pval_median", "HR_outlier", "pval_outlier", "perc_risk_outlier", "fdr_pval_outlier", "median_p_Eriks_data",
          "subid_p","subid_fdr",  "median_fdr_Eriks_data", "HR_microarray", "pval_microarray", "conc_microarray", "fdr_microarray"))
alldat = as.data.frame(alldat)
alldat[,z] = apply(alldat[,z], 2, function(x){round(x, digits=4)})
alldat = as.data.table(alldat)
alldat$method = NULL

saveRDS(alldat, file="ion_channels_4_methods_survival_analysis_KI_Feb7.rds")
write.csv(alldat, file="ion_channels_4_methods_survival_analysis_KI_Feb7.csv", quote=F, row.names=F)

alldat = as.data.table(alldat)
alldat = alldat[order(pval_outlier, subid_p, pval_median)]
alldat$rank = 1:nrow(alldat)

cands = c("GJB2", "CATSPER1", "AQP9", "SCN9A")
cands_dat = filter(alldat, gene %in% cands)
write.csv(cands_dat, file="ion_channels_four_exo_cands_KI_Feb7.csv", quote=F, row.names=F)


####################################################################################################
#MAKE FIGURE FOR REPORT 
####################################################################################################







