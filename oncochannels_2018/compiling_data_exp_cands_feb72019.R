setwd("/Users/kisaev/remote/IC_exp_cands")

library(data.table)
library(dplyr)
library(stringr)
library(metap)

#merging p-value function
merge_p_values <- function(scores, method=c("Fisher", "Brown", "logitp",
                                            "meanp", "sump", "sumz", "sumlog")) {
  # Validation on scores
  if (is.list(scores)) scores <- unlist(scores, recursive=FALSE)
  if (!(is.vector(scores) || is.matrix(scores))) stop("scores must be a matrix or list")
  if (any(is.na(scores))) stop("scores may not contain missing values")
  if (!is.numeric(scores)) stop("scores must be numeric")
  if (any(scores < 0 | scores > 1)) stop("All values in scores must be in [0,1]")
  
  method <- match.arg(method)
  if(method == "Fisher") method <- "sumlog"
  
  if (is.vector(scores)) {
    if (method == "Brown") stop("Brown's method cannot be used with a single list of p-values")
    
    # Some metap functions don't like p-values that are 0 or 1 so make them (0, 1) to avoid errors
    scores <- sapply(scores, function(x) if (x == 0) 1e-16 else if (x==1) 1-1e-16 else x)
    func <- function(x) getFromNamespace(method, 'metap')(x)$p
    return(func(scores))
  }
  
  # scores is a matrix
  if (ncol(scores) == 1) return (scores[, 1, drop=TRUE])
  
  if (method == "Brown") {
    cov.matrix <- calculateCovariances(t(scores))
    return(apply(scores, 1, brownsMethod, cov.matrix=cov.matrix))
  }
  
  scores <- apply(scores, c(1,2), function(x) if (x == 0) 1e-16 else if (x==1) 1-1e-16 else x)
  func <- function(x) getFromNamespace(method, 'metap')(x)$p
  return (apply(scores, 1, func))
}

brownsMethod <- function(p.values, data.matrix=NULL, cov.matrix=NULL) {
  if (missing(data.matrix) && missing(cov.matrix)) {
    stop ("Either data.matrix or cov.matrix must be supplied")
  }
  if (!(missing(data.matrix) || missing(cov.matrix))) {
    message("Both data.matrix and cov.matrix were supplied. Ignoring data.matrix")
  }
  if (missing(cov.matrix)) cov.matrix <- calculateCovariances(data.matrix)
  
  N <- ncol(cov.matrix)
  expected <- 2 * N
  cov.sum <- 2 * sum(cov.matrix[lower.tri(cov.matrix, diag=FALSE)])
  var <- (4 * N) + cov.sum
  sf <- var / (2 * expected)
  
  df <- (2 * expected^2) / var
  if (df > 2 * N) {
    df <- 2 * N
    sf <- 1
  }
  
  x <- 2 * sum(-log(p.values), na.rm=TRUE)
  p.brown <- pchisq(df=df, q=x/sf, lower.tail=FALSE)
  p.brown
}


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
#alldat[,z] = apply(alldat[,z], 2, function(x){round(x, digits=4)})
alldat = as.data.table(alldat)
alldat$method = NULL

alldat$HR_match = ""
z = which((alldat$HR_median >1) & (alldat$HR_outlier >1))
alldat$HR_match[z] = "yes"
alldat = alldat[z,]

saveRDS(alldat, file="ion_channels_4_methods_survival_analysis_KI_Feb22.rds")
write.csv(alldat, file="ion_channels_4_methods_survival_analysis_KI_Feb22.csv", quote=F, row.names=F)

alldat = as.data.table(alldat)
alldat = alldat[order(pval_outlier, subid_p, pval_median)]
alldat$rank = 1:nrow(alldat)

cands = c("GJB2", "CATSPER1", "AQP9", "SCN9A")
cands_dat = filter(alldat, gene %in% cands)
write.csv(cands_dat, file="ion_channels_four_exo_cands_KI_Feb22.csv", quote=F, row.names=F)

####################################################################################################
#MAKE FIGURE FOR REPORT 
####################################################################################################

#keep just the cols with the pvalues
ps = c("gene", "pval_median", "pval_outlier", "subid_p", "pval_microarray")
alldat = alldat[,..ps]
z = which(is.na(alldat$pval_microarray))
alldat$pval_microarray[z] = 1

alldat = as.data.frame(alldat)
rownames(alldat) = alldat$gene
alldat$gene = NULL

#get fisher's merged p-value
#alldat$fish = apply(alldat, 1, function(x){merge_p_values(x)})
#alldat$fish = round(alldat$fish, digits=7)

alldat = as.matrix(alldat)
browns = merge_p_values(alldat, method="Brown")
alldat = as.data.frame(alldat)
alldat$browns = browns

#reorder
alldat$gene = rownames(alldat)

cols =  colnames(alldat)[1:5]
alldat = alldat[,c("gene", cols)]

alldat = as.data.table(alldat)
alldat = alldat[order(browns)]
z = which(alldat$gene %in% cands)
alldat$cands = ""
alldat$cands[z] = "yes"

write.csv(alldat, file="ion_channels_merged_pvalues_browns_KI_onlyhazardours_220219.csv", quote=F, row.names=F)

#add fdr
alldat$fdr = p.adjust(alldat$browns, method="fdr")

write.csv(alldat, file="ion_channels_merged_pvalues_browns_KI_onlyhazardours_withFDR_010319.csv", quote=F, row.names=F)

