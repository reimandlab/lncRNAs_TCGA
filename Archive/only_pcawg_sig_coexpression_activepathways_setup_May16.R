#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 28th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#List of canddidates and cox results
#allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
#List of canddidates and cox results
#allCands <- readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$pval = as.numeric(allCands$pval)

#subset to those that validated in PCAWG
allCands = as.data.table(allCands)
allCands = filter(allCands, pval <= 0.15)
allCands$Combo = NULL
#check if their HR matches HR in TCGA
allCands$Cancer = NULL
allCands$freq = NULL

allCands = allCands[!duplicated(allCands), ]

lncs = as.list(as.character(unique(allCands$gene[allCands$data == "PCAWG"])))

check_match = function(lnc){
	z = which(allCands$gene == lnc)
	res = as.data.table(allCands[z,])
	
	if(dim(res)[1] > 2){
	
	canc = res$cancer[which(duplicated(res$cancer))]
	res = filter(res, cancer == canc)
	}

	test = which(as.numeric(res$HR) >= 1)
	test = length(test)
	if(test ==2){
		match = "match"
	}
	if(test ==0){
		match = "match"
	}
	if(test ==1){
		match = "nomatch"
	}
	canc = unique(res$canc)
	return(c(lnc, match, canc))
}

matches = llply(lncs, check_match)
matches = as.data.table(as.data.frame(do.call("rbind", matches)))
matches = filter(matches, V2 == "match")
colnames(matches) = c("lnc", "match", "cancer")

allCands = as.data.table(allCands)
allCands = filter(allCands, gene %in% matches$lnc, cancer %in% matches$cancer) #######----only PCAWG singificant 

#6 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
z = which(str_detect(colnames(all), "ENSG"))	
all[,z] <- log1p(all[,z])

#2. Get lncRNA - median within each tissue type
tissues <- unique(allCands$cancer)

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- all[all$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#Function 2
#for each cancer 

get_lnc_canc = function(dat){
	canc = dat$Cancer[1]
	lncs = as.character(unique(subset(allCands, cancer == canc)$gene))
	
	evaluate_each_lnc = function(lnc){
		pcgs = colnames(pcg)[2:19351]
		dat_keep = dat[,which(colnames(dat) %in% c("patient", lnc, pcgs))]
		rownames(dat_keep) = dat_keep$patient
		dat_keep$patient = NULL
		#figure out which patients are high risk and which patients low risk
		dat_keep$median <- ""
 		median2 <- quantile(as.numeric(dat_keep[,1]), 0.5)

	  	 if(median2 ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(dat_keep[,1] > 0)
		    l2 = which(dat_keep[,1] ==0)
		    dat_keep$median[l1] = 1
		    dat_keep$median[l2] = 0
		    }

		  if(!(median2 ==0)){
		    l1 = which(dat_keep[,1] >= median2)
		    l2 = which(dat_keep[,1] < median2)
		    dat_keep$median[l1] = 1
		    dat_keep$median[l2] = 0
		}

		#which one is high risk --> need surivval data
		dat_keep$patient = rownames(dat_keep)
		
  		dat_keep$median[dat_keep$median ==0] = "Low"
  		dat_keep$median[dat_keep$median==1] = "High"

  		#cox ph
  		z = which((allCands$gene == lnc) & (allCands$cancer == canc))

  		HR = as.numeric(allCands$HR[z])
  		
  		if(HR <1){
  			risk = "Low"
  			dat_keep$risk = ""
  			dat_keep$risk[dat_keep$median=="High"] ="no"
  			dat_keep$risk[dat_keep$median=="Low"] ="yes"
  		}
  		if(HR >1){
  			risk = "High"
  			dat_keep$risk = ""
  			dat_keep$risk[dat_keep$median=="High"] ="yes"
  			dat_keep$risk[dat_keep$median=="Low"] ="no"
  		}

  		dat_keep$lnc = colnames(dat_keep)[1]
  		dat_keep$canc = canc
  		colnames(dat_keep)[1] = "lncRNA"

  	return(dat_keep)	
	}#end function evaluate_each_lnc

	results_lncs = llply(lncs, evaluate_each_lnc, .progress="text")
	results_lncs1 = as.data.frame(do.call("rbind", results_lncs))
}

all_canc_lnc_data = llply(tissues_data, get_lnc_canc, .progress="text")


#Function 3
#wtihin each cancer 
#calculate for each lncRNAs which PCGs have enriched expression in the high risk group

get_pcg_enrichment = function(dat){
	lncs = unique(dat$lnc)
	
	get_pcgs_high = function(lncrna){
		newdat = subset(dat, lnc == lncrna)
		#which PCGs have higher expression in the high risk group
		pcgs = colnames(newdat)[which(str_detect(colnames(newdat), "ENSG"))]
		pcgs_exps = newdat[,c(pcgs)]
		#medians = apply(pcgs_exps, 2, median)
		#pcgs = names(medians[medians > 2])
		print(length(pcgs))
		#lnc_pcg_results = as.data.frame(matrix(ncol=5)) ; colnames(lnc_pcg_results) = c("lnc", "pcg", "canc", "mean_diff", "pvalue")
		#pcgs = pcgs[1:10]
		get_correlation = function(pcg){
			p = pcg
			z = which(colnames(newdat) %in% p)
			lncpcg = newdat[,c(z, 1:2, 19353:ncol(newdat))]
			colnames(lncpcg)[1] = "pcgExp"
			fit <- lm(pcgExp ~ risk, data=lncpcg)
			#get p-value and generate boxplot with wilcoxon p-value 
			fit_pval = summary(fit)$coefficients[2,4]
			#which group is it higher in? 
			mean_diff = mean(lncpcg$pcgExp[lncpcg$risk == "yes"])/mean(lncpcg$pcgExp[lncpcg$risk == "no"])
			#if higher than 1 --> more expressed in risk group, less than 1 --> more expressed in low risk group
			#g = ggboxplot(lncpcg, x = "risk", y="pcgExp", color="median", title=paste(lncpcg$lnc[1], p, lncpcg$canc[1]))
			#g = g + stat_compare_means()
			#print(g)
			row = c(lncpcg$lnc[1], p, lncpcg$canc[1], mean_diff, fit_pval)
			return(row)
			#names(row) = colnames(lnc_pcg_results)
			#lnc_pcg_results = rbind(lnc_pcg_results, row)
		
		}#end get_correlation function 
		
		pcgs_results = llply(pcgs, get_correlation, .progress="text")
		pcgs_results1 = as.data.frame(do.call("rbind", pcgs_results))
		colnames(pcgs_results1) = c("lnc", "pcg", "canc", "mean_diff", "pvalue")
		return(pcgs_results1)

	}  #end get_pcgs_high function
	
	results_lncs = llply(lncs, get_pcgs_high, .progress="text")
	results_lncs1 = as.data.frame(do.call("rbind", results_lncs))
	#results for all lncRNA-PCG correlations in a single cancer type 
	return(results_lncs1)

}

#all_canc_lnc_data = all_canc_lnc_data[1:2] ###TEST CASE -------------------------------------------------------------
all_canc_lnc_data = llply(all_canc_lnc_data, get_pcg_enrichment, .progress="text")


all_canc_lnc_data1 = as.data.frame(do.call("rbind", all_canc_lnc_data))


saveRDS(all_canc_lnc_data1, file="7lncRNAs_pcawg_val_coexpression_analysis_june11th_allCands.rds")


#divide into high risk and low risk lncRNAs
high_risk = subset(all_canc_lnc_data1, mean_diff >=1.2) #should set higher mean difference threshold?
low_risk = subset(all_canc_lnc_data1, mean_diff <=0.8) #should set higher mean difference threshold? 

library(reshape2)

#pcgs enriched in high risk lncRNAs 
high_risk_matrix = acast(high_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
      value.var = 'pvalue', fill = 'na')

#pcgs enriched in low risk lncRNAS
low_risk_matrix = acast(low_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
      value.var = 'pvalue', fill = 'na')

#columns are lncRNAs and rows are PCGs

saveRDS(high_risk_matrix, file="high_risk_matrix_lncRNA_candidates_June11_only7cands.rds")
saveRDS(low_risk_matrix, file="low_risk_matrix_lncRNA_candidates_June11_only7cands.rds")



#--------------------------------CANCER SPECIFIC------------------------------------------------------------------------


#Seperate all_canc_lnc_data1 by cancer type so that can look at PCGs
#enriched in high and low risk groups by cancer type 
#most cancer types have multiple candidates 

out <- split(all_canc_lnc_data1 , forcats::fct_inorder(factor(all_canc_lnc_data1$canc)))

convert_high_risk_matrix = function(dtt){

	high_risk = subset(dtt, mean_diff >=1.5) #should set higher mean difference threshold?
	#remove means that are NaN or Inf
	z = which(high_risk$mean_diff ==NaN)
	if(!(length(z)==0)){
	high_risk = high_risk[-z,]}

	z = which(high_risk$mean_diff ==Inf)
	if(!(length(z)==0)){
	high_risk = high_risk[-z,]}

	library(reshape2)

	#pcgs enriched in high risk lncRNAs 
	high_risk_matrix = acast(high_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
    	  value.var = 'pvalue', fill = 'na')

	return(high_risk_matrix)
	
	}

order_cancers = unique(all_canc_lnc_data1$canc)

high_risk_pcgs_by_cancer = llply(out, convert_high_risk_matrix, .progress="text")

#by cancer type --> PCGs with higher expression in risk group 
saveRDS(high_risk_pcgs_by_cancer, file= "high_risk_matrix_lncRNAs_by_cancer_type_June7th.rds")


convert_low_risk_matrix = function(dtt){

	low_risk = subset(dtt, mean_diff <= 0.75) #should set higher mean difference threshold?
	#remove means that are NaN or Inf
	z = which(low_risk$mean_diff ==NaN)
	if(!(length(z)==0)){
	low_risk = low_risk[-z,]}

	z = which(low_risk$mean_diff ==Inf)
	if(!(length(z)==0)){
	low_risk = low_risk[-z,]}

	library(reshape2)

	#pcgs enriched in high risk lncRNAs 
	low_risk_matrix = acast(low_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
    	  value.var = 'pvalue', fill = 'na')

	return(low_risk_matrix)
}

low_risk_pcgs_by_cancer = llply(out, convert_low_risk_matrix, .progress="text")

#by cancer type --> PCGs with higher expression in lower risk group 
saveRDS(low_risk_pcgs_by_cancer, file= "low_risk_matrix_lncRNAs_by_cancer_type_June7th.rds")

saveRDS(order_cancers, file="order_cancers_risk_PCGs_matrices_June7th.rds")








