###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)

rna = as.data.frame(rna)

###---------------------------------------------------------------

#function that tests each lncRNA's survival 

#1. remove discrepancy 
z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]

#2. list of cancers to apply function to 
cancers = as.list(unique(rna$Cancer))

#3. function that splits data into cancers 
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. function that calculates survival for each gene 

canc_survival_genes = function(dato){
	genes = as.list(colnames(dato)[2:(ncol(dato)-34)])
	#genes = genes[1:100]
	canc_data_genes_analyze = dato 
	
	get_survival = function(gene){
	print(gene)
  	results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
  	z = which(colnames(canc_data_genes_analyze) == gene)
  	dat = canc_data_genes_analyze[,c(1,z,(ncol(canc_data_genes_analyze)-33):ncol(canc_data_genes_analyze))]
  	dat$DFI.time = as.numeric(dat$DFI.time)
  	dat$DFI = as.numeric(dat$OS)
  	#remove NAs
  	z = which(is.na(dat$DFI.time))
  	if(!(length(z) ==0)){
  	dat = dat[-z,]}
	med_gene = median(as.numeric(dat[,2]))  	
	if(med_gene == 0){
		med_gene = mean(as.numeric(dat[,2]))
	}
	if(!(med_gene == 0)){
	dat$med = ""
	for(y in 1:nrow(dat)){
		check = dat[y,2] > med_gene
		if(check){
			dat$med[y] = 1
		}
		if(!(check)){
			dat$med[y] = 0
		}
	}
	res.cox <- coxph(Surv(DFI.time, DFI) ~ med, data = dat)
  	row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox)
  	return(row)
  	}
	}

	genes_survival = llply(genes, get_survival, .progress="text")
	genes_survival_res = ldply(genes_survival, rbind)
	#fdr
	colnames(genes_survival_res) = c("gene", "coef", "HR", "pval", "low95", "upper95")
	genes_survival_res$fdr = p.adjust(as.numeric(genes_survival_res$pval), method="fdr")
	genes_survival_res$canc = dato$Cancer[1]
	return(genes_survival_res)
}

all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)

###---------------------------------------------------------------

#plot scatter plot - HR versus p-value draw line for FDR = 0.05
all_cancers_genes_surv_comb$pval = -log10(as.numeric(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb$fdr = -log10(all_cancers_genes_surv_comb$fdr)
all_cancers_genes_surv_comb$HR = as.numeric(all_cancers_genes_surv_comb$HR)

lineval = -log10(0.05)

#facet by cancer type 
pdf("HR_vs_pval_disease_free_survival_all_cancers_volcano_plots.pdf", width=20, height=20)
ggscatter(all_cancers_genes_surv_comb, x = "HR", y = "pval", facet.by = "canc") +
geom_hline(yintercept = lineval) + geom_vline(xintercept = 1)
dev.off()































