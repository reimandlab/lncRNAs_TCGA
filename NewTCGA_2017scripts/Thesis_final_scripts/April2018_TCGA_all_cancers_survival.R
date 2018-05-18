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
dim(rna)
dim(pcg)
dim(norm)
dim(met)

all_genes = as.data.frame(unique(c(colnames(rna), colnames(pcg))))
z = which(str_detect(all_genes[,1], "ENSG"))
all_genes = all_genes[z,]
all_genes = as.data.frame(all_genes)
colnames(all_genes)[1] = "gene"

all_genes$type = ""
z = which(all_genes$gene %in% colnames(rna))
all_genes$type[z] = "lncRNA"
z = which(all_genes$gene %in% colnames(pcg))
all_genes$type[z] = "pcg"

saveRDS(all_genes, file="all_genes_used_in_TCGA_april17.rds")

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

det_lncs = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
det_lncs =filter(det_lncs, status =="detectable")

canc_survival_genes = function(dato){
	#genes = as.list(colnames(dato)[2:(ncol(dato)-34)])
	genes = unique(det_lncs$lncRNA[which(det_lncs$cancer %in% dato$Cancer[1])])	

	#genes = genes[1:100]
	canc_data_genes_analyze = dato 
	
	get_survival = function(gene){
	print(gene)
  	results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
  	z = which(colnames(canc_data_genes_analyze) == gene)
  	dat = canc_data_genes_analyze[,c(1,z,(ncol(canc_data_genes_analyze)-33):ncol(canc_data_genes_analyze))]
  	dat$OS.time = as.numeric(dat$OS.time)
  	dat$OS = as.numeric(dat$OS)
  	#remove NAs
  	z = which(is.na(dat$OS.time))
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
	res.cox <- coxph(Surv(OS.time, OS) ~ med, data = dat)
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
	genes_survival_res = as.data.table(genes_survival_res)
	genes_survival_res = genes_survival_res[order(fdr)]
	return(genes_survival_res)
}

all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)


saveRDS(all_cancers_genes_surv_comb, file="lncRNAs_for_plotting_HAzard_Ratios_Pvalues_May18.rds")


###---------------------------------------------------------------

#plot scatter plot - HR versus p-value draw line for FDR = 0.05
all_cancers_genes_surv_comb$pval = -log10(as.numeric(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb$fdr = -log10(all_cancers_genes_surv_comb$fdr)
all_cancers_genes_surv_comb$HR = as.numeric(all_cancers_genes_surv_comb$HR)

z = which(is.na(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

z1 = which(all_cancers_genes_surv_comb$fdr == "Inf")
z2 = which(all_cancers_genes_surv_comb$upper95 == "Inf")
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-c(z1,z2),]

z = which(all_cancers_genes_surv_comb$HR > 5)
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

lineval = -log10(0.05)

all_cancers_genes_surv_comb$fdrsig = ""
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr >= lineval] = "FDRsig"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr < lineval] = "FDRnotSig"
z = which((all_cancers_genes_surv_comb$fdr < lineval) & (all_cancers_genes_surv_comb$fdr > -log10(0.1)))
all_cancers_genes_surv_comb$fdrsig[z] = "FDR00.1"

#facet by cancer type 
pdf("HR_vs_pval_survival_all_cancers_volcano_plots.pdf", width=12, height=12)
ggscatter(all_cancers_genes_surv_comb, x = "HR", y = "pval", color = "fdrsig", facet.by = "canc", palette=mypal[c(2,1,3)], ylab="-log10(p-value)") +
geom_hline(yintercept=lineval, linetype="dashed", color = "red") + geom_vline(xintercept = 1, linetype="dashed", color = "blue")
dev.off()

saveRDS(all_cancers_genes_surv_comb, file="all_cancers_all_genes_univariate_survival_results_April16.rds")





























