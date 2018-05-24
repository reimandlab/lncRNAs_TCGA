###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

###Data
gtex = readRDS("allGTEX_lncRNAs_scored_May23.rds")
tcga = readRDS("TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")

#summary of lncRNAs detected in each cancer 
#lncs_det = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
#lncs_det_info = readRDS("summary_detectable_lncs_howmanycancers_typesLNCRNAS.rds")

#********get distribution of lncRNA rank ? ****************************************

###Functions

#1. Divide data into matching tissues, one dataframe wtih both GTEx and TCGA 
#for one tissue

tcga_canc = unique(tcga$tis)

gtex$tis = str_sub(gtex$tis, 1, 4)
gtex_canc = unique(gtex$tis)

tis_match = as.data.frame(matrix(ncol=2)) ; 
colnames(tis_match) = c("cancer", "tis")

for(i in 1:length(gtex_canc)){
	 z = (which(str_detect(tcga_canc, gtex_canc[[i]])))
	 if(!(length(z)==0)){
	 	if(length(z)==1){
		 	canc = tcga_canc[z]
		 	tis = gtex_canc[i]
		 	row = c(canc, tis)
		 	names(row) = names(tis_match)
		 	tis_match = rbind(tis_match, row)
		 	}
	 	if(length(z)>1){
	 		for(k in 1:length(z)){
	 			canc = tcga_canc[z][k]
			 	tis = gtex_canc[i]
			 	row = c(canc, tis)
			 	names(row) = names(tis_match)
			 	tis_match = rbind(tis_match, row)
		 	}
	 	}
	 }

}

tis_match = tis_match[-1,]

#2. type of cancers with tissues available -> seperate into dataframes
cancers = unique(tis_match$cancer)

get_data = function(cancer){
	z = which(tcga$tis == cancer)
	dat_canc = tcga[z,]
	t = tis_match$tis[which(tis_match$cancer ==cancer)]
	z = which(gtex$tis == t)
	dat_gt = gtex[z,]
	all = rbind(dat_canc, dat_gt) #fix order of columns 
	return(all)
}

all_datas = llply(cancers, get_data, .progress="text")

#3. get median rank for each gene within tumour and gtex 
allCands <- readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")
allCands = filter(allCands, AnalysisType == "noFDR")

#plot resepctive candidates wtihin this cancer type 

lncs_gtex = unique(gtex$gene)
lncs_tcga = unique(tcga$gene)
lncs_both = lncs_tcga[which(lncs_tcga %in% lncs_gtex)]

#4. get fold change and p-value for each gene 
get_fc = function(dataframee){
	#which lncRNAs actually covered by both datasets 
	genes = as.list(lncs_both)
	gene_analyze = function(gene){
		x = dataframee[dataframee$data=="TCGA", ]
		x = x$score[x$gene == gene]
		y = dataframee[dataframee$data=="GTEX", ]
		y = y$score[y$gene == gene]
		t = wilcox.test(x, y, alternative = "two.sided")
		p = t$p.value
		fc = mean(x)/mean(y)
		med_diff = median(x)-median(y)
		row = c(gene, fc, p, med_diff) #cancer/gtex
		names(row) = c("gene", "fc_mean", "pval_wilcoxon", "median_difference")
		return(row)
	}
	#genes=genes[1:10]
	all_genes_results = llply(genes, gene_analyze, .progress="text")
	all_genes_results2 = ldply(all_genes_results, rbind)
	all_genes_results2$fc_mean = as.numeric(all_genes_results2$fc_mean)
	all_genes_results2$pval_wilcoxon = as.numeric(all_genes_results2$pval_wilcoxon)
	all_genes_results2$fc_mean = log2(all_genes_results2$fc_mean)
	all_genes_results2$fdr = p.adjust(all_genes_results2$pval_wilcoxon, method="fdr")
	all_genes_results2$fdrtag = ""
	all_genes_results2$fdrtag[all_genes_results2$fdr<=0.05] = "FDRsig"
	all_genes_results2$fdrtag[all_genes_results2$fdr>0.05] = "NotFDRsig"
	#subest to only those with 25% difference in corresponding median gene ranks
	all_genes_results2 = all_genes_results2[which(abs(as.numeric(all_genes_results2$median_difference)) >=0.1),]
	all_genes_results2$pval_wilcoxon = -log10(as.numeric(all_genes_results2$pval_wilcoxon))
	print(ggdensity(dataframee, x = "score",
  			 color = "data", palette = c("#00AFBB", "#E7B800"), title=dataframee$tis[1]))

	return(all_genes_results2)
}

pdf("volcano_plots_foldchange_of_lncRNA_ranks_TCGA_vs_GTEX_May24.pdf")
results_analysis = llply(all_datas, get_fc, .progress="text")
dev.off()



