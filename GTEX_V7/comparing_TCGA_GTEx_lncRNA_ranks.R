###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

###Data
gtex = readRDS("all_lncRNAs_exp_scores_inGTEX_all_tissues_May3.rds")
tcga = readRDS("TCGA_all_TCGA_cancers_scored_byindexMay4.rds")

#summary of lncRNAs detected in each cancer 
lncs_det = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
lncs_det_info = readRDS("summary_detectable_lncs_howmanycancers_typesLNCRNAS.rds")

#********get distribution of lncRNA rank ? ****************************************

###Functions

#1. Divide data into matching tissues, one dataframe wtih both GTEx and TCGA 
#for one tissue

tcga_canc = unique(tcga$canc)

gtex$canc = str_sub(gtex$canc, 1, 4)
gtex_canc = unique(gtex$canc)

tis_match = as.data.frame(matrix(ncol=2)) ; colnames(tis_match) = c("cancer", "tis")

for(i in 1:length(gtex_canc)){
	 z = (which(str_detect(tcga_canc, gtex_canc[i])))
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
	z = which(tcga$canc == cancer)
	dat_canc = tcga[z,]
	t = tis_match$tis[which(tis_match$cancer ==cancer)]
	z = which(gtex$canc == t)
	dat_gt = gtex[z,]
	all = rbind(dat_canc, dat_gt)
	return(all)
}

all_datas = llply(cancers, get_data)

#3. get median rank for each gene within tumour and gtex 
allCands <- readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")

#plot resepctive candidates wtihin this cancer type 


#4. get fold change and p-value for each gene 

get_fc = function(dataframee){
	genes = as.list(unique(dataframee$gene))
	gene_analyze = function(gene){
		x = dataframee[dataframee$data=="TCGA", ]
		x = x$score[x$gene == gene]
		y = dataframee[dataframee$data=="GTEX", ]
		y = y$score[y$gene == gene]
		t = t.test(x, y, alternative = "two.sided")
		p = t$p.value
		fc = t$estimate[1]/t$estimate[2]
		row = c(gene, fc, p) #cancer/gtex
		names(row) = c("gene", "fc", "pval")
		return(row)
	}
	all_genes_results = llply(genes, gene_analyze, .progress="text")
	all_genes_results2 = ldply(all_genes_results, rbind)
	all_genes_results2$fc = log2(as.numeric(all_genes_results2$fc))
	all_genes_results2$fdr = p.adjust(all_genes_results2$pval, method="fdr")
	all_genes_results2$fdrtag = ""
	all_genes_results2$fdrtag[all_genes_results2$fdr<=0.05] = "FDRsig"
	all_genes_results2$fdrtag[all_genes_results2$fdr>0.05] = "NotFDRsig"
	all_genes_results2$pval = -log10(as.numeric(all_genes_results2$pval))
	ggscatter(all_genes_results2, x = "fc", y = "pval", color = "fdrtag", size=0.25, title=dataframee$tis[1])
}

pdf("volcano_plots_foldchange_of_lncRNA_ranks_TCGA_vs_GTEX.pdf")
llply(all_datas, get_fc, .progress="text")
dev.off()


