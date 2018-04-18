###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

###Data
gtex = readRDS("all_lncRNAs_exp_scores_inGTEX_all_tissues_April17.rds")
tcga = readRDS("all_lncRNAs_exp_scores_inTCGA_all_tissues_April17.rds")

###Functions

#1. Divide data into matching tissues, one dataframe wtih both GTEx and TCGA 
#for one tissue

tcga_tis = unique(tcga$tis)

gtex$tis = str_sub(gtex$tis, 1, 4)
gtex_tis = unique(gtex$tis)

tis_match = as.data.frame(matrix(ncol=2)) ; colnames(tis_match) = c("cancer", "tis")

for(i in 1:length(gtex_tis)){
	 z = (which(str_detect(tcga_tis, gtex_tis[i])))
	 if(!(length(z)==0)){
	 	if(length(z)==1){
		 	canc = tcga_tis[z]
		 	tis = gtex_tis[i]
		 	row = c(canc, tis)
		 	names(row) = names(tis_match)
		 	tis_match = rbind(tis_match, row)
		 	}
	 	if(length(z)>1){
	 		for(k in 1:length(z)){
	 			canc = tcga_tis[z][k]
			 	tis = gtex_tis[i]
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
	all = rbind(dat_canc, dat_gt)
	return(all)
}

all_datas = llply(cancers, get_data)

#3. get fold change and p-value for each gene 

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


