###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

###Data
gtex = readRDS("allGTEX_lncRNAs_scored_Feb2619.rds")
tcga = readRDS("TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

#cands -- lncRNAs 

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

#summary of lncRNAs detected in each cancer 
#lncs_det = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
#lncs_det_info = readRDS("summary_detectable_lncs_howmanycancers_typesLNCRNAS.rds")

#********get distribution of lncRNA rank ? ****************************************

###Functions

#1. Divide data into matching tissues, one dataframe wtih both GTEx and TCGA 
#for one tissue

tcga_canc = unique(tcga$tis)

gtex$simp_tis = gtex$tis
gtex$simp_tis = str_sub(gtex$simp_tis, 1, 4)
gtex_canc = unique(gtex$tis)

tis_match = as.data.frame(matrix(ncol=2)) ; 
colnames(tis_match) = c("cancer", "tis")

for(i in 1:length(gtex_canc)){
	 z = (which(str_detect(tcga_canc, gtex$simp_tis[which(gtex$tis == gtex_canc[[i]])][1])))
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

#add glioblastoma manually
gbm = subset(tis_match, cancer == "Brain Lower Grade Glioma")
gbm$cancer = "Glioblastoma multiforme"
tis_match = rbind(tis_match, gbm) #54 comparisons in total 


#2. type of cancers with tissues available -> seperate into dataframes
cancers = as.data.frame(unique(tis_match$cancer))
colnames(cancers)[1] = "canc"
canc_conv = readRDS("cancers_conv_july23.rds")
cancers = merge(cancers, canc_conv, by="canc") #23 of 32 cancer types had gtex comparison tissue 

#remove cancer types with less than 50 patients
z = which(cancers$type %in% c("KICH", "CHOL", "DLBC", "UCS"))
cancers = cancers[-z,]
cancers = cancers$canc
cancers = tis_match$tis
tis_match$combo = paste(tis_match$cancer, tis_match$tis, sep = "_")
cancers = tis_match$combo

get_data = function(cancer){
	
	print(cancer)
	tis = unlist(strsplit(cancer, "_"))[2]
	canc = unlist(strsplit(cancer, "_"))[1]

	z = which(gtex$tis == tis)
	dat_gt = gtex[z,]

	z = which(tcga$tis == canc)
	dat_canc = tcga[z,]
	dat_canc$simp_tis = dat_canc$tis
	dat_canc = dat_canc[,c("exp", "gene", "score", "patient", "data", "tis", "simp_tis")]

	all = rbind(dat_canc, dat_gt) #fix order of columns 
	return(all)
}

all_datas = llply(cancers, get_data, .progress="text")

#3. get median rank for each gene within tumour and gtex 

#plot resepctive candidates wtihin this cancer type 
lncs_gtex = unique(gtex$gene)
lncs_tcga = unique(tcga$gene)
lncs_both = lncs_tcga[which(lncs_tcga %in% lncs_gtex)] #271 ion channels 

library(plyr)
library(doParallel)

registerDoParallel(cores=2)

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
		canc = dataframee$tis[dataframee$data=="TCGA"][1]
		row = c(gene, fc, p, med_diff, canc, dataframee$tis[dataframee$data=="GTEX"][1]) #cancer/gtex
		names(row) = c("gene", "fc_mean", "pval_wilcoxon", "median_difference", "canc", "tis")
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
	#all_genes_results2 = all_genes_results2[which(abs(as.numeric(all_genes_results2$median_difference)) >=0.1),]
	#all_genes_results2$pval_wilcoxon = -log10(as.numeric(all_genes_results2$pval_wilcoxon))
	#print(ggdensity(dataframee, x = "score",
  	#		 color = "data", palette = c("#00AFBB", "#E7B800"), title=dataframee$tis[1]))

	return(all_genes_results2)
}

#pdf("volcano_plots_foldchange_of_lncRNA_ranks_TCGA_vs_GTEX_May24.pdf")
results_analysis = llply(all_datas, get_fc, .parallel=TRUE)
#dev.off()
#saveRDS(results_analysis, file="results_analysis_July24.rds")
#saveRDS(results_analysis, file="results_analysis_Dec30.rds")
saveRDS(results_analysis, file="results_analysis_Feb26.rds")


