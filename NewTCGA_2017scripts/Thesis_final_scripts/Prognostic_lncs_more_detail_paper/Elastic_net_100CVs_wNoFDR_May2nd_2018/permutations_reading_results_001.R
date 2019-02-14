
#ran elastic net simulations for each cancer types 100 times 
#this means 100 runs of 100 elastic net cross-validations 
#should be 100 files per cancer types if all ran correctly 

#results files are here: 
#/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/elastic_net_permutations_results

library(data.table)
library(dplyr)
library(stringr)
library(plyr)
library(ggplot2)
library(ggpubr)

#get cancer names 

all_files = list.files()

f = list.files()[[1]]
unlist(strsplit(f, "_"))[1]

cancers = unique(sapply(all_files, function(x){unlist(strsplit(x, "_"))[1]}))


#[[1]]--------------------------------------------------------------------------
#seperate files into cancer types

get_canc_files = function(canc){
	z = which(str_detect(all_files, canc))
	canc_f = all_files[z]
	return(canc_f)
}


cancer_files = llply(cancers, get_canc_files)


#[[2]]--------------------------------------------------------------------------
#get number of genes in each fold for each cancer type 

get_canc_lncs = function(canc_files){

	cancer = unlist(strsplit(canc_files[[1]], "_"))[1]
	print(cancer)
	z = which(str_detect(canc_files, "lncRNAs_selected_permutations_1000CV_1000_no_fdr_ELASTICNET.rds"))
	lncs = canc_files[z]
	#read all files and summarize how many lncRNAs were selected more than 50 times in each file 
	#and make barplot and save in data file

	read_f = function(filee){
		f = readRDS(filee)
		genes_list = as.data.table(table(unlist(f)))
		genes_list = genes_list[order(N)]
		genes_list$round = unlist(strsplit(filee, "_"))[2]
		#genes_list = as.data.table(filter(genes_list, N >=50))
		return(genes_list)
	}

	lncs_sum = llply(lncs, read_f)
	lncs_sum = as.data.table(ldply(lncs_sum))
	lncs_sum$round = as.numeric(lncs_sum$round)
	lncs_sum = lncs_sum[order(round)]

	if(!(dim(lncs_sum)[1]==0)){

	t = as.data.table(table(lncs_sum$round))

	g = ggbarplot(t, x = "V1", y = "N",
          ylab = "# lncRNAs selected more than 50%",
          xlab = "Round, n=100") + theme_bw()
	g = ggpar(g, font.tickslab = c(6,"plain", "black"),
 	xtickslab.rt = 90, yticks.by = 1, title = cancer)
	print(g)

	lncs_sum$cancer = cancer
	return(lncs_sum)

}
}

pdf("summary_barplots_all_cancers_feb13.pdf", width=10, height=5)
results = llply(cancer_files, get_canc_lncs, .progress="text")
dev.off()

results1 = as.data.table(ldply(results))
colnames(results1)[1] = "lncRNA"