set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files

#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_name")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

cancers = c("Brain Lower Grade Glioma", "Rectum adenocarcinoma", "Breast invasive carcinoma")

get_plot = function(canc){

	print(canc)

	lncs = as.data.table(filter(allCands, cancer==canc))$gene
	cancer = canc_conv$type[canc_conv$Cancer==canc]	

	print(lncs)

	#pdf(paste("/u/kisaev/", cancer, "PFI_OS_curves.pdf", sep="_"))

	splots <- list()

	for(i in 1:length(lncs)){
		print(i)
		if(i==1){
		splots[[i]] <- get_km_plot_os(lncs[i], cancer)
		splots[[i+1]] <- get_km_plot_pfi(lncs[i], cancer)
		}
		if(!(i==1)){
		splots[[i+1]] <- get_km_plot_os(lncs[i], cancer)
		splots[[i+2]] <- get_km_plot_pfi(lncs[i], cancer)
		}
	}

	#dev.off()
	# Arrange multiple ggsurvplots and print the output
	res = arrange_ggsurvplots(splots, print = TRUE,
  	ncol = 2, nrow = length(lncs))
  	ggsave(paste("/u/kisaev/", cancer, "PFI_OS_curves.pdf", sep="_"), res, width=12, height=12)

}


get_plot(cancers[1])
get_plot(cancers[2])
get_plot(cancers[3])
