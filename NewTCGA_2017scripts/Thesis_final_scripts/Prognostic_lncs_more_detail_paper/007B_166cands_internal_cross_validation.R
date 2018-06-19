#4. perform 1000CV survival LASSO on each cancer 

source("universal_LASSO_survival_script.R")

set.seed(911)

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

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
z <- which(fantom$CAT_geneName %in% rm)
rm <- fantom$CAT_geneName[z]
fantom <- fantom[-z,]

#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#----------------------------------------Analyze Results---------------------------------------------------------------
#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

results = readRDS(file="lncRNAs_100_internal_CVs_individual_cands_june19.rds")
results = do.call(rbind.data.frame, results)

for(i in 1:length(unique(results$Cancer))){
	canc_data = subset(results, results$Cancer %in% unique(results$Cancer)[i])

	#x - axis is the type column 
	#y - axis is the C-index
	#facet by lncRNA 
	colnames(canc_data)[3] = "cindex"
	canc_data$cindex = as.numeric(canc_data$cindex)
	canc_data =  as.data.table(canc_data)
	canc_data = filter(canc_data, type %in% c("ClinicalVariables", "lncRNAonly"))
	g = ggboxplot(canc_data, x="lncRNA", y="cindex", color="type") + stat_compare_means() + 
	xlab("Predictor")+
	theme_bw() + 
	geom_hline(yintercept=0.5, linetype="dashed", color = "red")
	ggpar(g, font.tickslab = c(8,"plain", "black"),
 		xtickslab.rt = 45)
	
	dev.off()


}
