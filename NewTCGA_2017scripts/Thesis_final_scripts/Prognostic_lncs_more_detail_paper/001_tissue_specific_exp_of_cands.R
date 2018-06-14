#4. perform 1000CV survival LASSO on each cancer 

source("universal_LASSO_survival_script.R")

set.seed(911)

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")


#all RNA, PCG and clinical data is lwoded in already? 
#clinical columns attached to rna file **

#Four things we want to do here for now:

#1. check correlation of each lncRNA candidate with available clinical variables for cancer type

#2. calculate c-index distribution for each lncRNA prognostic marker comapred to clinical variables available 

#3. visualize results of each lncRNA candidate via forest plot to show differences in HR and p-value between all predictors 

#4. train model using all lncRNA predictors obtained from TCGA analysis and test performance of model on PCAWG, 50% of PCAWG at a time randomly selected --> train model and get c-index distribution , compare to age? whatever clinical data available? 


#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

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
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]


#Analysis--------------------------------------------

#For each lncRNA candidate evaluate how its expression varies from its risk cancer type
#and other cancer types, is it able to differenitate between them? 

get_tissue_specificity = function(lnc){

	#first just plot boxplot of expression ordered least to highest and 
	#print cancer type 
	z1 = which(colnames(rna) %in% lnc)
	z2 = which(colnames(rna) %in% c("type", "Cancer"))
	lnc_rna = rna[,c(z1, z2)]
	lnc_rna[,1] = log1p(lnc_rna[,1])
	lnc_canc = allCands$cancer[allCands$gene == lnc]
	lnc_canc = unique(lnc_rna$type[lnc_rna$Cancer == lnc_canc])

	print(lnc)
	print(lnc_canc)

	colnames(lnc_rna)[1] = "lncRNAexpression"

	#get order - get cancer type median
	means = as.data.table(aggregate(lnc_rna[, 1], list(lnc_rna$type), mean))
	means = means[order(x)]

	order= unlist((means[,1]))

	#is cancer with highest mean expression the candidate cancer type?
	check = order[[32]] == lnc_canc

	if(length(lnc_canc) >1){
		lnc_canc_combo = paste(lnc_canc, collapse="_")
	}

	if(length(lnc_canc)==1){
		lnc_canc_combo = lnc_canc
	}

	lnc_rna$type = factor(lnc_rna$type, levels = order)

	#boxplot

    g = ggboxplot(lnc_rna, x = "type", y="lncRNAexpression", fill="type", color="black") + 
    geom_boxplot(data=lnc_rna[lnc_rna$type==lnc_canc,],
                        aes(x = type, y = lncRNAexpression),fill="seashell", color="turquoise4")+
	theme_bw() + ggtitle(paste(lnc, "risk cancer is", lnc_canc_combo)) + xlab("Cancer")
	ggpar(g,
 		font.tickslab = c(7,"plain", "black"),
 		xtickslab.rt = 45, legend="none")
}


lncs = as.list(unique(allCands$gene)) #for those in multiple cancers just 
#label in the title or something that it's a cand in multiple cancers 

pdf("166_candidate_lncRNAs_expression_Across_cancers.pdf")
llply(lncs, get_tissue_specificity, .progress="text")
dev.off()

##--------Summary Heatmap via Geom Tile--------------------------------------------------------------------------------------


#one column is lncRNA name 
#one column is cancer type median expression
#one column for prognostic 

get_tissue_specificity = function(lnc){

	#first just plot boxplot of expression ordered least to highest and 
	#print cancer type 
	z1 = which(colnames(rna) %in% lnc)
	z2 = which(colnames(rna) %in% c("type", "Cancer"))
	lnc_rna = rna[,c(z1, z2)]
	lnc_rna[,1] = log1p(lnc_rna[,1])
	lnc_canc = allCands$cancer[allCands$gene == lnc]
	lnc_canc = unique(lnc_rna$type[lnc_rna$Cancer == lnc_canc])

	print(lnc)
	print(lnc_canc)

	colnames(lnc_rna)[1] = "lncRNAexpression"

	#get order - get cancer type median
	means = as.data.frame(aggregate(lnc_rna[, 1], list(lnc_rna$type), mad))
	colnames(means) = c("Cancer", "MeanExpression")

	means$lncrna = lnc

	#if(length(lnc_canc) >1){
	#lnc_canc = paste(lnc_canc, collapse="_")
	#}

	check_canc = function(canc){
		z = which(lnc_canc %in% canc)
		if(length(z) ==0){
			return(NA)
		}
		if(!(length(z)==0)){
			return(lnc_canc[z])
		}
	} 

	means$risk_canc =""
	means$risk_canc = unlist(llply(means$Cancer, check_canc))

	return(means)

}


heatmap_set_up = llply(lncs, get_tissue_specificity, .progress="text")
heatmap_set_up = ldply(heatmap_set_up, data.frame)
#heatmap_set_up$MeanExpression = floor(heatmap_set_up$MeanExpression)
heatmap_set_up[,3] = as.character(heatmap_set_up[,3])

colnames(fantom)[1] = "lncrna"
fantom = merge(fantom, heatmap_set_up, by="lncrna")


#order cancers by highest to least number of lncrna candidates
cancers = as.data.table(table(fantom$risk_canc))
cancers = cancers[order(-N)]

order= cancers$V1
fantom$Cancer = factor(fantom$Cancer, levels = order)

pdf("166_lncrna_candidates_heatmap_mean_Expression.pdf", width=8, height=9)

g = ggplot(fantom, aes(Cancer, CAT_geneName)) +
  geom_tile(aes(fill = MeanExpression)) +
  geom_text(aes(label = risk_canc), size=1.5) +
    scale_fill_gradient(low = "white", high = "orange", na.value = 'transparent') +
    xlab("Cancer") + ylab("lncRNA") + theme_bw() 
ggpar(g,
 font.tickslab = c(4, "plain", "black"), 
 xtickslab.rt = 45, legend.title="MAD value")

dev.off()































































































