###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(VennDiagram)

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

###Data
gtex = readRDS("allGTEX_lncRNAs_scored_May23.rds")
tcga = readRDS("TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")
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


#expression data 
rna = readRDS("TCGA_rna_expression_Data_forgtex_analysis.rds")

#tcga gtex comparison 
ranked_comp = readRDS("TCGA_GTEX_lncRNAs_ranked_wDifferences_May25.rds")

#########compare ranks-------------------------------------------------------------------

#which of the candidates are differentially expressed? 

results_cands = as.data.frame(matrix(ncol=ncol(ranked_comp)))
colnames(results_cands) = colnames(ranked_comp)

for(i in 1:nrow(allCands)){
	lnc = as.character(allCands$gene[i])
	canc = allCands$cancer[i]
	z = which((ranked_comp$gene == lnc) & (ranked_comp$canc == canc))
	print(z)
	row = ranked_comp[z,]
	names(row) = colnames(results_cands)
	results_cands = rbind(results_cands, row)
}

results_cands = results_cands[-1,]

#22 lncRNAs from candidates coming from cancer types that have available 
#matched gtex tissues are significantly differetnially expressed 
#across 9 cancer types 

colnames(allCands)[c(1, 5)] = c("lnc", "Cancer")
colnames(results_cands)[c(1, 7)] = c("lnc", "Cancer")
allCands = merge(allCands, results_cands, by=c("lnc", "Cancer"))


#check if biological match 
#if HR > 1 & significantly upregulated in cancer --> potential OG? 
#if HR < 1 & significantly downregulated in cacner --> potential TS? 

allCands$reg[allCands$median_difference >0] = "Upregulated_Cancer"
allCands$reg[allCands$median_difference <0] = "Downregulated_Cancer"

allCands$biological_match = ""
z = which((allCands$HR >1) & (allCands$reg == "Upregulated_Cancer"))
allCands$biological_match[z] = "PredictedOG"

z = which((allCands$HR <1) & (allCands$reg == "Downregulated_Cancer"))
allCands$biological_match[z] = "PredictedTS"
saveRDS(allCands, file="lncRNA_cands_with_GTEx_results_June22.rds")

lncs = as.list(as.character(unique(allCands$lnc)))

compare_exp_boxplots = function(lnc){

	#get tumour expression and seperate by high and low expression
	z = which(allCands$lnc == lnc)
	canc = allCands$Cancer[z]
	#subset tumour exp
	canc_exp = subset(rna, Cancer==canc)
	canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(lnc, "Cancer", "patient"))]
	#sep by median
	median2 = median(as.numeric(canc_exp[,2]))
	canc_exp$tag = ""
	if(median2 == 0){
		z = which(canc_exp[,2] ==0)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}
	if(!(median2==0)){
		z = which(canc_exp[,2] < median2)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}

	#risk 
	risk = allCands$HR[allCands$lnc == lnc]
	if(as.numeric(risk) > 1){
		risk_exp = "High_exp"
	}
	if(as.numeric(risk) < 1){
		risk_exp = "Low_exp"
	}

	canc_exp$risk_type = risk_exp
	canc_exp$exp_type = "Tumour"

	#subset tcga to just lncrna
	tcga_lnc = subset(tcga, gene == lnc)

	#get ranks for patients and compare to ranks of lncRNA in GTEx 
	canc_exp = merge(tcga_lnc, canc_exp, by=c("patient"))	

	#get gtex ranks for same tissue
	tis = tis_match$tis[tis_match$cancer %in% canc_exp$Cancer]

	canc_exp = canc_exp[,c(1:6, 9:11)]

	z = which(gtex$tis == tis)
	norm_exp = gtex[z,]
	norm_exp = subset(norm_exp, gene == lnc)

	norm_exp = norm_exp[,c(4, 1, 2, 3, 6, 5)]
	norm_exp$tag = "norm"
	norm_exp$risk_type = "norm"
	norm_exp$exp_type = "norm"
	norm_exp = as.data.frame(norm_exp)

	#all expression data needed for boxplot
	all_exp = rbind(norm_exp, canc_exp)
	#all_exp[,1] = log1p(all_exp[,1])

	colnames(all_exp)[4] = "lncRNA_score"

	#boxplot
	ggboxplot(all_exp, ylab="lncRNA Score", x="tag", y="lncRNA_score", palette = mypal[c(3,2,1)], add = "jitter", fill = "tag", order=c("norm", "Low", "High"), 
		title= paste(lnc, canc, "risk=", risk_exp))+ 
		stat_compare_means(label = "p.signif", 
                     ref.group = "norm") + theme_light()  

}

pdf("lncRNA_expression_tumours_GTEX_matched_normals_cancds_June5.pdf", width=10)
llply(lncs, compare_exp_boxplots, .progress="text")
dev.off()
