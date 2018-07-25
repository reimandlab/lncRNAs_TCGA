###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(VennDiagram)


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

###Results
results = readRDS("results_analysis_July24.rds")

new_results = as.data.frame(matrix(ncol=8)) ; colnames(new_results) = c("gene", "fc_mean", "pval_wilcoxon", 
	"median_difference", "fdr", "fdrtag", "canc", "tis")


for(i in 1:length(results)){
	df = results[[i]]
	df$canc = cancers[i]
	df$tis = tis_match$tis[which(tis_match$cancer %in% cancers[i])]
	new_results = rbind(new_results, df)
}

new_results = new_results[-1,]

#panel of violin plots for each cancer-tissue show distribution of significantly differnt ranked 

new_results$median_difference = as.numeric(new_results$median_difference)
all_results_pre_filtering = new_results
new_results = subset(new_results, fdrtag == "FDRsig")

#try only median rank difference 25%
new_results = as.data.table(new_results)
saveRDS(new_results, file="TCGA_GTEX_lncRNAs_ranked_wDifferences_May25.rds")


#get mean order
means = as.data.table(aggregate(new_results[,4], list(new_results$canc), mean))		
	means = means[order(median_difference)]
	order= means$Group.1
pdf("summary_gtex_tcga_med_ranksdifferences_May24.pdf", width=9, height=7)
g = ggviolin(new_results, x="canc", y="median_difference", ylab="Difference in Median Ranks", order=order, fill="tis", add = "mean_sd") + theme_light() + 
geom_jitter(position=position_jitter(width=0.05,height=0),
         alpha=0.3,
         size=0.5)

ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65) + geom_hline(yintercept=c(-0.25, 0.25), linetype="dashed", color = "red") 

dev.off()

###how many lncRNAs in each cancer expressed more in cancers compared to normal 
new_results = filter(new_results, abs(median_difference) >= 0.15 )

#how many of these are candidate lncRNAs? 
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])


#how many candidates were actually evaluated?
z = which((all_results_pre_filtering$gene %in% allCands$gene) & (all_results_pre_filtering$canc %in% allCands$cancer))
genes = unique(all_results_pre_filtering$gene[z])

#how many candidates are significantly different with median change >= 0.15?
z = which((new_results$gene %in% allCands$gene) & (new_results$canc %in% allCands$cancer))
genes = unique(new_results$gene[z])
new_results = new_results[z,]

#how many up/down in cancers

all_cancers = as.data.table(table(new_results$canc))
all_cancers = all_cancers[order(N)]

pdf("summary_num_sig_different_ranked_lncRNAs_greater_than025.pdf")
g = ggbarplot(all_cancers, x="V1", y="N", color = "black", fill="steelblue", xlab="Cancer", ylab="#lncRNAs with DiffMedianRank") + theme_light()
g1 = ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65)
print(g1)
dev.off()

#broken down 

tums = filter(new_results, median_difference > 0) #1059 unique lncRNAs upregulated in cancers
canc_tums = as.data.table(table(tums$canc))
canc_tums = canc_tums[order(N)]

pdf("summary_num_sig_different_ranked_lncRNAs_greater_than025_enriched_in_tumours.pdf")
g = ggbarplot(canc_tums, x="V1", y="N", color = "black", fill="steelblue", xlab="Cancer", title= "Higher rank in TCGA", ylab="#lncRNAs with DiffMedianRank") + theme_light()
g2 = ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65)
print(g2)
dev.off()

norms = filter(new_results, median_difference <0) #86 unique lncRNAs downregulated in cancers 
canc_norms = as.data.table(table(norms$canc))
canc_norms = canc_norms[order(N)]

pdf("summary_num_sig_different_ranked_lncRNAs_greater_than025_enriched_in_normals.pdf")
g = ggbarplot(canc_norms, x="V1", y="N", color = "black", fill="steelblue", xlab="Cancer", title= "Higher rank in GTEx", ylab="#lncRNAs with DiffMedianRank") + theme_light()
g3 = ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65)
dev.off()

library(patchwork)
pdf("summary_num_sig_different_ranked_lncRNAs_greater_than025_enriched_all_three_plots.pdf")
g2 + g3 - g1 + plot_layout(ncol = 1)
dev.off()

#combine the three plots 


#161 lncRNAs in both so might be uprgulated in some cancers and downregulated in other cancers 
#compare results to results from TCGA
gtex_results = new_results
gtex_results$diff = ""
gtex_results$diff[gtex_results$median_difference > 0] = "Upregulated_Cancer"
gtex_results$diff[gtex_results$median_difference < 0] = "Downregulated_Cancer"
gtex_results$diff[gtex_results$diff == ""] = "NotDiff"

tcga_results = readRDS("cancer_normals_lncRNA_diff_expression_analysis_June5.rds")
tcga_results$reg[tcga_results$logFC >0] = "Upregulated_Cancer"
tcga_results$reg[tcga_results$logFC <0] = "Downregulated_Cancer"

#keep only fdr significant 
tcga_results = as.data.table(tcga_results)

z = which(tcga_results$Cancer %in% gtex_results$canc)
tcga_results = tcga_results[z,]

z = which(gtex_results$canc %in% tcga_results$Cancer)
gtex_results = gtex_results[z,]

#look at individual cancer type from both files and see overlap 
#how many up in gtex also up in tcga 
#how many low also low in tcga 

cancers = as.list(unique(tcga_results$Cancer))

check_overlap = function(canc){
	lncs_tcga_pos = tcga_results$lnc[which((tcga_results$Cancer == canc) & (tcga_results$reg == "Upregulated_Cancer"))]
	lncs_gtex_pos = gtex_results$gene[which((gtex_results$canc == canc) & (gtex_results$diff == "Upregulated_Cancer"))]
	overlap_pos = lncs_tcga_pos[which(lncs_tcga_pos %in% lncs_gtex_pos)]


	lncs_tcga_neg = tcga_results$lnc[which((tcga_results$Cancer == canc) & (tcga_results$reg == "Downregulated_Cancer"))]
	lncs_gtex_neg = gtex_results$gene[which((gtex_results$canc == canc) & (gtex_results$diff == "Downregulated_Cancer"))]
	overlap_neg = lncs_tcga_neg[which(lncs_tcga_neg %in% lncs_gtex_neg)]


	grid.newpage()
	draw.pairwise.venn(length(unique(lncs_tcga_pos)), length(unique(lncs_gtex_pos)), length(unique(overlap_pos)), category = c("TCGA upregulated", "GTEx upregulated"), lty = rep("blank", 
    2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.025, 2), tilte=canc)

	grid.newpage()
	draw.pairwise.venn(length(unique(lncs_tcga_neg)), length(unique(lncs_gtex_neg)), length(unique(overlap_neg)), category = c("TCGA downregulated", "GTEx downregulated"), lty = rep("blank", 
    2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.025, 2), tilte=canc)

	row = c(canc, length(unique(lncs_tcga_pos)), length(unique(lncs_gtex_pos)), length(unique(overlap_pos)), 
	 length(unique(lncs_tcga_neg)), length(unique(lncs_gtex_neg)), length(unique(overlap_neg)))
	names(row) = c("Cancer", "num_tcga_up", "num_gtex_up", "num_both_up", "num_tcga_down", "num_gtex_down", 'num_both_down')

	return(row)

}

pdf("ven_diagram_overlaps_tcga_gtex.pdf")
llply(cancers, check_overlap)
dev.off()















