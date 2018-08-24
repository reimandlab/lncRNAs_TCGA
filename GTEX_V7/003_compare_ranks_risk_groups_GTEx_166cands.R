###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(VennDiagram)
library(EnvStats)

#how many of these are candidate lncRNAs? 
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, data=="TCGA", fdr_pval <=0.05) #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

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


#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

###Data
gtex = readRDS("allGTEX_lncRNAs_scored_Aug21.rds")
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
ranked_comp = readRDS("TCGA_GTEX_lncRNAs_ranked_wDifferences_July23_noFDR.rds")

#########compare ranks-------------------------------------------------------------------

#which of the candidates are differentially expressed? 

allCands$combo = paste(allCands$gene, allCands$cancer, sep="_")
ranked_comp$combo = paste(ranked_comp$gene, ranked_comp$canc, sep="_")

results_cands = as.data.frame(matrix(ncol=ncol(ranked_comp)))
colnames(results_cands) = colnames(ranked_comp)

for(i in 1:nrow(allCands)){
	lnc = as.character(allCands$gene[i])
	canc = allCands$cancer[i]
	z = which(ranked_comp$combo %in% allCands$combo[i])
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
allCands$combo = paste(allCands$lnc, allCands$Cancer, sep="_")
colnames(results_cands)[c(1, 5)] = c("lnc", "Cancer")
allCands = merge(allCands, results_cands, by="combo")

#check if biological match 
#if HR > 1 & significantly upregulated in cancer --> potential OG? 
#if HR < 1 & significantly downregulated in cacner --> potential TS? 

allCands$reg[allCands$median_difference >=0.25] = "Upregulated_Cancer"
allCands$reg[allCands$median_difference <=-0.25] = "Downregulated_Cancer"

allCands$biological_match = ""
z = which((allCands$HR >1) & (allCands$reg == "Upregulated_Cancer"))
allCands$biological_match[z] = "PredictedOG"

z = which((allCands$HR <1) & (allCands$reg == "Downregulated_Cancer"))
allCands$biological_match[z] = "PredictedTS"
saveRDS(allCands, file="lncRNA_cands_with_GTEx_results_August21.rds")

lncs = as.list(as.character(unique(allCands$combo)))

compare_exp_boxplots = function(lnc){

	#get tumour expression and seperate by high and low expression
	z = which(allCands$combo == lnc)
	canc = unlist(strsplit(lnc, "_"))[2]
	lnc = unlist(strsplit(lnc, "_"))[1]
	#subset tumour exp
	canc_exp = subset(rna, Cancer==canc)
	canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(lnc, "Cancer", "patient"))]
	#sep by median
	median2 = median(as.numeric(canc_exp[,2]))
	canc_exp$Group = ""
	if(median2 == 0){
		z = which(canc_exp[,2] ==0)
		canc_exp$Group[z] = "Low_exp"
		canc_exp$Group[-z] = "High_exp"
	}
	if(!(median2==0)){
		z = which(canc_exp[,2] < median2)
		canc_exp$Group[z] = "Low_exp"
		canc_exp$Group[-z] = "High_exp"
	}

	#risk 
	risk = allCands$HR[allCands$lnc.x == lnc]
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
	norm_exp$Group = "GTEx"
	norm_exp$risk_type = "norm"
	norm_exp$exp_type = "norm"
	norm_exp = as.data.frame(norm_exp)

	#all expression data needed for boxplot
	all_exp = rbind(norm_exp, canc_exp)
	#all_exp[,1] = log1p(all_exp[,1])

	colnames(all_exp)[4] = "lncRNA_score"

	#get median rank difference between GTEx and high risk group
	#compare only GTEx to high risk
	risk = unique(all_exp$risk_type)[2]
	all_exp = as.data.table(filter(all_exp, Group %in% c("GTEx", risk)))

		x = all_exp[all_exp$data=="TCGA", ]
		x = x$lncRNA_score
		y = all_exp[all_exp$data=="GTEX", ]
		y = y$lncRNA_score

		t = wilcox.test(x, y, alternative = "two.sided")
		p = t$p.value
		fc = mean(x)/mean(y)
		med_diff = median(x)-median(y)

		res = c(lnc, canc, risk, med_diff, fc, p)
		names(res) = c("lnc", "canc", "risk", "median_diff", "fc", "wilcox_p")

		#boxplot
		g = ggboxplot(all_exp, ylab="lncRNA Score", x="Group", y="lncRNA_score", palette = mypal[c(3,1)], add = "jitter", fill = "Group", 
		order=c("GTEx", unique(all_exp$Group)[2]), 
		title= paste(lnc, canc, "risk=", risk))+ 
		stat_compare_means(label = "p.signif", 
                     ref.group = "GTEx") + stat_n_text() + theme_classic()
		g = ggpar(g, font.tickslab = c(10, "plain", "black"))
		print(g)

		return(res)

}

#DO NOT RUN - ALREADY DONE
#pdf("lncRNA_expression_tumours_GTEX_matched_normals_cancds_August24.pdf")
#results = llply(lncs, compare_exp_boxplots, .progress="text")
#dev.off()

#-----------SUMMARIZE---------------------------------------------------------

results2 = ldply(results)
results2$fdr = ""
results2$fdr = p.adjust(as.numeric(results2$wilcox_p), method="fdr")
results2$median_diff = as.numeric(results2$median_diff)
results2$status[results2$risk == "High_exp"] = "Unfavourable"
results2$status[results2$risk == "Low_exp"] = "Favourable"

results2$fdr_sig[results2$fdr <= 0.05] = "FDRsig"
results2$fdr_sig[results2$fdr > 0.05] = "NotFDRsig"

#check if median difference matches type of risk
results2$match = ""
results2$match[(results2$risk == "High_exp") & (results2$median_diff >= 0.2)] = "OG"
results2$match[(results2$risk == "Low_exp") & (results2$median_diff <= -0.2)] = "TS"

#get lnc names
get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

results2$name = llply(results2$lnc, get_name)
library(ggrepel)

#get cancer type 
canc_conv = readRDS("cancers_conv_july23.rds")
results2 = merge(results2, canc_conv, by="canc")

#keep only sig
results2 = as.data.table(filter(results2, fdr_sig == "FDRsig"))

pdf("new_figure5A_aug24.pdf", width=10, height=7)
ggplot(results2, aes(x=status, y=median_diff)) +
  geom_point() + geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-0.2, linetype="dashed", color = "red") + 
  scale_color_gradient2(low="grey",
                     high="blue", space ="Lab" ) + xlab("lncRNA expression type") + ylab("Median rank difference") +
  geom_label_repel(data=filter(results2, median_diff >= 0.2, fdr <= 0.05, status == "Unfavourable"), aes(label=name, fill=type), size=2) +
  geom_label_repel(data=filter(results2, median_diff <= -0.2, fdr <= 0.05, status == "Favourable"), aes(label=name, fill=type), size=2) +
  scale_fill_brewer(palette="Paired")
dev.off()


#figure 5B 
#summarize how many lncRNAs fall into which bins 
results2$match[(results2$risk == "High_exp") & (results2$median_diff <= 0)] = "High_tum_exp"
results2$match[(results2$risk == "Low_exp") & (results2$median_diff >= 0)] = "Low_tum_exp"

sum = as.data.table(table(results2$match))



















