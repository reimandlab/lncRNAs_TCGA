###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(VennDiagram)
library(EnvStats)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

#cands -- ion channels 
allCands = read.csv("ION_CHANNELS_targets_and_families.csv")
allCands = merge(allCands, ucsc, by = "HGNC.symbol")

###Data
gtex = readRDS("allGTEX_ionchannels_scored_Jan23.rds")
tcga = readRDS("TCGA_all_ionchannels_cancers_scored_byindexDec30.rds")

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

#2. type of cancers with tissues available -> seperate into dataframes
cancers = as.data.frame(unique(tis_match$cancer))
colnames(cancers)[1] = "canc"
canc_conv = readRDS("cancers_conv_july23.rds")
cancers = merge(cancers, canc_conv, by="canc")

#remove cancer types with less than 50 patients
z = which(cancers$type %in% c("KICH", "CHOL", "DLBC", "UCS"))
cancers = cancers[-z,]
cancers = cancers$canc
cancers = tis_match$tis
tis_match$combo = paste(tis_match$cancer, tis_match$tis, sep = "_")
cancers = tis_match$combo

#expression data 
rna = readRDS("19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")

#tcga gtex comparison - MAIN FILE 
ranked_comp = readRDS("TCGA_GTEX_ionchannels_ranked_wDifferences_Dec30_noFDR.rds")

#########compare ranks-------------------------------------------------------------------

#which of the candidates are differentially expressed? 

allCands$combo = allCands$hg19.ensGene.name2
ranked_comp$combo = ranked_comp$gene

results_cands = ranked_comp

#check if biological match 
#if HR > 1 & significantly upregulated in cancer --> potential OG? 
#if HR < 1 & significantly downregulated in cacner --> potential TS? 

results_cands$reg[results_cands$median_difference >=0.15] = "Upregulated_Cancer"
results_cands$reg[results_cands$median_difference <=-0.15] = "Downregulated_Cancer"

surv_results = readRDS("TCGA_ION_CHANNEL_results_Sept21.rds")
surv_results$combo = paste(surv_results$gene, surv_results$cancer, sep="_")
results_cands$combo = paste(results_cands$gene, results_cands$canc, sep="_")

combined = merge(surv_results, results_cands, by="combo")
combined$gene.y = NULL
colnames(combined)[2] = "gene"
combined$HR = as.numeric(combined$HR)

combined$biological_match = ""
z = which((combined$HR >1) & (combined$reg == "Upregulated_Cancer"))
combined$biological_match[z] = "PredictedOG"

z = which((combined$HR <1) & (combined$reg == "Downregulated_Cancer"))
combined$biological_match[z] = "PredictedTS"
saveRDS(combined, file="ION_CHANNELS_cands_with_GTEx_results_Jan23.rds")

combined$comboo = paste(combined$gene, combined$canc, combined$tis, sep="_")

#lets look only at the ones that have some significnat prognostic signal 
combined = as.data.table(filter(combined, pval < 0.05))

lncs = as.list(as.character(unique(combined$comboo)))
saveRDS(combined, file="ICs_wSurvival_GTEx_info_Jan24.rds")

compare_exp_boxplots = function(lnc){
	print(lnc)
	combo = lnc
	#get tumour expression and seperate by high and low expression
	z = which(combined$combo == lnc)
	cance = unlist(strsplit(lnc, "_"))[2]
	lnc = unlist(strsplit(lnc, "_"))[1]
	tiss = unlist(strsplit(combo, "_"))[3]
	#subset tumour exp
	cancc = canc_conv$type[canc_conv$canc == cance]
	canc_exp = subset(rna, type==cancc)
	canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(lnc, "type", "patient"))]
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
	risk = combined$HR[combined$combo == paste(lnc, cance, sep="_")]
	if(as.numeric(risk) > 1){
		risk_exp = "High_exp"
	}
	if(as.numeric(risk) < 1){
		risk_exp = "Low_exp"
	}

	canc_exp$risk_type = risk_exp
	canc_exp$exp_type = "Tumour"

	#subset tcga to just lncrna
	tcga$combo = paste(tcga$gene, tcga$tis, sep="_")
	tcga_lnc = subset(tcga, combo == paste(lnc, cance, sep="_"))

	#get ranks for patients and compare to ranks of lncRNA in GTEx 
	canc_exp = merge(tcga_lnc, canc_exp, by=c("patient"))	

	#get gtex ranks for same tissue (make sure it's the specific one)###***** <-------
	tis = tiss
	print(tis)

	#canc_exp = canc_exp[,c(1:6, 9:11)]
	z = which(gtex$tis == tis)
	norm_exp = gtex[z,]
	norm_exp = subset(norm_exp, gene == lnc)

	#norm_exp = norm_exp[,c(4, 1, 2, 3, 6, 5)]
	norm_exp$Group = "GTEx"
	norm_exp$risk_type = "norm"
	norm_exp$exp_type = "norm"
	norm_exp = as.data.frame(norm_exp)

	canc_exp = canc_exp[,c("patient", "gene", "exp", "score", "tis", "data", "Group", "risk_type", "exp_type")]
	norm_exp = norm_exp[,c("patient", "gene", "exp", "score", "tis", "data", "Group", "risk_type", "exp_type")]

	#all expression data needed for boxplot
	all_exp = rbind(norm_exp, canc_exp)
	#all_exp[,1] = log1p(all_exp[,1])

	colnames(all_exp)[4] = "lncRNA_score"

	#get median rank difference between GTEx and high risk group
	#compare only GTEx to high risk
	risk = unique(all_exp$risk_type)[2]
	all_exp_for_plot = all_exp
	all_exp = as.data.table(filter(all_exp, Group %in% c("GTEx", risk)))

		x = all_exp[all_exp$data=="TCGA", ]
		x = x$lncRNA_score
		y = all_exp[all_exp$data=="GTEX", ]
		y = y$lncRNA_score

		t = wilcox.test(x, y, alternative = "two.sided")
		p = t$p.value
		fc = mean(x)/mean(y)
		med_diff = median(x)-median(y)

		res = c(lnc, canc, risk, med_diff, fc, p, tiss)
		names(res) = c("lnc", "canc", "risk", "median_diff", "fc", "wilcox_p", "tissue")
		all_exp_for_plot$Group = factor(all_exp_for_plot$Group, levels = c("GTEx", "High_exp", "Low_exp"))

		#boxplot
		g = ggboxplot(all_exp_for_plot, ylab="lncRNA Score", x="Group", y="lncRNA_score", palette = mypal[c(3,1,2)], add = "jitter", fill = "Group", 
		order=c("GTEx", "High_exp", "Low_exp"), 
		title= paste(lnc, tiss, cancc, "\nrisk=", risk))+ 
		stat_compare_means(label = "p.signif", 
                     ref.group = "GTEx") + stat_n_text() + theme_classic()
		g = ggpar(g, font.tickslab = c(14, "plain", "black"), legend="none")
		print(g)

		return(res)

}

#DO NOT RUN - ALREADY DONE
pdf("IonChannels_expression_tumours_GTEX_matched_normals_cancds_Jan23.pdf")
results = llply(lncs, compare_exp_boxplots, .progress="text")
dev.off()

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
results2$match[(results2$risk == "High_exp") & (results2$median_diff >= 0.1)] = "OG"
results2$match[(results2$risk == "Low_exp") & (results2$median_diff <= -0.1)] = "TS"

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
results2 = as.data.table(results2)
results2[,11]=as.character(results2[,11])
saveRDS(results2, file="110_lncRNAs_wgtex_data_nov16.rds")
write.csv(results2, file="110_lncRNAs_wgtex_data_nov16.csv", quote=F,row.names=F)

#keep only sig
#results2 = as.data.table(filter(results2, fdr_sig == "FDRsig"))

pdf("new_figure5A_aug24.pdf", width=7, height=5)
ggplot(results2, aes(x=status, y=median_diff)) +
  geom_point() + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  scale_color_gradient2(low="grey",
                     high="blue", space ="Lab" ) + xlab("lncRNA expression type") + ylab("Median rank difference \nRisk Group - GTEx") +
  geom_label_repel(data=filter(results2, median_diff >= 0.2, fdr <= 0.05, status == "Unfavourable"), aes(label=name, fill=type), size=1.5) +
  geom_label_repel(data=filter(results2, median_diff <= -0.2, fdr <= 0.05, status == "Favourable"), aes(label=name, fill=type), size=1.5) +
  scale_fill_brewer(palette="Paired")
dev.off()


#figure 5B 
#summarize how many lncRNAs fall into which bins 
results2$match[(results2$risk == "High_exp") & (results2$median_diff <= 0)] = "High_tum_exp"
results2$match[(results2$risk == "Low_exp") & (results2$median_diff >= 0)] = "Low_tum_exp"

sum = as.data.table(table(results2$match))



















