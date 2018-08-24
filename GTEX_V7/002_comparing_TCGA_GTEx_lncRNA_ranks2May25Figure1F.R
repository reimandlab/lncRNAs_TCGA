###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(VennDiagram)
library(patchwork)
library(ggrepel)

###Data
gtex = readRDS("allGTEX_lncRNAs_scored_Aug21.rds")
tcga = readRDS("TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")

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

#candidates
val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
#val_cands = read.csv("112_lncRNA_cancers_combos_22_cancer_types_aug8.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 


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
cancers = as.data.frame(unique(tis_match$cancer))
colnames(cancers)[1] = "canc"
canc_conv = readRDS("cancers_conv_july23.rds")
cancers = merge(cancers, canc_conv, by="canc")

#remove cancer types with less than 50 patients
z = which(cancers$type %in% c("KICH", "CHOL", "DLBC", "UCS"))
cancers = cancers[-z,]
cancers = cancers$canc


#####------START PLOT----------------------------------------------------------------------------------------

###Results
results = readRDS("results_analysis_Aug21.rds")

new_results = as.data.frame(matrix(ncol=8)) ; colnames(new_results) = c("gene", "fc_mean", "pval_wilcoxon", 
	"median_difference", "canc", "fdr", "fdrtag", "tis")

for(i in 1:length(results)){
	df = results[[i]]
	canc = df$canc[1]
	df$tis = tis_match$tis[which(tis_match$cancer %in% canc)]
	new_results = rbind(new_results, df)
}

new_results = new_results[-1,]

#panel of violin plots for each cancer-tissue show distribution of significantly differnt ranked 
new_results$median_difference = as.numeric(new_results$median_difference)
new_results = as.data.table(new_results)
saveRDS(new_results, file="TCGA_GTEX_lncRNAs_ranked_wDifferences_July23_noFDR.rds")

all_results_pre_filtering = new_results
new_results = subset(new_results, fdrtag == "FDRsig")

#change long name to short name 
canc_conv = readRDS("cancers_conv_july23.rds")
new_results = merge(new_results, canc_conv, by="canc")

#remove cancer types with less than 50 patients
z = which(new_results$type %in% c("KICH", "CHOL", "DLBC", "UCS"))
if(!(length(z)==0)){new_results = new_results[-z,]}

#only look at those with abs median difference > 0.25
new_results = filter(new_results, abs(median_difference) >= 0.25)
new_results = as.data.table(new_results)

#get mean order
means = as.data.table(aggregate(new_results[,3], list(new_results$type), mean))		
means = means[order(fc_mean)]
order= means$Group.1
new_results$fold_change_sign = ""

#label lncRNA whether it is significantly up or downregulated in Cancer 
new_results$fold_change_sign[new_results$median_difference >=0.25] = "Up in Cancer"
new_results$fold_change_sign[new_results$median_difference <= -0.25] = "Up in Normal \nTissue"

#add gene name 
get_name = function(gene){
	name = fantom$CAT_geneName[which(fantom$CAT_geneID == gene)]
	return(name)
}

new_results$name = ""
new_results$name = unlist(llply(new_results$gene, get_name))
new_results = as.data.frame(new_results)

#Label top 5 genes in each cancer type if there is enough space
#if not label top 3 
types = unique(new_results$type)
get_best_lncs = function(tis){
	dat = new_results[which(new_results$type == tis),]
	dat$lnc_tag = ""
	dat = as.data.table(dat)
	dat = dat[order(fc_mean)]
	dat$lnc_tag[1:100] = "MostUpNormal"
	dat = dat[order(-fc_mean)]
	dat$lnc_tag[1:100] = "MostUpCancer"
	return(dat)
}

datas = llply(types, get_best_lncs)
datas = ldply(datas)
new_results = datas

#which ones global or tissue specific
num_times = as.data.table(table(new_results$name, new_results$fold_change_sign))
num_times = filter(num_times, V2 %in% c("Up in Cancer", "Up in Normal \nTissue"), N >0)
num_times = as.data.table(num_times)
num_times = num_times[order(N)]
num_times$spef = ""
num_times$spef[num_times$N==1] = "Cancer Specific"
num_times$spef[num_times$N > 1] = ">1 Cancer"

#which ones are both up and down regulated 
sum_dups = as.data.table(table(num_times$V1))
sum_dups = sum_dups[order(N)]
sum_dups = as.data.table(filter(sum_dups, N > 1))

z = which(num_times$V1 %in% sum_dups$V1)
nondups = num_times[-z,]

only_up = as.data.table(filter(nondups, V2 == "Up in Cancer"))
only_down = as.data.table(filter(nondups, V2 == "Up in Normal \nTissue"))

#remove those that appear in both up and downregulated
#These are 2,598 lncRNAs that are FDR significant and have a FC of at least 2 
z = which(num_times$V1 %in% sum_dups$V1)
num_times = num_times[-z,]

colnames(num_times)[1:3] = c("name", "fold_change_sign", "freq")
new_results = as.data.table(new_results)
new_results_extra_detail = merge(new_results, num_times, by=c("name", "fold_change_sign"))
new_results_extra_detail$fdr = -log10(new_results_extra_detail$fdr)
num_times$freq = as.factor(num_times$freq)
num_times$fold_change_sign = factor(num_times$fold_change_sign, levels = c("Up in Cancer", "Up in Normal \nTissue"))

#how many times do the 6 lncRNA candidates appear (that validated)? 

#summarize how many cancer specific vs how many multiple cancer types

#FINAL PLOT 1F
pdf("summary_gtex_tcga_updown_regulated_cancer_speicifc_or_not.pdf")
#x-axis fold_change, y-axis freq of cancer types it is apparent in, color -fold change sign
g = ggplot(num_times, aes(x=freq, fill=fold_change_sign)) + geom_histogram(stat="count") + 
scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999")) + 
theme_bw() + xlab("Number of Cancer Types") + ylab("Count") + ggtitle("Summary of Up/Downregulated lncRNAs \nbetween GTEx & TCGA")
ggpar(g, legend.title = "Difference in \nRank Expression")
dev.off()


#FINAL PLOT 1E
#pdf("summary_gtex_tcga_med_ranksdifferences_july23_mean_diff.pdf", width=9, height=7)
violins = ggviolin(new_results, x="type", y="fc_mean", ylab="Difference in Mean Ranks", order=order, add = "mean_sd") + theme_bw() + 
geom_jitter(position=position_jitter(width=0.05,height=0),
         alpha=0.3,
         size=0.5, aes(colour = fold_change_sign)) + 
         scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
          xlab("Cancer Type")
#geom_label_repel(data=filter(new_results, lnc_tag == "MostUpNormal"), aes(label=name, fill= spef), size=1.5)+
#geom_label_repel(data=filter(new_results, lnc_tag == "MostUpCancer"), aes(label=name, fill = spef), size=1.5)

violins = ggpar(violins, legend.title="Difference in \nRanked \nExpression",
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45) + geom_hline(yintercept=c(1, -1), linetype="dashed", color = "red") 

#instead of violins summarize in barplot

#order by most diff lncRNAs 
ordert = as.data.table(table(new_results$type))
ordert = ordert[order(N)]
new_results$type = factor(new_results$type, levels = ordert$V1)

bars = ggplot(new_results, aes(x=type, fill=fold_change_sign)) + geom_histogram(stat="count") + 
scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999")) + 
theme_bw() + xlab("Cancer Types") + ylab("Count")
bars = ggpar(bars, legend.title = "Difference in \nRank Expression", xtickslab.rt = 65)

#add gtex tissue covariate
tis_data = new_results[,c(8,9)]
tis_data = as.data.table(tis_data)

new_results$tis[new_results$tis == "Adre"] = "Adrenal \nGland"
new_results$tis[new_results$tis == "Brai"] = "Brain"
new_results$tis[new_results$tis == "Brea"] = "Breast"
new_results$tis[new_results$tis == "Uter"] = "Uterus"
new_results$tis[new_results$tis == "Cerv"] = "Cervix"
new_results$tis[new_results$tis == "Pros"] = "Prostate"
new_results$tis[new_results$tis == "Esop"] = "Esophagus"
new_results$tis[new_results$tis == "Kidn"] = "Kidney"
new_results$tis[new_results$tis == "Live"] = "Liver"
new_results$tis[new_results$tis == "Ovar"] = "Ovary"
new_results$tis[new_results$tis == "Panc"] = "Pancreas"
new_results$tis[new_results$tis == "Thyr"] = "Thyroid"

ordert = as.data.table(table(new_results$type, new_results$tis))

#WHAT IS GOING ON ?????? keeps throwing eerrors when running 
#line below; 

ordert = ordert[order(N)]
ordert = as.data.table(filter(ordert, N >0))
tis = unique(ordert$V2)

tis_data = tis_data[!duplicated(tis_data)]
tis_data$tisue = factor(tis_data$tis, levels = tis)

mypal = c("#B83CE0","#79E490", "#D4E7CD", "#E1B148", "#9D9887" ,"#D962C0" ,
	"#E6CACD", "#866BDF", "#D8DD8D" ,"#CDA6DC" ,"#71DCC6", "#ABE64D", "#89CBE1",
	 "#698BC5" ,"#DB6B73")

tis_data$type = factor(tis_data$type, levels = ordert$V1)

tissues = ggplot(tis_data, aes(type, 0.2)) +
    geom_tile(aes(fill = tis)) + geom_text(aes(label = tis), size=2.5) +
    theme_classic() + scale_fill_manual(values=mypal)

tissues = ggpar(tissues, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlab("GTEx Tissue")

bars + tissues + plot_layout(ncol = 1, heights = c(10, 1))

dev.off()

#which are cancer specific
new_results = merge(new_results, num_times, by = c("name", "fold_change_sign"))

#check prostate
pros = as.data.table(filter(new_results, type=="PRAD"))
#318 in total (that are only either UP or DOWN regulated not both across cancers)
pros_spef = as.data.table(filter(pros, freq ==1))
#122 specific lncRNAs 
pros_spef = pros_spef[order(abs(median_difference))]

#function to plot diff in distributions of rank bewteen tcga 
#and gtex
lnc = "ENSG00000225937"
canc = "Prostate adenocarcinoma"
get_dis_ranks = function(lnc, canc){
	t = subset(tcga, (gene == lnc) & (tis == canc))
	tissue = tis_match$tis[tis_match$cancer == canc]
	g = subset(gtex, (gene==lnc) & (tis==tissue))
	combo = rbind(t, g)
	combo$data = factor(combo$data, levels = c("GTEX", "TCGA"))
	g = ggviolin(combo, x="data", y="score",color = "data", palette = c("#E7B800","#FC4E07"),add = c("jitter", "mean_sd")) + ggtitle(paste(get_name(lnc), canc)) +
	stat_compare_means()
	print(g)
	}
pdf("pca3_gtex_vs_tumours.pdf")
get_dis_ranks(lnc, canc)
dev.off()

###how many lncRNAs in each cancer expressed more in cancers compared to normal 
z1 = which(new_results$median_difference >=0.25)
z2 = which(new_results$median_difference <= 0.25)
new_results = new_results[c(z1,z2),]

#how many of these are candidate lncRNAs? 
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, data=="TCGA", fdr_pval <=0.05) #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#how many candidates were actually evaluated? (22)
new_results = as.data.table(new_results)
new_results$combo = paste(new_results$gene, new_results$canc, sep="_")
allCands$combo = paste(allCands$gene, allCands$cancer, sep="_")

#How many were evaluated
all_results_pre_filtering$combo = paste(all_results_pre_filtering$gene, all_results_pre_filtering$canc, sep="_")
z = which(all_results_pre_filtering$combo %in% allCands$combo)
genes = unique(all_results_pre_filtering$gene[z]) #were only able to evaluate 108 unique lncRNAs 
#108 unique lncRNAs-cancers were evaluated across 10 unique cancer types 

z = which(new_results$combo %in% allCands$combo)
#how many candidates are significantly different with fold change >2 --> 31 lncRNAs 
genes = unique(new_results$gene[z]) 

#save all that were evlautaed 
z = which(all_results_pre_filtering$combo %in% allCands$combo)
all_results_pre_filtering = all_results_pre_filtering[z,]
saveRDS(all_results_pre_filtering, file="108_candidates_lncRNAs_evaluated_gtex_analysis.rds")
#only robust lncRNAs
#saveRDS(all_results_pre_filtering, file="68_candidates_lncRNAs_evaluated_gtex_analysis.rds")

z = which(new_results$combo %in% allCands$combo)
new_results = new_results[z,]
saveRDS(new_results, file="31_candidates_lncRNAs_significant_gtex_analysis.rds")
#saveRDS(new_results, file="16_candidates_lncRNAs_significant_gtex_analysis.rds")

print("done sir")


##---------------------------FIGURE 5A-----------------------------------------------------------------------------

#PLOT what the 68/112 candidates look like 

### --> all_results_pre_filtering <-- 

###Summary of 68 lncRNAs-cancer combo, the fold change between risk group and gtex group 
#plus p-value as colour, ordered by increase fold change, with tissue type covariate. 

#add HR 
cands_gtex = merge(all_results_pre_filtering, allCands, by=c("combo", "gene"))
cands_gtex$HR = as.numeric(cands_gtex$HR)
cands_gtex$HR = log2(cands_gtex$HR)
cands_gtex = merge(cands_gtex, canc_conv, by = "canc")
cands_gtex$fdr = -log10(cands_gtex$fdr)

pdf("final_figure_5A_aug22.pdf")

ggplot(cands_gtex, aes(x=median_difference, y=HR, shape=fdrtag)) + geom_vline(xintercept=0.25, linetype="dashed", color = "red") + 
geom_vline(xintercept=-0.25, linetype="dashed", color = "red") + 
  geom_point() + geom_hline(yintercept=0, linetype="dashed", color = "red") +
  scale_color_gradient2(low="grey",
                     high="blue", space ="Lab" ) + xlab("Fold Change") + ylab("log2(HR)") +
  geom_label_repel(data=filter(cands_gtex, median_difference >=0.25, fdr > -log10(0.05), HR > 0), aes(label=CAT_geneName, fill=type), size=2) +
  geom_label_repel(data=filter(cands_gtex, median_difference <= -0.25, fdr > -log10(0.05), HR < 0), aes(label=CAT_geneName, fill=type), size=2) +
  scale_fill_brewer(palette="Paired")
dev.off()



pdf("final_figure_5A_aug22_using_fold_change.pdf")
ggplot(cands_gtex, aes(x=fc_mean, y=HR, shape=fdrtag)) +
  geom_point() + geom_hline(yintercept=0, linetype="dashed", color = "red") +
  scale_color_gradient2(low="grey",
                     high="blue", space ="Lab" ) + xlab("Fold Change") + ylab("log2(HR)") +
  geom_label_repel(data=filter(cands_gtex, median_difference >=log(2), fdr > -log10(0.05), HR > 0), aes(label=CAT_geneName, fill=type), size=2) +
  geom_label_repel(data=filter(cands_gtex, median_difference <= -log(0.5), fdr > -log10(0.05), HR < 0), aes(label=CAT_geneName, fill=type), size=2) +
  scale_fill_brewer(palette="Paired")
dev.off()


#check how many cancer types a lncRNA is up/down compared to GTEx
colnames(num_times)[1] = "CAT_geneName"
cands_times = merge(cands_gtex, num_times, by="CAT_geneName")

##----------------------------DONE---------------------------------------------------------------------------------
























































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















