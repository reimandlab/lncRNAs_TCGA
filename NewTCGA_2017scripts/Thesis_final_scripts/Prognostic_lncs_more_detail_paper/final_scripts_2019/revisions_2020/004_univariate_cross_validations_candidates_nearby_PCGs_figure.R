set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files

#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs

#which cancer types are the non-unique lncRNAs from?
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_symbol")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
allCands$combo=paste(allCands$gene, allCands$cancer, sep="_")
allCands$gene_name = allCands$gene_symbol

#------DATA-----------------------------------------------------

single_lncs = fread("summary_cross_validation_clinical_random_lncRNAs_Nov1.txt")

#single pcgs
single_pcgs = readRDS("nearby_pcgs_1000_internal_CVs_individual_cands_june19.rds")
single_pcgs = as.data.table(ldply(single_pcgs))
colnames(single_pcgs)=c("pcg", "cancer", "cindex", "type", "round")
single_pcgs = as.data.table(filter(single_pcgs, type == "lncRNAonly"))

single_lncs = as.data.table(filter(single_lncs, type == "lncRNA_canc"))
single_lncs = single_lncs[,c("lncRNA", "cancer", "cindex", "type", "round")]

r = readRDS("127_cis_antisense_pairs_survival_results_10kb_nov16.rds")

#get ensgs
get_ensg_pcg = function(gene){
	ensg = filter(hg38, symbol == gene)$ensgene[1]
	return(ensg)
}

get_ensg_lnc = function(gene){
	fantom=as.data.table(fantom)
	z = which(fantom$CAT_geneName == gene)
	ensg = fantom$gene[z]
	return(ensg)
}

r$pcg = sapply(r$pcg, get_ensg_pcg)
r$lnc = sapply(r$lnc, get_ensg_lnc)
colnames(canc_conv)[1]="canc"
r=merge(r, canc_conv)

#------ANALYSIS-------------------------------------------------

#for each cancer type
#generate barplot
#x-axis: clinica, lnc1, lnc2, lnc3, lnc4...., all lncs, all lncs+clin
#y-axis: median c-index

combos = paste(r$pcg, r$lnc, r$Cancer, sep="_")

get_bar = function(combo){

	print(combo)
	#canc
	canc_type = unlist(strsplit(combo, "_"))[3]

	#get lnc
	canc_lnc = unlist(strsplit(combo, "_"))[2]

	#get pcg
	canc_pcg = unlist(strsplit(combo, "_"))[1]

	lnc_cinds = filter(single_lncs, lncRNA == canc_lnc, cancer==canc_type, !(is.na(cindex)))
	pcg_cinds = filter(single_pcgs, pcg == canc_pcg, cancer==canc_type, !(is.na(cindex)))

	if(!(dim(pcg_cinds)[1]==0)){
	med_lnc = median(as.numeric(lnc_cinds$cindex))
	med_pcg = median(as.numeric(pcg_cinds$cindex))

	pval = as.numeric(tidy(wilcox.test(as.numeric(pcg_cinds$cindex), as.numeric(lnc_cinds$cindex)))[2])

	row = c(canc_type, canc_lnc, canc_pcg, med_lnc, med_pcg, pval)
	return(row)

}}

all_res = as.data.table(ldply(llply(combos, get_bar)))
colnames(all_res) = c("canc_type", "canc_lnc", "canc_pcg", "med_lnc", "med_pcg", "pval")
all_res$fdr = p.adjust(all_res$pval)
all_res$lnc_better = ""
z=which((all_res$med_lnc > all_res$med_pcg) & (all_res$fdr < 0.05))
all_res$lnc_better[z] = "yes"

all_res$med_lnc = as.numeric(all_res$med_lnc)
all_res$med_pcg = as.numeric(all_res$med_pcg)
colnames(all_res)[1] = "Cancer"
all_res = merge(canc_conv, all_res)
colnames(all_res)[2]="type"

saveRDS(all_res, file="lncRNA_vs_nearby_pcgs_cindices_cross_validations.rds")

pdf("/u/kisaev/Jan2021/figure2D_lncs_pcgs_antisense_10kb_cross_val_cindices.pdf", width=5,height=5)
g = ggplot(all_res, aes(med_pcg, med_lnc)) +
  geom_point(aes(colour=type,
     shape=lnc_better), size=1.75, show.legend = FALSE) +
     scale_shape_manual(values = c(5, 17)) +
     colScale+
    xlab("Neighbour PCG Concordance") + ylab("lncRNA Concordance") +
    theme(axis.text = element_text(size=13),
    legend.position = "none") +
     xlim(0.3,0.9) + ylim(0.3,0.9) + geom_abline(intercept=0) +
     geom_vline(xintercept=0.5, linetype="dashed", color = "red") +
     geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
     theme_bw()
     #geom_text_repel(data = subset(prog_pcgs, lncConcordance > 0.7 | pcgConcordance > 0.65), size=2, nudge_y = 0.1,
      #direction = "x",segment.color = "grey50",
      #segment.size = 0.05)
g
dev.off()

#combine data for supplementary tables
r$combo_merge = paste(r$lnc, r$pcg, r$canc, sep="_")
all_res$combo_merge = paste(all_res$canc_lnc, all_res$canc_pcg, all_res$type, sep="_")
both = merge(r, all_res, by="combo_merge")
write.table(both, file="/u/kisaev/Jan2021/nearby_pcgs_results.txt", sep="}", quote=F, row.names=F)
