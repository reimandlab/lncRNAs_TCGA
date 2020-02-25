#------------------------------------------------------------------------------
#This script takes all the PCGs that are significantly differentially expressed
#between lncRNA risk groups and calculates their correlation 
#plots heatmap of this data as well as lncRNA and pcg prognostic
#relationship 
#------------------------------------------------------------------------------

set.seed(911)

#load all libraries and functions 
source("check_lnc_exp_cancers.R")
source("survival_script_march14_already_labelled_highlow.R")
#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#old = readRDS("old_168_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
#val_cands = as.data.table(val_cands)
#val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
#val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
#val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
#all <- merge(rna, pcg, by = c("patient", "Cancer"))
#all = all[,1:25170]

#canc conversion
#canc_conv = rna[,c(which(colnames(rna) %in% c("Cancer", "type")))]
#canc_conv = canc_conv[!duplicated(canc_conv),]
#colnames(canc_conv) = c("type", "cancer")

#------DATA-----------------------------------------------------

multivariate_lncs = fread("summary_cross_validation_clinical_multivariate_models_lncRNAs_Nov1.csv")
single_lncs = fread("summary_cross_validation_clinical_random_lncRNAs_Nov1.txt")

#------ANALYSIS-------------------------------------------------

#for each cancer type 
#generate barplot 
#x-axis: clinica, lnc1, lnc2, lnc3, lnc4...., all lncs, all lncs+clin 
#y-axis: median c-index 

cancers = unique(allCands$canc)

get_bar = function(canc){
	
	multi_dat = as.data.table(filter(multivariate_lncs, cancer == canc))
	single_dat = as.data.table(filter(single_lncs, cancer == canc))

	#make input for plot 
	lncs = single_dat[,c("lncRNA", "median_lncRNA")]
	lncs$type = "lncRNA"

	lncs_all_combo=c("lncs_all_combo", multi_dat$clin_combo_cind_med)
	lncs_all=c("lncs_all", multi_dat$median_lncRNA)
	clin=c("clinical", multi_dat$median_clinical)
	all = as.data.frame(rbind(lncs_all_combo,  lncs_all, clin))
	colnames(all) = c("lncRNA", "median_lncRNA")
	rownames(all) = NULL
	all$type = c("lncs_clin", "lncs_all", "clin")
	all=rbind(lncs, all)
	all$median_lncRNA = as.numeric(all$median_lncRNA)
	all = all[order(-median_lncRNA)]
	all$lncRNA = factor(all$lncRNA, levels=all$lncRNA)
	g = ggbarplot(all, x="lncRNA", y="median_lncRNA", fill="type") +ylim(0,1)+
		geom_hline(yintercept=0.5, linetype="dashed", color = "red")+ylab("Median c-index")
	g = ggpar(g, x.text.angle=45, font.xtickslab = c(7, "plain", "black"))+ggtitle(canc)
	print(g)
}

pdf("cancers_multivaraite_models.pdf")
llply(cancers, get_bar)
dev.off()













