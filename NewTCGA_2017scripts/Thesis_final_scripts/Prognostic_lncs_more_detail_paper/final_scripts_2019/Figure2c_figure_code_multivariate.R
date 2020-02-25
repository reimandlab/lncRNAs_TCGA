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

################################################################################
#cindices from 166 candidates 
lnc_cands = readRDS("lncRNAs_1000_internal_CVs_individual_cands_multivariate_feb25.rds")
r = ldply(lnc_cands)
r = as.data.table(r)
#z = which(r$type == "lncRNA&clin")
#r = r[-z,]

z = which(r$type == "lncRNAonly")
r = r[z,]

#remove those that are no longer considred in our analysis
#r$combo = paste(r$lncRNA, r$Cancer, sep="_")
#z1 = which(r$combo %in% allCands$combo)
#z2 = which(!(str_detect(r$lncRNA, "ENSG")))
#r = r[c(z1,z2),]

list_r = split(r, by = "Cancer")

get_sum = function(canc){
	#get summ for c-indices for each lnc and clin
	colnames(canc)[3] = "cindex"
	z = which(is.na(canc$cindex))
	if(!(length(z)==0)){
	canc = canc[-z,]}
	canc$cindex = as.numeric(canc$cindex)
	lncs <- group_by(canc, lncRNA)
	s = as.data.table(summarise(lncs, mean_cindex = mean(cindex), median_cindex = median(cindex),sd=sd(cindex)))
	clin = s[nrow(s),]
	s$imp = s$median_cindex - clin$median_cindex
	return(s)
}

lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
#z= which(str_detect(lists$lncRNA, "ENSG"))
lncs_cands_lists = lists
lncs_cands_all_cindices = r

#---------------------------------------------------------------------------------

################################################################################
#cindices from clinical variables
lnc_rand = readRDS("lncRNAs_1000_internal_CVs_individual_cands_multivariate_feb25.rds")
r = ldply(lnc_rand)
r = as.data.table(r)
z = which(r$type == "ClinicalVariables")
r = r[z,]

#remove those that are no longer considred in our analysis
#r$combo = paste(r$lncRNA, r$Cancer, sep="_")
#z1 = which(r$combo %in% allCands$combo)
#z2 = which(!(str_detect(r$lncRNA, "ENSG")))
#r = r[c(z1,z2),]

list_r = split(r, by = "Cancer")

get_sum_clin = function(canc){
	#get summ for c-indices for each lnc and clin
	colnames(canc)[3] = "cindex"
	z = which(is.na(canc$cindex))
	if(!(length(z)==0)){
	canc = canc[-z,]}
	canc$cindex = as.numeric(canc$cindex)
	lncs <- group_by(canc, lncRNA)
	s = as.data.table(summarise(lncs, mean_cindex = mean(cindex), median_cindex = median(cindex), sd=sd(cindex)))
	clin = s[nrow(s),]
	s$imp = s$median_cindex - clin$median_cindex
	return(s)
}

lists = llply(list_r, get_sum_clin)
lists = ldply(lists)
lists = as.data.table(lists)
#z= which(str_detect(lists$lncRNA, "ENSG"))
clinical_lists = lists
clinical_cindices = r

################################################################################
#cindices from lncs AND clinical 
lnc_cands = readRDS("lncRNAs_1000_internal_CVs_individual_cands_multivariate_feb25.rds")
r = ldply(lnc_cands)
r = as.data.table(r)
#z = which(r$type == "lncRNA&clin")
#r = r[-z,]

z = which(r$type == "lncRNA&clin")
r = r[z,]

#remove those that are no longer considred in our analysis
#r$combo = paste(r$lncRNA, r$Cancer, sep="_")
#z1 = which(r$combo %in% allCands$combo)
#z2 = which(!(str_detect(r$lncRNA, "ENSG")))
#r = r[c(z1,z2),]

list_r = split(r, by = "Cancer")

get_sum = function(canc){
	#get summ for c-indices for each lnc and clin
	colnames(canc)[3] = "cindex"
	z = which(is.na(canc$cindex))
	if(!(length(z)==0)){
	canc = canc[-z,]}
	canc$cindex = as.numeric(canc$cindex)
	lncs <- group_by(canc, lncRNA)
	s = as.data.table(summarise(lncs, mean_cindex = mean(cindex), median_cindex = median(cindex),sd=sd(cindex)))
	clin = s[nrow(s),]
	s$imp = s$median_cindex - clin$median_cindex
	return(s)
}

lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
#z= which(str_detect(lists$lncRNA, "ENSG"))
lncs_cands_lists = lists
lnc_clin_cindices = r

#---------------------------------------------------------------------------------

#for each lncRNA-cancer pair 
#compare distribution of cindices via two-sided U test between lncRNA candidate
#and random lncRNAs 

allCands = filter(allCands, data == "TCGA")
combos = (unique(allCands$cancer))

get_comparison = function(comboo){
	
	#get lncRNA candidate c-index 
	lnc_cind = as.data.table(filter(lncs_cands_all_cindices, Cancer == comboo))
	
	#how many runs works
	runs = dim(lnc_cind)[1]
	lnc_cind$type = "lncRNA_canc"
	
	#get random lncs 
	#make sure they dont include any of the candidates 
	canc = comboo
	
	#only genes
	colnames(lnc_cind)[3] = "cindex"
	z = which(is.na(lnc_cind$cindex))
	if(!(length(z)==0)){
	lnc_cind = lnc_cind[-z,]}

	#get clinical cindices 
	clinical_cinds = as.data.table(filter(clinical_cindices, Cancer == canc))
	colnames(clinical_cinds)[3] = "cindex"
	z = which(is.na(clinical_cinds$cindex))
	if(!(length(z)==0)){
	clinical_cinds = clinical_cinds[-z,]}
	clinical_cinds$type = "clinical_variables"
#	clinical_cinds$combo = paste(clinical_cinds$lncRNA, clinical_cinds$Cancer, sep="_")

	#combo cinds
	combo_cinds = as.data.table(filter(lnc_clin_cindices, Cancer == canc))
	colnames(combo_cinds)[3] = "cindex"
	z = which(is.na(combo_cinds$cindex))
	if(!(length(z)==0)){
	combo_cinds = combo_cinds[-z,]}
	combo_cinds$type = "lncRNA&clin"
#	clinical_cinds$combo = paste(clinical_cinds$lnc

	lnc_cind$cindex = as.numeric(lnc_cind$cindex)
	clinical_cinds$cindex = as.numeric(clinical_cinds$cindex)
	combo_cinds$cindex = as.numeric(combo_cinds$cindex)

	lnc_med = median(lnc_cind$cindex)
	clin_med = median(clinical_cinds$cindex)
	combo_med = median(combo_cinds$cindex)
	
	#rbind lnc candidate and random lncs 
	all_cindices = rbind(lnc_cind, clinical_cinds, combo_cinds)
	all_cindices$cindex = as.numeric(all_cindices$cindex)
	
	#get pvalue 
	#w = wilcox.test(all_cindices$cindex[all_cindices$type=="lncRNA_canc"], all_cindices$cindex[all_cindices$type=="lncRNA_random"], alternative="greater")
	#wp_lnc_random = glance(w)[2]
	#diff_meds_lnc_random = median(all_cindices$cindex[all_cindices$type=="lncRNA_canc"]) - median(all_cindices$cindex[all_cindices$type=="lncRNA_random"])

	w = wilcox.test(all_cindices$cindex[all_cindices$type=="lncRNA_canc"], all_cindices$cindex[all_cindices$type=="clinical_variables"], alternative="greater")
	wp_lnc_clinical = glance(w)[2]
	diff_meds_lnc_clinical = median(all_cindices$cindex[all_cindices$type=="lncRNA_canc"]) - median(all_cindices$cindex[all_cindices$type=="clinical_variables"])
	
	w = wilcox.test(all_cindices$cindex[all_cindices$type=="lncRNA&clin"], all_cindices$cindex[all_cindices$type=="clinical_variables"], alternative="greater")
	wp_lnc_clinical_combo = glance(w)[2]
	diff_meds_lnc_clinical_combo = median(all_cindices$cindex[all_cindices$type=="lncRNA&clin"]) - median(all_cindices$cindex[all_cindices$type=="clinical_variables"])

	all_cindices$type = factor(all_cindices$type, levels=c("clinical_variables", "lncRNA_canc", "lncRNA&clin"))

	# Visualize: Specify the comparisons you want
	#my_comparisons <- list( c("clinical_variables", "lncRNA_canc"), c("lncRNA_canc", "lncRNA_random"),
	# c("clinical_variables", "lncRNA_random"),c("clinical_variables", "lncRNA&clin"))
	
	#boxplots
	#p = ggplot(data = all_cindices, aes(x = type, y = cindex)) +
    #geom_jitter(alpha = 0.3, color = "tomato") + 
    #geom_boxplot() + labs(title=paste(gene, canc, "\nc-index vs random lncRNAs and clinical variables"))
	##  Add p-value
	#p = p + stat_compare_means(comparisons = my_comparisons) + geom_hline(yintercept=0.5, linetype="dashed", color = "red")
	res = unlist(c(canc, runs, wp_lnc_clinical, diff_meds_lnc_clinical, 
	wp_lnc_clinical_combo, diff_meds_lnc_clinical_combo, 
			lnc_med, clin_med, combo_med))
	#print(p)
	return(res)
}


#pdf("all_lncRNA_cands_vs_random_lncRNAs.pdf")
random_lncs_vs_cand = llply(combos, get_comparison, .progress="text")
#dev.off()
random_lncs_vs_cand1 = ldply(random_lncs_vs_cand)
random_lncs_vs_cand1 = as.data.table(random_lncs_vs_cand1)
colnames(random_lncs_vs_cand1) = c("cancer", "num_rounds", "wp_lnc_clinical", 
	"diff_meds_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical_combo", "median_lncRNA", "median_clinical", "clin_combo_cind_med")
#random_lncs_vs_cand1$wp_lnc_random = as.numeric(random_lncs_vs_cand1$wp_lnc_random)
random_lncs_vs_cand1$wp_lnc_clinical = as.numeric(random_lncs_vs_cand1$wp_lnc_clinical)
random_lncs_vs_cand1$wp_lnc_clinical_combo = as.numeric(random_lncs_vs_cand1$wp_lnc_clinical_combo)

random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical)]
random_lncs_vs_cand1$wp_lnc_clinical_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical, method="fdr")
random_lncs_vs_cand1$wp_lnc_clinical_combo_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical_combo, method="fdr")

write.csv(random_lncs_vs_cand1, file="summary_cross_validation_clinical_multivariate_models_lncRNAs_Nov1.csv", quote=F, row.names=F)

random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical, diff_meds_lnc_clinical, wp_lnc_clinical)]

#try new Figure 2E
#x-axis median difference between lncRNA and clinical variables
#y-axis p-value 
random_lncs_vs_cand1$wp_lnc_clinical_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical, method="fdr")
random_lncs_vs_cand1$wp_lnc_clinical_fdr = -log10(random_lncs_vs_cand1$wp_lnc_clinical_fdr)
random_lncs_vs_cand1$wp_lnc_clinical_combo = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical_combo, method="fdr")

z = which(random_lncs_vs_cand1$wp_lnc_clinical_fdr == "Inf")
random_lncs_vs_cand1$wp_lnc_clinical_fdr[z] = -log10(0.0000001)

#random_lncs_vs_cand1$wp_lnc_random_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_random, method="fdr")
#random_lncs_vs_cand1$wp_lnc_random_fdr = -log10(random_lncs_vs_cand1$wp_lnc_random_fdr)

#plot 1 - lncs vs clinical 
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "cancer"
random_lncs_vs_cand1 = merge(random_lncs_vs_cand1, canc_conv, by="cancer")
random_lncs_vs_cand1$diff_meds_lnc_clinical = as.numeric(random_lncs_vs_cand1$diff_meds_lnc_clinical)
#random_lncs_vs_cand1$diff_meds_lnc_random = as.numeric(random_lncs_vs_cand1$diff_meds_lnc_random)

random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical, diff_meds_lnc_clinical, wp_lnc_clinical)]
random_lncs_vs_cand1$type = factor(random_lncs_vs_cand1$type, levels=unique(random_lncs_vs_cand1$type))

#part a
pdf("figure2_e_lncRNA_cands_vs_clinical_variables_multivariate.pdf", width=6, height=6)
ggplot(random_lncs_vs_cand1, aes(x=diff_meds_lnc_clinical, y=wp_lnc_clinical_fdr)) + geom_point(aes(fill=type), 
       colour="black", pch=21, size=2)+
labs(fill = "Cancer") + 
geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black")+
geom_vline(xintercept=0, linetype="dashed", color = "black")+
labs(x="(median lncRNAs candidate c-index) - (median clinical variables c-index)", y="-log10(adjusted p-val)")+
scale_fill_manual(values = mypal5[1:22]) + theme(legend.position="bottom", text = element_text(size=11))
dev.off()


#-----FIGURE 2E-------------------------------a------------------
#summary of how many individual lncRNAs were selected in each cancer type
#these are the ones that were studied throughout the paper
#show how many of them had a 5% improvement over clinical 
allCands = filter(allCands, data == "TCGA")
head(lncs_perc)

allCands$better_than_clin = ""
z = which(allCands$combo %in% lncs_perc$combo)
allCands$better_than_clin[z] = "yes"
allCands$better_than_clin[-z] = "no"
allCands$multivar_sig = ""

z = which(allCands$fdr <= 0.05)
allCands$multivar_sig[z] = "yes"
allCands$multivar_sig[-z] = "no"
allCands = merge(allCands, canc_conv, by = c("cancer"))

#get summary 
summ = as.data.table(table(allCands$type))#, #allCands$better_than_clin))
#colnames(summ) = c("Cancer", "BetterThanClin", "NumCandidates")
colnames(summ) = c("Cancer", "NumCandidates")
summ = summ[order(-NumCandidates)]

#order of cancer types
od = as.data.table(table(allCands$type))
od = od[order(-N)]

summ$Cancer = factor(summ$Cancer, levels = unique(od$V1))
#summ$BetterThanClin = factor(summ$BetterThanClin, levels = c("yes", "no"))

#barplot summary
# Create the barplot
#pdf("figure2E_aug13.pdf")
pdf("figure2E_oct26.pdf")
ggplot(data=summ, aes(x=Cancer, y=NumCandidates)) +
  theme_bw()+
  geom_bar(stat="identity", color="black")+ 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=c('#E69F00', '#999999')) + 
  xlab("Cancer") + ylab("# of lncRNAs selected in more than 50% of cross validations")
  dev.off()

