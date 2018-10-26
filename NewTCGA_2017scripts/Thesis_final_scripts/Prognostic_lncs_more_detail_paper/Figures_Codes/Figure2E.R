#------------------------------------------------------------------------------
#This script takes all the PCGs that are significantly differentially expressed
#between lncRNA risk groups and calculates their correlation 
#plots heatmap of this data as well as lncRNA and pcg prognostic
#relationship 
#------------------------------------------------------------------------------

#source code
source("check_lnc_exp_cancers.R")

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#canc conversion
canc_conv = rna[,c(which(colnames(rna) %in% c("Cancer", "type")))]
canc_conv = canc_conv[!duplicated(canc_conv),]
colnames(canc_conv) = c("type", "cancer")

#------DATA-----------------------------------------------------

#lncRNAs with 5% improvemenet over clinical 
lncs_perc = readRDS("173_combos_robust_5perc_increase_internal_validation_survival_lncRNAs_aug9.rds")
lncs_perc$better = ""
z = which(lncs_perc$diff >= 0.1)

################################################################################
#cindices from 166 candidates 
lnc_cands = readRDS("lncRNAs_100_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_cands)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[-z,]
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
z= which(str_detect(lists$lncRNA, "ENSG"))
lncs_cands_lists = lists[z,]
#---------------------------------------------------------------------------------

################################################################################
#cindices from clinical variables
lnc_rand = readRDS("lncRNAs_RANDOM_100_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_rand)
r = as.data.table(r)
z = which(r$type == "ClinicalVariables")
r = r[z,]

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
#---------------------------------------------------------------------------------



################################################################################
#cindices from 100 random lncRNAs 
lnc_rand = readRDS("lncRNAs_RANDOM_100_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_rand)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[-z,]
list_r = split(r, by = "Cancer")
lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
z= which(str_detect(lists$lncRNA, "ENSG"))
lncs_random_lists = lists[z,]
#---------------------------------------------------------------------------------



################################################################################
#cindices from 1000 random PCGs
pcg_rand = readRDS("PCGs__RANDOM_100_internal_CVs_individual_cands_june19.rds")
r = ldply(pcg_rand)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[-z,]
list_r = split(r, by = "Cancer")
lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
z= which(str_detect(lists$lncRNA, "ENSG"))
pcgs_random_lists = lists[z,]
#---------------------------------------------------------------------------------


#for each cancer type 
#plot summary clin vs lncRNA cands vs random lncs vs random PCGs

full_summary_cinds = function(canc){

	#[1] lncs cands 
	lncs_cands = as.data.table(filter(lncs_cands_lists, .id==canc))
	if(length(unique(lncs_cands$lncRNA)) >1){
		lncs_cands$overall_mean = mean(lncs_cands$mean_cindex)
		lncs_cands$overall_sd = sd(lncs_cands$mean_cindex)
	}
	if(!(length(unique(lncs_cands$lncRNA)) >1)){
		lncs_cands$overall_mean = (lncs_cands$mean_cindex)
		lncs_cands$overall_sd = (lncs_cands$sd)
	}

	lncs_cands$type = "lncRNA_cands"

	#[2] clin variables
	clin_cands = as.data.table(filter(clinical_lists, .id==canc))
	if(length(unique(clin_cands$lncRNA)) >1){
		clin_cands$overall_mean = mean(clin_cands$mean_cindex)
		clin_cands$overall_sd = sd(clin_cands$mean_cindex)
	}
	if(!(length(unique(clin_cands$lncRNA)) >1)){
		clin_cands$overall_mean = (clin_cands$mean_cindex)
		clin_cands$overall_sd = (clin_cands$sd)
	}

	clin_cands$type = "clin_cands"

	#[3] random lncs
	rand_lncs = as.data.table(filter(lncs_random_lists, .id==canc))
	if(length(unique(rand_lncs$lncRNA)) >1){
		rand_lncs$overall_mean = mean(rand_lncs$mean_cindex)
		rand_lncs$overall_sd = sd(rand_lncs$mean_cindex)
	}
	if(!(length(unique(rand_lncs$lncRNA)) >1)){
		rand_lncs$overall_mean = (rand_lncs$mean_cindex)
		rand_lncs$overall_sd = (rand_lncs$sd)
	}

	rand_lncs$type = "rand_lncs"

	#[4] random PCGs
	rand_pcgs = as.data.table(filter(pcgs_random_lists, .id==canc))
	if(length(unique(rand_pcgs$lncRNA)) >1){
		rand_pcgs$overall_mean = mean(rand_pcgs$mean_cindex)
		rand_pcgs$overall_sd = sd(rand_pcgs$mean_cindex)
	}
	if(!(length(unique(rand_pcgs$lncRNA)) >1)){
		rand_pcgs$overall_mean = (rand_pcgs$mean_cindex)
		rand_pcgs$overall_sd = (rand_pcgs$sd)
	}

	rand_pcgs$type = "rand_pcgs"

	#all variables
	all_cs = rbind(lncs_cands, clin_cands, rand_lncs, rand_pcgs)

	#for plotting
	all_cs = all_cs[,c("type", "overall_mean", "overall_sd")]
	all_cs = unique(all_cs)

	#make plot
	# Make the graph with the 95% confidence interval
	p = ggplot(all_cs, aes(x=type, y=overall_mean, group=1)) +
    geom_line() +
    geom_errorbar(width=.1, aes(ymin=overall_mean-overall_sd, ymax=overall_mean+overall_sd)) +
    geom_point(shape=21, size=3, fill="white")+
    ylim(0,1)
    return(p)

}


#RUN 
cancs = as.list(unique(allCands$cancer))
plots_cancs = llply(cancs, full_summary_cinds)
plots_cancs
dev.off()


#-----FIGURE 2E-------------------------------------------------
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

