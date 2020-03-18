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
	
	multi_dat = as.data.table(filter(multivariate_lncs, Cancer == canc))
	single_dat = as.data.table(filter(single_lncs, Cancer == canc, !(type=="clinical_variables")))

	#make input for plot 
	multi_dat$analysis="multi"
	single_dat$analysis="single"
	
	multi_dat = multi_dat[,c("lncRNA",  "Cancer",  "cindex", "type", "round", "analysis", "wp_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical", "diff_meds_lnc_clinical_combo" , 
		"lnc_med", "clin_med", "combo_med")]
	single_dat = single_dat[,c("lncRNA",  "Cancer",  "cindex", "type", "round", "analysis", "wp_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical", "diff_meds_lnc_clinical_combo", 
		"lnc_med", "clin_med", "combo_med")]

	all = rbind(multi_dat, single_dat)
	all = as.data.table(filter(all, !(is.na(cindex))))

	z = which(str_detect(all$lncRNA, "ENSG"))
	all$lncRNA[z] = sapply(all$lncRNA[z], get_name)

	all$lncRNA[all$type=="clinical_variables"] = "clin"
	all$type[all$type=="clinical_variables"] = "clinical"
	all$type[all$type=="lncRNA_canc"] = "lncRNA"

	all$label =  paste(all$lncRNA, all$type)
	all$label[all$label=="all_lncs lncRNA"] = "all lncRNAs"
	all$label[all$label=="clin clinical"] = "clinical variables"
	all$label[all$label=="all_lncs lncRNA&clin"] = "all lncRNAs + clinical"

	#full lnc model vs full lnc model + clin 
	full_lnc_clin = tidy(wilcox.test(all$cindex[all$label =="all lncRNAs + clinical"], all$cindex[all$label == "all lncRNAs"]))[2]
	full_lnc_clin_diff = median(all$cindex[all$label =="all lncRNAs"]) - median(all$cindex[all$label == "all lncRNAs + clinical"])

	all$full_lnc_clin_pval = full_lnc_clin
	all$full_lnc_clin_diff = full_lnc_clin_diff
	
	#get median c-index by labl
	meds = as.data.table(all %>% 
		group_by(label) %>% summarise(median_cindex = median(cindex)))

	#test individual lncRNAs that are better than model with all lncRNAs s
	lncs_test = filter(meds, median_cindex > meds$median_cindex[meds$label == "all lncRNAs"])$label
	if(!(length(lncs_test)==0)){
		pvals = unlist(sapply(lncs_test, function(x){tidy(wilcox.test(all$cindex[which(all$label == x)], all$cindex[all$label == "all lncRNAs"]))[2]}))
		med_diffs = unlist(sapply(lncs_test, function(x){median(all$cindex[which(all$label == x)]) - median(all$cindex[all$label == "all lncRNAs"])}))

		names(lncs_test) = lncs_test
		names(pvals) = lncs_test
		names(med_diffs) = lncs_test

		super_lncs = rbind(lncs_test, pvals, med_diffs)
		super_lncs=super_lncs[-1,]
		super_lncs=as.data.frame(super_lncs)
		super_lncs$cancer = all$Cancer[1]
		print(super_lncs)
	}

	#which is the best model?
	meds = meds[order(median_cindex)]
	all$label = factor(all$label, levels=meds$label)


	# Change error plot type and add mean points
	g = ggerrorplot(all, x = "label", y = "cindex", size=1, 
            desc_stat = "median_mad", palette="npg", 
            error.plot = "errorbar",            # Change error plot type
            add = "median" , color="type" ,                   # Add mean points
            ) + theme_bw() + ggtitle(all$Cancer[1])
	g=ggpar(g, font.xtickslab=c(7, "plain", "black"), font.ytickslab=c(6, "plain", "black")) + ylim(c(0,1))+
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") + stat_compare_means(label = "p.signif", 
                     ref.group = "clinical variables", hide.ns = TRUE)  + coord_flip()
	print(g)

	meds = meds[order(-median_cindex)]

	all$label = factor(all$label, levels=meds$label)
	all$cindex=NULL
	all$round=NULL
	all = unique(all)
	all$diff_meds_lnc_clinical_combo = all$combo_med - all$lnc_med

	return(all)

}

pdf("/u/kisaev/cancers_multivaraite_models.pdf", width=6, height=5)
res=as.data.table(ldply(llply(cancers, get_bar)))
dev.off()
res$full_lnc_clin_fdr = p.adjust(res$full_lnc_clin_pval, method="fdr")
write.csv(res, file="/u/kisaev/all_lncRNA_cands_single_multi_models_performance.csv", quote=F, row.names=F)

#for each cancer, extract model with highest  c-index
#get risk score for each patient 
#use continuous risk score as predictor
#split patients into KM plot by  score 

get_super_lncs = function(canc){
	
	multi_dat = as.data.table(filter(multivariate_lncs, Cancer == canc))
	single_dat = as.data.table(filter(single_lncs, Cancer == canc, !(type=="clinical_variables")))

	#make input for plot 
	multi_dat$analysis="multi"
	single_dat$analysis="single"
	
	multi_dat = multi_dat[,c("lncRNA",  "Cancer",  "cindex", "type", "round", "analysis", "wp_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical", "diff_meds_lnc_clinical_combo" , 
		"lnc_med", "clin_med", "combo_med")]
	single_dat = single_dat[,c("lncRNA",  "Cancer",  "cindex", "type", "round", "analysis", "wp_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical", "diff_meds_lnc_clinical_combo", 
		"lnc_med", "clin_med", "combo_med")]

	all = rbind(multi_dat, single_dat)
	all = as.data.table(filter(all, !(is.na(cindex))))

	z = which(str_detect(all$lncRNA, "ENSG"))
	all$lncRNA[z] = sapply(all$lncRNA[z], get_name)

	all$lncRNA[all$type=="clinical_variables"] = "clin"
	all$type[all$type=="clinical_variables"] = "clinical"
	all$type[all$type=="lncRNA_canc"] = "lncRNA"

	all$label =  paste(all$lncRNA, all$type)
	all$label[all$label=="all_lncs lncRNA"] = "all lncRNAs"
	all$label[all$label=="clin clinical"] = "clinical variables"
	all$label[all$label=="all_lncs lncRNA&clin"] = "all lncRNAs + clinical"

	#full lnc model vs full lnc model + clin 
	full_lnc_clin = tidy(wilcox.test(all$cindex[all$label =="all lncRNAs + clinical"], all$cindex[all$label == "all lncRNAs"]))[2]
	full_lnc_clin_diff = median(all$cindex[all$label =="all lncRNAs"]) - median(all$cindex[all$label == "all lncRNAs + clinical"])

	all$full_lnc_clin_pval = full_lnc_clin
	all$full_lnc_clin_diff = full_lnc_clin_diff
	
	#get median c-index by labl
	meds = as.data.table(all %>% 
		group_by(label) %>% summarise(median_cindex = median(cindex)))

	#test individual lncRNAs that are better than model with all lncRNAs s
	lncs_test = filter(meds, median_cindex > meds$median_cindex[meds$label == "all lncRNAs"])$label
	if(!(length(lncs_test)==0)){
		pvals = unlist(sapply(lncs_test, function(x){tidy(wilcox.test(all$cindex[which(all$label == x)], all$cindex[all$label == "all lncRNAs"]))[2]}))
		med_diffs = unlist(sapply(lncs_test, function(x){median(all$cindex[which(all$label == x)]) - median(all$cindex[all$label == "all lncRNAs"])}))

		names(lncs_test) = lncs_test
		names(pvals) = lncs_test
		names(med_diffs) = lncs_test

		super_lncs=as.data.frame(matrix(ncol=3, nrow=length(lncs_test)))
		colnames(super_lncs)=c("lnc", "pval", "med_diff")
		super_lncs$lnc = lncs_test
		super_lncs$pval = pvals
		super_lncs$med_diff = med_diffs
		super_lncs$cancer = all$Cancer[1]
	}

	return(super_lncs)

}

super_lncs = as.data.table(ldply(llply(cancers, get_super_lncs)))
super_lncs = as.data.table(filter(super_lncs, !(lnc == "all lncRNAs + clinical")))
z=which(duplicated(super_lncs$lnc))
super_lncs = super_lncs[-z,]
super_lncs$fdr = p.adjust(super_lncs$pval, method="fdr")

fig = unique(res[,c("Cancer", "full_lnc_clin_fdr", "full_lnc_clin_diff")])
fig=merge(fig, canc_conv)

pdf("/u/kisaev/figure2X_multivariate_multivariate_plus_clin.pdf", width=6, height=6)
ggplot(data=fig, aes(x=full_lnc_clin_diff, y=-log10(full_lnc_clin_fdr), label=type)) + geom_point( 
       color="red", size=2) + geom_text_repel()+
geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black")+
geom_vline(xintercept=0, linetype="dashed", color = "black")+ theme_bw()+
labs(x="(median cindex multivariate lncRNA model) - (median cindex lncRNAs + clinical)", y="-log10(adjusted p-val)")+
theme(legend.position="bottom", text = element_text(size=11))
dev.off()

#summary of super lncRNAs 
figu2 = as.data.table(filter(res, label %in% super_lncs$lnc))
colnames(figu2)[12] = "lnc"
figu2=merge(super_lncs, figu2, by="lnc")
figu2=merge(figu2, canc_conv, by="Cancer")
figu2$name = paste(figu2$lncRNA, figu2$type.y)

pdf("/u/kisaev/figure2XX_multivariate_vs_lnc_only.pdf", width=6, height=6)
figu2$plot = ""
figu2$plot[figu2$med_diff > 0.05] = figu2$name[figu2$med_diff > 0.05]
write.csv(figu2, file="/u/kisaev/input_figure2XX_multivariate_vs_lnc_only.csv", quote=F, row.names=F)
ggplot(data=figu2, aes(x=med_diff, y=-log10(fdr), label=plot)) + geom_point( 
       shape = 21, colour = "black", fill = "red", size=2) + geom_text_repel()+
geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black")+
geom_vline(xintercept=0, linetype="dashed", color = "black")+ theme_bw()+
labs(x="(median cindex lncRNA model) - (median cindex multivariate lncRNAs)", y="-log10(adjusted p-val)")+
theme(legend.position="bottom", text = element_text(size=11))
dev.off()





