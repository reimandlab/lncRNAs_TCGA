set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_name")]
allCands = allCands[!duplicated(allCands), ]

#-----------------------------------------------------------------------------------------------------------------------------------
#1 - read in data  
#-----------------------------------------------------------------------------------------------------------------------------------

multivariate = as.data.table(ldply(readRDS("lncRNAs_clinical_variables_risk_score_models_multivariate.rds")))
univariate = as.data.table(ldply(readRDS("lncRNAs_clinical_variables_risk_score_models_univariate.rds")))

#-----------------------------------------------------------------------------------------------------------------------------------
#2 - prepare data for plot  
#-----------------------------------------------------------------------------------------------------------------------------------

cancers = unique(univariate$Cancer)

get_multi_plot = function(canc){

	multi_dat = as.data.table(filter(multivariate, Cancer == canc))
	single_dat = as.data.table(filter(univariate, Cancer == canc, !(type=="clinical_variables")))

	#make input for plot 
	multi_dat$analysis="multi"
	single_dat$analysis="single"
	
	colnames(multi_dat)[3] = "cindex"
	colnames(single_dat)[3] = "cindex"

	all = rbind(multi_dat, single_dat)
	all = as.data.table(filter(all, !(is.na(cindex))))

	z = which(str_detect(all$lncRNA, "ENSG"))
	all$lncRNA[z] = sapply(all$lncRNA[z], get_name)

	all$lncRNA[all$type=="ClinicalVariables"] = "clin"
	all$type[all$type=="ClinicalVariables"] = "clinical"

	all$label =  paste(all$lncRNA, all$type)
	all$label[all$label=="all_lncs lncRNAonly"] = "all lncRNAs"
	all$label[all$label=="clin clinical"] = "clinical variables"
	all$label[all$label=="all_lncs lncRNA&clin"] = "all lncRNAs + clinical"

	#which is the best model?
	all = all[order(cindex)]
	all$label = factor(all$label, levels=all$label)
	all$type[all$lncRNA == "all_lncs"] = "all_lncs"
	all$type[all$label == "all lncRNAs + clinical"] = "all_lncs_clin"

	#within each category choose best model
	choose = as.data.table(filter(all, type %in% c("lncRNA&clin" , "lncRNAonly")))
	not_choose = as.data.table(filter(all, !(type %in% c("lncRNA&clin" , "lncRNAonly"))))

	types = c("lncRNA&clin" , "lncRNAonly")
	get_best_type = function(t){
		dat = as.data.table(filter(choose, type == t))
		dat = dat[order(-cindex)][1,]
		return(dat)
	}

	choose_dat = as.data.table(ldply(llply(types, get_best_type)))
	all = rbind(choose_dat, not_choose)
	all$HR = as.numeric(all$HR)
	all$ci05 = as.numeric(all$ci05)	
	all$ci95 = as.numeric(all$ci95)
	all$pval = as.numeric(all$pval)
	all$HR = log2(all$HR)
	all$ci05 = log2(all$ci05)
	all$ci95 = log2(all$ci95)
	all=all[order(-type)]
	all$label = factor(all$label, levels=all$label)
	all$type_canc = canc_conv$type[canc_conv$Cancer == all$Cancer[1]]
	all$cindex=as.numeric(all$cindex)

	# Change error plot type and add mean points
	g = ggplot(all, aes(x=label, y=cindex, colour=type, group=type)) + 
    geom_point(size=3)+ theme(plot.title = element_text(hjust = 0.5))+ theme_bw()+ ggtitle(all$type_canc[1])
	g=ggpar(g, font.xtickslab=c(5, "plain", "black"), font.ytickslab=c(5, "plain", "black")) + 
	ylim(c(0,1))+
    coord_flip()+ylab("c-index")+rremove("legend")+
	geom_hline(yintercept=0.5, linetype="dashed", color = "black")

    #geom_errorbar(aes(ymin=ci05, ymax=ci95), colour="black", width=.1) +
    return(g)

}

list_plots = llply(cancers, get_multi_plot)
library(gridExtra)
n <- length(list_plots)
nCol <- floor(sqrt(n))

pdf("/u/kisaev/cancers_multivaraite_models_via_risk_scores.pdf", width=11, height=9)
do.call("grid.arrange", c(list_plots, ncol=nCol))
dev.off()

