set.seed(911)

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

#load all libraries and functions 
source("check_lnc_exp_cancers.R")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#----------------------------------------Analyze Results---------------------------------------------------------------
#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

results = readRDS(file="lncRNAs_100_internal_CVs_individual_cands_june19.rds")
results = do.call(rbind.data.frame, results)

pdf("166_unique_lncRNAs_cindices.pdf", width=10, height=6)

results_lncs = as.data.frame(matrix(ncol=2))
colnames(results_lncs) = c("Cancer", "NumHigherthan5percent")

all_meds = as.data.frame(matrix(ncol=4))
colnames(all_meds) = c("lncRNA", "Median", "diff", "Cancer")

cisums_results = as.data.frame(matrix(ncol=9))
colnames(cisums_results) = c("Cancer", "lncRNA", "type", "n", "sd", 
	"mean", "median", "low", "high")

all_canc_data_robust = as.data.frame(matrix(ncol=7))
colnames(all_canc_data_robust) = c("lncRNA", "Cancer", "cindex", "type", "round", "Median", "diff")

for(i in 1:length(unique(results$Cancer))){
	canc_data = subset(results, results$Cancer %in% unique(results$Cancer)[i])

	#x - axis is the type column 
	#y - axis is the C-index
	#facet by lncRNA 
	colnames(canc_data)[3] = "cindex"
	canc_data$cindex = as.numeric(canc_data$cindex)
	canc_data =  as.data.table(canc_data)
	canc_data = filter(canc_data, type %in% c("ClinicalVariables", "lncRNAonly"))

	z = which(str_detect(canc_data$lncRNA, "ENSG"))
	canc_data$lncRNA[-z] = "Clinical"

	z = which(is.na(canc_data$cindex))
	if(!(length(z)==0)){
	canc_data = canc_data[-z,]}

	#first get confidence intervals 
	#need to calcualte confidence interval 
	cisum = as.data.table(canc_data %>% dplyr::group_by(Cancer, lncRNA, type) %>% 
	  dplyr::summarize(n = n(), sd = sd(cindex), mean = mean(cindex), median=median(cindex), low=quantile(cindex, 0.025), high=quantile(cindex, 0.975)))
	cisum = cisum[order(-median, sd)]
	#cisum = as.data.table(dplyr::filter(cisum, type == "lncRNAonly", low >=0.5, high>=0.5))
	cisum = as.data.table(dplyr::filter(cisum, type == "lncRNAonly"))
	as.data.table(dplyr::filter(cisum, median >= 0.5))

	cisums_results = rbind(cisums_results, cisum)

	#keep only those lncRNAs with CIs not overlapping 0.5 
	#z = which(canc_data$lncRNA %in% c(cisum$lncRNA, "Clinical"))
	#canc_data = canc_data[z,]

	#change ENSG to gene names 
	#for(k in 1:length(unique(canc_data$lncRNA))){
	#	z = which(fantom$CAT_geneID %in% as.character(unique(canc_data$lncRNA)[k]))
	#	if(!(length(z)==0)){
	#		newlnc = fantom$CAT_geneName[z]
	#		z = which(canc_data$lncRNA == unique(canc_data$lncRNA)[k])
	#		canc_data$lncRNA[z] = newlnc
	#	}
	#}

	#keep clinical on the leftmost side and order the lncRNA from lowest to highest median c-index
	z = which(str_detect(canc_data$lncRNA, "Clinical"))
	meds = as.data.table(aggregate(canc_data$cindex[-z], by=list(canc_data$lncRNA[-z]), FUN=median))
	colnames(meds) = c("lncRNA", "Median")
	meds = meds[order(Median)]
	order = c("Clinical", meds$lncRNA)
	canc_data$lncRNA <- factor(canc_data$lncRNA, levels =order)

	#add labels = median increase in c-index 
	clinc_median = median(canc_data$cindex[canc_data$lncRNA %in%"Clinical"])
	#lncRNAs meds
	meds$diff = round(((meds$Median/clinc_median)-1)*100, digits=2)
	row = c("Clinical", clinc_median, round(clinc_median, digits=2))
	names(row) = c("lncRNA", "Median", "diff")
	row = t(as.data.frame(row))
	#meds = as.data.frame(meds)
	meds = rbind(row, meds)
	meds$Median = as.numeric(meds$Median)

	canc_data = merge(canc_data, meds, by="lncRNA")
	#meds = unique(canc_data$diff)
	g = ggboxplot(canc_data, x="lncRNA", y="cindex", color="black", fill="type", palette=c("lightpink3", "lightsteelblue4")) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
	geom_label(data=meds, aes(x=lncRNA ,y = Median, label=Median), col='tomato', size=2)+
	xlab("Predictor")+ theme_classic() + ylim(c(0,1)) + 
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") + 
	ylim(0,1)
	
	print(ggpar(g, font.tickslab = c(8,"plain", "black"),
 		xtickslab.rt = 45) + ggtitle(canc_data$Cancer[1]))
	
	#SUMMARIZE how many lncRNAs better than clinical 
	#how many lncRNAs with percent increase of c-index greater than 5%?
	canc_data = as.data.table(canc_data)
	canc_data$diff = as.numeric(canc_data$diff) 
	
	all_canc_data_robust = rbind(all_canc_data_robust, canc_data)

	higher = filter(canc_data, diff >= 5)

	meds = filter(meds, !(lncRNA %in% "Clinical"))
	meds$Cancer = canc_data$Cancer[1]
	all_meds = rbind(all_meds, meds)

	#num 5% higher 
	row = c(canc_data$Cancer[1], length(unique(higher$lncRNA)))
	names(row) = colnames(results_lncs)
	results_lncs = rbind(results_lncs, row)
}

dev.off()

#final table of robust lncRNAs 
cisums_results = cisums_results[-1,]
cisums_results$combo = paste(cisums_results$lncRNA, cisums_results$Cancer, sep="_")

#112 unique combos that are robust 
#109 unique lncRNAs 
saveRDS(cisums_results, file="112_combos_robust_internal_validation_survival_lncRNAs_aug8.rds")

#all cindices for all robust lncRNAs
all_canc_data_robust = all_canc_data_robust[-1,]
saveRDS(all_canc_data_robust, file="112_cindicies_combos_robust_internal_validation_survival_lncRNAs_aug8.rds")

#summary
results_lncs = results_lncs[-1,]
results_lncs = as.data.table(results_lncs)
results_lncs = filter(results_lncs, NumHigherthan5percent > 0 )
results_lncs = as.data.table(results_lncs)
results_lncs$NumHigherthan5percent = as.numeric(results_lncs$NumHigherthan5percent)
results_lncs = results_lncs[order(NumHigherthan5percent)]

#how many of the lncRNAs with greater c-index also remain significant when fitting multivariate models
multi = readRDS("TCGA_results_multivariate_results_June22.rds")
multi = as.data.table(multi)
multi = filter(multi, fdr_pval <0.05)

head(all_meds)
all_meds = all_meds[-1,]
all_meds = as.data.table(all_meds)
all_meds$diff = as.numeric(all_meds$diff)
all_meds = filter(all_meds, diff >=5)

#how many of these significant in multivaraite models
all_meds = as.data.table(all_meds)
all_meds = filter(all_meds, lncRNA %in% multi$gene)
all_meds$combo = paste(all_meds$lncRNA, all_meds$Cancer, sep="_")
#merge confidence interval and % median c-index increase 
all_meds = merge(all_meds, cisums_results, by=c("lncRNA", "Cancer", "combo")
all_meds = as.data.table(all_meds)
saveRDS(all_meds, file="148_combos_robust_5perc_increase_internal_validation_survival_lncRNAs_aug9.rds")









