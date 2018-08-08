source("source_script_006_166_cands_bootrstapping_on_PCAWG.R")

#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#----------------------------------------Analyze Results---------------------------------------------------------------
#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#results = readRDS(file="lncRNAs_100_external_PCAWG_CVs_individual_cands_june19.rds")
results = readRDS(file="lncRNAs_100_external_PCAWG_CVs_individual_cands_aug8.rds")

results = ldply(results, data.frame)

#pdf("166_EXTERNAL_validation_unique_lncRNAs_cindices_wclinical_june27.pdf", width=9, height=9)
pdf("112_EXTERNAL_validation_unique_lncRNAs_cindices_wclinical_aug8.pdf", width=9, height=9)

results_lncs = as.data.frame(matrix(ncol=2))
colnames(results_lncs) = c("Cancer", "NumHigherthan5percent")

all_meds = as.data.frame(matrix(ncol=4))
colnames(all_meds) = c("lncRNA", "Median", "diff", "Cancer")

for(i in 1:length(unique(results$Cancer))){
	canc_data = subset(results, results$Cancer %in% unique(results$Cancer)[i])
	#x - axis is the type column 
	#y - axis is the C-index
	#facet by lncRNA 
	colnames(canc_data)[3] = "cindex"
	canc_data$cindex = as.numeric(canc_data$cindex)
	
	#If KIRC remove the weird candidate with just one patient risk 
	if(canc_data$Cancer[1] == "Kidney renal clear cell carcinoma"){
		z = which(canc_data$lncRNA == "ENSG00000250360")
		canc_data = canc_data[-z,]
	}

	canc_data =  as.data.table(canc_data)
	all_canc_data = canc_data

	#-------------------------------Clinical VS univaraite lncRNAs------------------------------------------------------
	#-------------------------------------------------------------------------------------------------------------------

	canc_data = filter(canc_data, type %in% c("ClinicalVariables", "lncRNAonly"))

	z = which(str_detect(canc_data$lncRNA, "ENSG"))
	canc_data$lncRNA[-z] = "Clinical"

	z = which(is.na(canc_data$cindex))
	if(!(length(z)==0)){
	canc_data = canc_data[-z,]}

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
	meds$Median = round(meds$Median, digits=2)

	canc_data = merge(canc_data, meds, by="lncRNA")

	#pcawg validated lnc
	z = which(matches$cancer %in% canc_data$Cancer)
	canc_data$PCAWG_sig = ""

	if(!(length(z)==0)){
		lncs_pcawg = matches$gene[z]
		canc_data$PCAWG_sig = ""
		z = which(canc_data$lncRNA %in% lncs_pcawg)
		canc_data$PCAWG_sig[z] = "Yes"
		canc_data$PCAWG_sig[-z] = "No"}

	if(length(z)==0){
		canc_data$PCAWG_sig = "No"
	}

	z = which(canc_data$lncRNA == "Clinical")
	canc_data$PCAWG_sig[z] = "JustClin"

	sig_pcawg = ggplot(canc_data, aes(lncRNA, 1)) +
  	geom_tile(aes(fill = PCAWG_sig)) + theme_bw() + scale_fill_brewer(palette="Pastel1")
  	sig_pcawg = ggpar(sig_pcawg, legend = "right")
	sig_pcawg = sig_pcawg + 
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
	theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

	#meds = unique(canc_data$diff)
	g = ggboxplot(canc_data, x="lncRNA", y="cindex", color="black", fill="type", palette=c("lavenderblush", "aliceblue")) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
	geom_label(data=meds, aes(x=lncRNA ,y = Median, label=Median), col='tomato', size=2)+
	xlab("Predictor")+ theme_bw() + ylim(c(0,1)) +
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") + stat_n_text(size=2)
	
	clin_lncRNA_uni = (ggpar(g, font.tickslab = c(8,"plain", "black"), ylim = c(0, 1), 
 		xtickslab.rt = 45, legend = "right") + ggtitle(canc_data$Cancer[1]))
	
	part1 = clin_lncRNA_uni

	#-------------------------------Clinical VS Combo Clin + lncRNAs----------------------------------------------------
	#-------------------------------------------------------------------------------------------------------------------
	combo_data = filter(all_canc_data, type %in% c("ClinicalVariables", "lncRNA&clin"))
	z = which(str_detect(combo_data$lncRNA, "ENSG"))
	combo_data$lncRNA[-z] = "Clinical"

	z = which(is.na(combo_data$cindex))
	if(!(length(z)==0)){
	combo_data = combo_data[-z,]}

	#keep clinical on the leftmost side and order the lncRNA from lowest to highest median c-index
	z = which(str_detect(combo_data$lncRNA, "Clinical"))
	meds = as.data.table(aggregate(combo_data$cindex[-z], by=list(combo_data$lncRNA[-z]), FUN=median))
	colnames(meds) = c("lncRNA", "Median")
	meds = meds[order(Median)]
	order = c("Clinical", meds$lncRNA)
	combo_data$lncRNA <- factor(combo_data$lncRNA, levels =order)

	#add labels = median increase in c-index 
	clinc_median = median(combo_data$cindex[combo_data$lncRNA %in%"Clinical"])
	#lncRNAs meds
	meds$diff = round(((meds$Median/clinc_median)-1)*100, digits=2)
	row = c("Clinical", clinc_median, round(clinc_median, digits=2))
	names(row) = c("lncRNA", "Median", "diff")
	row = t(as.data.frame(row))
	#meds = as.data.frame(meds)
	meds = rbind(row, meds)
	meds$Median = as.numeric(meds$Median)
	meds$Median = round(meds$Median, digits=2)

	combo_data = merge(combo_data, meds, by="lncRNA")

	#pcawg validated lnc
	z = which(matches$cancer %in% combo_data$Cancer)
	combo_data$PCAWG_sig = ""

	if(!(length(z)==0)){
		lncs_pcawg = matches$gene[z]
		combo_data$PCAWG_sig = ""
		z = which(combo_data$lncRNA %in% lncs_pcawg)
		combo_data$PCAWG_sig[z] = "Yes"
		combo_data$PCAWG_sig[-z] = "No"}

	if(length(z)==0){
		combo_data$PCAWG_sig = "No"
	}

	z = which(combo_data$lncRNA == "Clinical")
	combo_data$PCAWG_sig[z] = "JustClin"

	sig_pcawg_combo = ggplot(combo_data, aes(lncRNA, 1)) +
  	geom_tile(aes(fill = PCAWG_sig)) + theme_bw() + scale_fill_brewer(palette="Pastel1")
  	sig_pcawg_combo = ggpar(sig_pcawg_combo, legend = "none")
	sig_pcawg_combo = sig_pcawg_combo + 
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
	theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

	g = ggboxplot(combo_data, x="lncRNA", y="cindex", color="black", fill="type", palette=c("lavenderblush", "aliceblue")) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
	geom_label(data=meds, aes(x=lncRNA ,y = Median, label=Median), col='tomato', size=2)+
	xlab("Predictor")+ theme_bw() + ylim(c(0,1)) + 
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") + stat_n_text(size=2) 
	
	clin_lncRNA_combo = (ggpar(g, font.tickslab = c(8,"plain", "black"),ylim = c(0, 1), 
 		xtickslab.rt = 45, legend = "right") + ggtitle(combo_data$Cancer[1]))

	part2 = clin_lncRNA_combo

	together = sig_pcawg + part1 + plot_layout(ncol=1, heights=c(1,10)) + sig_pcawg_combo + part2 + plot_layout(ncol=1, heights=c(1,10))
	print(together)

	#SUMMARIZE how many lncRNAs better than clinical 
	#how many lncRNAs with percent increase of c-index greater than 5%?
	canc_data = as.data.table(canc_data)
	canc_data$diff = as.numeric(canc_data$diff) 
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

	


