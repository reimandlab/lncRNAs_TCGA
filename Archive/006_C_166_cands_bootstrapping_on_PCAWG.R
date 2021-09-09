source("source_script_006_166_cands_bootrstapping_on_PCAWG.R")

#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#----------------------------------------Analyze Results---------------------------------------------------------------
#------––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#results = readRDS(file="lncRNAs_100_external_PCAWG_CVs_individual_cands_june19.rds")
results = readRDS(file="lncRNAs_100_external_PCAWG_CVs_individual_cands_aug8.rds")
results = ldply(results, data.frame)

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") 

matches = val_cands
library(EnvStats)


#pdf("166_EXTERNAL_validation_unique_lncRNAs_cindices_wclinical_june27.pdf", width=9, height=9)
pdf("166_EXTERNAL_validation_unique_lncRNAs_cindices_wclinical_aug8.pdf", width=9, height=10)

results_lncs = as.data.frame(matrix(ncol=2))
colnames(results_lncs) = c("Cancer", "NumHigherthan5percent")

all_meds = as.data.frame(matrix(ncol=5))
colnames(all_meds) = c("lncRNA", "Median", "diff", "lnc_name", "Cancer")

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

	print("pass 1")

	#change ENSG to gene names 		
	#for(k in 1:length(unique(canc_data$lncRNA))){
	#	z = which(fantom$CAT_geneID %in% as.character(unique(canc_data$lncRNA)[k]))
	#	if(!(length(z)==0)){
	#		newlnc = fantom$CAT_geneName[z]
	#		z = which(canc_data$lncRNA == unique(canc_data$lncRNA)[k])
	#		canc_data$lncRNA[z] = newlnc
	#	}
	#}

	get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
	}

	canc_data$lnc_name = llply(canc_data$lncRNA, get_name)
	z = which(is.na(canc_data$lnc_name))
	canc_data$lnc_name[z] = "Clinical"

	#keep clinical on the leftmost side and order the lncRNA from lowest to highest median c-index
	z = which(str_detect(canc_data$lncRNA, "Clinical"))
	meds = as.data.table(aggregate(canc_data$cindex[-z], by=list(canc_data$lncRNA[-z]), FUN=median))
	colnames(meds) = c("lncRNA", "Median")
	meds = meds[order(Median)]
	
	meds$lnc_name = llply(meds$lncRNA, get_name)
	z = which(is.na(meds$lnc_name))
	meds$lnc_name[z] = "Clinical"
	meds$lnc_name = as.character(meds$lnc_name)

	order = c("Clinical", meds$lncRNA)
	canc_data$lncRNA <- factor(canc_data$lncRNA, levels =order)
	order = c("Clinical", meds$lnc_name)
	canc_data$lnc_name = factor(canc_data$lnc_name, levels = order)

	#add labels = median increase in c-index 
	clinc_median = median(canc_data$cindex[canc_data$lncRNA %in%"Clinical"])
	#lncRNAs meds
	meds$diff = round(((meds$Median/clinc_median)-1)*100, digits=2)
	row = c("Clinical", clinc_median, "Clinical", round(clinc_median, digits=2))
	names(row) = c("lncRNA", "Median", "lnc_name", "diff")
	row = t(as.data.frame(row))
	#meds = as.data.frame(meds)
	meds = rbind(row, meds)
	meds$Median = as.numeric(meds$Median)
	meds$Median = round(meds$Median, digits=2)

	canc_data = merge(canc_data, meds, by=c("lncRNA", "lnc_name"))
	
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

	sig_pcawg = ggplot(canc_data, aes(lnc_name, 1)) +
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
	g = ggboxplot(canc_data, x="lnc_name", y="cindex", color="black", fill="type", palette=c("lightpink3", "lightsteelblue4")) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
	geom_label(data=meds, aes(x=lnc_name ,y = Median, label=Median), col='tomato', size=1.4)+
	xlab("Predictor")+ theme_bw() + ylim(c(0,1)) +
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") + stat_n_text(size=2)
	
	clin_lncRNA_uni = (ggpar(g, font.tickslab = c(8,"plain", "black"), ylim = c(0, 1), 
 		xtickslab.rt = 45, legend = "right") + ggtitle(canc_data$Cancer[1]))
	
	part1 = clin_lncRNA_uni
	print("pass part1")

	#-------------------------------Clinical VS Combo Clin + lncRNAs----------------------------------------------------
	#-------------------------------------------------------------------------------------------------------------------
	combo_data = filter(all_canc_data, type %in% c("ClinicalVariables", "lncRNA&clin"))
	z = which(str_detect(combo_data$lncRNA, "ENSG"))
	combo_data$lncRNA[-z] = "Clinical"

	z = which(is.na(combo_data$cindex))
	if(!(length(z)==0)){
	combo_data = combo_data[-z,]}

	get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
	}

	combo_data$lnc_name = llply(combo_data$lncRNA, get_name)
	z = which(is.na(combo_data$lnc_name))
	combo_data$lnc_name[z] = "Clinical"

	#keep clinical on the leftmost side and order the lncRNA from lowest to highest median c-index
	z = which(str_detect(combo_data$lncRNA, "Clinical"))
	meds = as.data.table(aggregate(combo_data$cindex[-z], by=list(combo_data$lncRNA[-z]), FUN=median))
	colnames(meds) = c("lncRNA", "Median")
	meds$lnc_name = llply(meds$lncRNA, get_name)
	z = which(is.na(meds$lnc_name))
	meds$lnc_name[z] = "Clinical"
	meds$lnc_name = as.character(meds$lnc_name)
	meds = meds[order(Median)]

	order = c("Clinical", meds$lncRNA)
	combo_data$lncRNA <- factor(combo_data$lncRNA, levels =order)
	order = c("Clinical", meds$lnc_name)
	combo_data$lnc_name <- factor(combo_data$lnc_name, levels =order)

	#add labels = median increase in c-index 
	clinc_median = median(combo_data$cindex[combo_data$lncRNA %in%"Clinical"])
	#lncRNAs meds
	meds$diff = round(((meds$Median/clinc_median)-1)*100, digits=2)
	row = c("Clinical", clinc_median, "Clinical", round(clinc_median, digits=2))
	names(row) = c("lncRNA", "Median", "lnc_name", "diff")
	row = t(as.data.frame(row))
	#meds = as.data.frame(meds)
	meds = rbind(row, meds)
	meds$Median = as.numeric(meds$Median)
	meds$Median = round(meds$Median, digits=2)

	combo_data = merge(combo_data, meds, by=c("lncRNA", "lnc_name"))

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

	sig_pcawg_combo = ggplot(combo_data, aes(lnc_name, 1)) +
  	geom_tile(aes(fill = PCAWG_sig)) + theme_bw() + scale_fill_brewer(palette="Pastel1")
  	sig_pcawg_combo = ggpar(sig_pcawg_combo, legend = "none")
	sig_pcawg_combo = sig_pcawg_combo + 
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
	theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

	g = ggboxplot(combo_data, x="lnc_name", y="cindex", color="black", fill="type", palette=c("lightpink3", "lightsteelblue4")) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
	geom_label(data=meds, aes(x=lnc_name ,y = Median, label=Median), col='tomato', size=1.4)+
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

head(all_meds)
all_meds = all_meds[-1,]
all_meds = as.data.table(all_meds)
all_meds$diff = as.numeric(all_meds$diff)
all_meds = filter(all_meds, Median >=0.6)
all_meds = filter(all_meds, diff >=5)
all_meds = as.data.table(all_meds)
all_meds$pcawg_val = ""
all_meds$lnc_name[which(all_meds$lncRNA %in% val_cands$gene)] = paste("*", all_meds$lnc_name, "*")
all_meds = as.data.table(all_meds)
filter(all_meds, !(pcawg_val == "pcawg_val"))
all_meds$combo = paste(all_meds$lncRNA, all_meds$Cancer, sep="_")

z = which((allCands$combo %in% all_meds$combo) & (allCands$data == "TCGA"))

#remove KIRP for now - something weird is going on with it 
#need to double check it 
z = which(all_meds$Cancer == "Kidney renal papillary cell carcinoma")
all_meds = all_meds[-z,]
all_meds$type = ""
for(i in 1:nrow(all_meds)){
	canc = all_meds$Cancer[i]
	type = rna$type[which(rna$Cancer==canc)][1]
	all_meds$type[i] = type
} 

#geom_tile
#x = cancer
#y = lncRNA 
#fill = diff
#shape = Median
#text = pcawg pval
all_meds = all_meds[order(type, Median)]
all_meds$type = factor(all_meds$type, levels=unique(all_meds$type))
all_meds$lnc_name = factor(all_meds$lnc_name, levels=unique(all_meds$lnc_name))

pdf("summary_best_external_val_PCAWG.pdf")
ggplot(all_meds, aes(type, lnc_name)) +
  geom_tile(aes(fill = diff)) + geom_point(aes(size = Median))+ 
    scale_fill_gradient(low = "darkcyan", high = "orange", na.value = 'transparent') +
    xlab("Cancer") + ylab("lncRNA") + theme_bw()
ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend.title="Median C-index \ndifference")

dev.off()










