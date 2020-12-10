set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs

#which cancer types are the non-unique lncRNAs from?
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_name")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#------DATA-----------------------------------------------------

################################################################################
#cindices from 166 candidates
lnc_cands = readRDS("lncRNAs_1000_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_cands)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[-z,]

#remove those that are no longer considred in our analysis
r$combo = paste(r$lncRNA, r$Cancer, sep="_")
z1 = which(r$combo %in% allCands$combo)
z2 = which(!(str_detect(r$lncRNA, "ENSG")))

r = r[c(z1,z2),]

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
lncs_cands_all_cindices = r

#---------------------------------------------------------------------------------

################################################################################
#cindices from clinical variables
lnc_rand = readRDS("lncRNAs_1000_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_rand)
r = as.data.table(r)
z = which(r$type == "ClinicalVariables")
r = r[z,]

#remove those that are no longer considred in our analysis
r$combo = paste(r$lncRNA, r$Cancer, sep="_")
z1 = which(r$combo %in% allCands$combo)
z2 = which(!(str_detect(r$lncRNA, "ENSG")))

r = r[c(z1,z2),]

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

#---------------------------------------------------------------------------------



################################################################################
#cindices from 100 random lncRNAs
#lnc_rand = readRDS("lncRNAs_RANDOM_100_internal_CVs_individual_cands_june19.rds") #<- still need to re-do this analysis
#r = ldply(lnc_rand)
#r = as.data.table(r)
#z = which(r$type == "lncRNA&clin")
#r = r[-z,]

#list_r = split(r, by = "Cancer")
#lists = llply(list_r, get_sum)
#lists = ldply(lists)
#lists = as.data.table(lists)
#z= which(str_detect(lists$lncRNA, "ENSG"))
#lncs_random_lists = lists[z,]
#lncs_random_cindices = r
#---------------------------------------------------------------------------------


################################################################################
#cindices from 100 combined models lncRNA + clinical
lnc_rand = readRDS("lncRNAs_1000_internal_CVs_individual_cands_june19.rds")
r = ldply(lnc_rand)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[z,]

#remove those that are no longer considred in our analysis
r$combo = paste(r$lncRNA, r$Cancer, sep="_")
z1 = which(r$combo %in% allCands$combo)
z2 = which(!(str_detect(r$lncRNA, "ENSG")))

r = r[c(z1,z2),]

list_r = split(r, by = "Cancer")
lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
z= which(str_detect(lists$lncRNA, "ENSG"))
lncsclinccombo_random_lists = lists[z,]
lncsclinccombo_random_cindices = r
#---------------------------------------------------------------------------------


################################################################################
#cindices from 1000 random PCGs
#pcg_rand = readRDS("PCGs__RANDOM_100_internal_CVs_individual_cands_june19.rds")
#r = ldply(pcg_rand)
#r = as.data.table(r)
#z = which(r$type == "lncRNA&clin")
#r = r[-z,]
#list_r = split(r, by = "Cancer")
#lists = llply(list_r, get_sum)
#lists = ldply(lists)
#lists = as.data.table(lists)
#z= which(str_detect(lists$lncRNA, "ENSG"))
#pcgs_random_lists = lists[z,]
#pcgs_random_cindices = r
#---------------------------------------------------------------------------------

combos = (unique(allCands$combo))
lncs_cands_all_cindices$combo = paste(lncs_cands_all_cindices$lncRNA, lncs_cands_all_cindices$Cancer, sep="_")
lncsclinccombo_random_cindices$combo = paste(lncsclinccombo_random_cindices$lncRNA, lncsclinccombo_random_cindices$Cancer, sep="_")

get_comparison = function(comboo){

	#get lncRNA candidate c-index
	lnc_cind = as.data.table(filter(lncs_cands_all_cindices, combo == comboo))

	#how many runs works
	runs = dim(lnc_cind)
	lnc_cind$type = "lncRNA_canc"

	#get random lncs
	#make sure they dont include any of the candidates
	canc = unlist(strsplit(comboo, "_"))[2]
	#rand_cind = as.data.table(filter(lncs_random_cindices, Cancer == canc))

	clin_combo_cind =  as.data.table(filter(lncsclinccombo_random_cindices, combo == comboo))

	#only genes
	#z = which(str_detect(rand_cind$lncRNA, "ENSG"))
	#rand_cind = rand_cind[z,]
	#colnames(rand_cind)[3] = "cindex"
	#z = which(is.na(rand_cind$cindex))
	#if(!(length(z)==0)){
	#rand_cind = rand_cind[-z,]}
	#z = which(rand_cind$lncRNA %in% allCands$gene)
	#if(!(length(z)==0)){
#		rand_cind = rand_cind[-z,]
#	}
#	rand_cind$type = "lncRNA_random"
#	rand_cind$combo = paste(rand_cind$lncRNA, rand_cind$Cancer, sep="_")
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
	clinical_cinds$combo = paste(clinical_cinds$lncRNA, clinical_cinds$Cancer, sep="_")

	lnc_cind$cindex = as.numeric(lnc_cind$cindex)
	clinical_cinds$cindex = as.numeric(clinical_cinds$cindex)
	#rand_cind$cindex = as.numeric(rand_cind$cindex)
	colnames(clin_combo_cind)[3] = "cindex"
	clin_combo_cind$cindex = as.numeric(clin_combo_cind$cindex)

	lnc_med = median(lnc_cind$cindex)
	clin_med = median(clinical_cinds$cindex)

	z = which(is.na(clin_combo_cind$cindex))
	if(!(length(z) == 0)){
		clin_combo_cind = clin_combo_cind[-z,]
	}
	clin_combo_cind_med = median(clin_combo_cind$cindex)

	#rbind lnc candidate and random lncs
	all_cindices = rbind(lnc_cind, clinical_cinds, clin_combo_cind)
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

	gene = get_name(unlist(strsplit(comboo, "_"))[1])

	all_cindices$type = factor(all_cindices$type, levels=c("clinical_variables", "lncRNA_canc", "lncRNA&clin"))
	all_cindices$gene= gene
	all_cindices$canc = canc

	all_cindices$wp_lnc_clinical = wp_lnc_clinical
	all_cindices$diff_meds_lnc_clinical = diff_meds_lnc_clinical
	all_cindices$wp_lnc_clinical_combo = wp_lnc_clinical_combo
	all_cindices$diff_meds_lnc_clinical_combo = diff_meds_lnc_clinical_combo
	all_cindices$lnc_med = lnc_med
	all_cindices$clin_med = clin_med
	all_cindices$combo_med = clin_combo_cind_med


	# Visualize: Specify the comparisons you want
	#my_comparisons <- list( c("clinical_variables", "lncRNA_canc"), c("lncRNA_canc", "lncRNA_random"),
	# c("clinical_variables", "lncRNA_random"),c("clinical_variables", "lncRNA&clin"))

	#boxplots
	#p = ggplot(data = all_cindices, aes(x = type, y = cindex)) +
    #geom_jitter(alpha = 0.3, color = "tomato") +
    #geom_boxplot() + labs(title=paste(gene, canc, "\nc-index vs random lncRNAs and clinical variables"))
	##  Add p-value
	#p = p + stat_compare_means(comparisons = my_comparisons) + geom_hline(yintercept=0.5, linetype="dashed", color = "red")
	#res = unlist(c(gene, canc, runs, wp_lnc_clinical, diff_meds_lnc_clinical,
	#wp_lnc_clinical_combo, diff_meds_lnc_clinical_combo,
#			lnc_med, clin_med, clin_combo_cind_med))
	#print(p)
#	return(res)
	return(all_cindices)
}


#pdf("all_lncRNA_cands_vs_random_lncRNAs.pdf")
random_lncs_vs_cand = llply(combos, get_comparison, .progress="text")
#dev.off()
random_lncs_vs_cand1 = ldply(random_lncs_vs_cand)
random_lncs_vs_cand1 = as.data.table(random_lncs_vs_cand1)
#random_lncs_vs_cand1$V4 = NULL
colnames(random_lncs_vs_cand1) = c("lncRNA", "cancer", "cindex", "type", "round", "combo", "gene", "canc", "wp_lnc_clinical",
	"diff_meds_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical_combo",
	"median_lncRNA", "median_clinical", "clin_combo_cind_med")
#random_lncs_vs_cand1$wp_lnc_random = as.numeric(random_lncs_vs_cand1$wp_lnc_random)
random_lncs_vs_cand1$wp_lnc_clinical = as.numeric(random_lncs_vs_cand1$wp_lnc_clinical)
random_lncs_vs_cand1$wp_lnc_clinical_combo = as.numeric(random_lncs_vs_cand1$wp_lnc_clinical_combo)

random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical)]
random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical, diff_meds_lnc_clinical, wp_lnc_clinical)]

#try new Figure 2E
#x-axis median difference between lncRNA and clinical variables
#y-axis p-value
full=random_lncs_vs_cand1
full$wp_lnc_clinical_fdr = p.adjust(full$wp_lnc_clinical, method="fdr")
full$wp_lnc_clinical_fdr = -log10(full$wp_lnc_clinical_fdr)
full$wp_lnc_clinical_combo = p.adjust(full$wp_lnc_clinical_combo, method="fdr")
z = which(full$wp_lnc_clinical_fdr == "Inf")
full$wp_lnc_clinical_fdr[z] = -log10(0.0000001)

write.table(full, file="summary_cross_validation_clinical_random_lncRNAs_Nov1.txt", quote=F, row.names=F, sep=";")

random_lncs_vs_cand1 = unique(random_lncs_vs_cand1[,c("cancer",
"gene", "canc", "wp_lnc_clinical",
	"diff_meds_lnc_clinical", "wp_lnc_clinical_combo", "diff_meds_lnc_clinical_combo",
	"median_lncRNA", "median_clinical", "clin_combo_cind_med")])

random_lncs_vs_cand1$wp_lnc_clinical_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical, method="fdr")
random_lncs_vs_cand1$wp_lnc_clinical_fdr = -log10(random_lncs_vs_cand1$wp_lnc_clinical_fdr)
random_lncs_vs_cand1$wp_lnc_clinical_combo = p.adjust(random_lncs_vs_cand1$wp_lnc_clinical_combo, method="fdr")
z = which(random_lncs_vs_cand1$wp_lnc_clinical_fdr == "Inf")
random_lncs_vs_cand1$wp_lnc_clinical_fdr[z] = -log10(0.0000001)

#z = which(random_lncs_vs_cand1$wp_lnc_clinical_fdr == "Inf")
#random_lncs_vs_cand1 = random_lncs_vs_cand1[-z,]

#random_lncs_vs_cand1$wp_lnc_random_fdr = p.adjust(random_lncs_vs_cand1$wp_lnc_random, method="fdr")
#random_lncs_vs_cand1$wp_lnc_random_fdr = -log10(random_lncs_vs_cand1$wp_lnc_random_fdr)

#plot 1 - lncs vs clinical
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv) = c("canc_type","cancer")
random_lncs_vs_cand1 = merge(random_lncs_vs_cand1, canc_conv, by="cancer")
random_lncs_vs_cand1$diff_meds_lnc_clinical = as.numeric(random_lncs_vs_cand1$diff_meds_lnc_clinical)
#random_lncs_vs_cand1$diff_meds_lnc_random = as.numeric(random_lncs_vs_cand1$diff_meds_lnc_random)
full = merge(full, canc_conv, by="cancer")
full$diff_meds_lnc_clinical = as.numeric(full$diff_meds_lnc_clinical)

random_lncs_vs_cand1 = random_lncs_vs_cand1[order(wp_lnc_clinical, diff_meds_lnc_clinical, wp_lnc_clinical)]
random_lncs_vs_cand1$type = factor(random_lncs_vs_cand1$type, levels=unique(random_lncs_vs_cand1$type))
random_lncs_vs_cand1$canc_type = factor(random_lncs_vs_cand1$canc_type, levels=unique(random_lncs_vs_cand1$canc_type))

full = full[order(wp_lnc_clinical, diff_meds_lnc_clinical, wp_lnc_clinical)]
full$type = factor(full$type, levels=unique(full$type))
full$canc_type = factor(full$canc_type, levels=unique(full$canc_type))

#add gene name
allCands$name = unlist(llply(allCands$gene, get_name))

z = which(random_lncs_vs_cand1$gene == "HOXA-AS4")
random_lncs_vs_cand1$gene[z] = "HOXA10-AS"
z = which(full$gene == "HOXA-AS4")
full$gene[z] = "HOXA10-AS"

#known lncRNAs
known_lncs = c("HOXB-AS2", "HOXA10-AS", "CCAT1", "BZRAP1-AS1", "AC114803.3", "AC097468.7")
colnames(random_lncs_vs_cand1)[2] = "name"
colnames(full)[7] = "name"

allCands$name[allCands$name == "HOXA-AS4"] = "HOXA10-AS"

random_lncs_vs_cand1 = merge(random_lncs_vs_cand1, allCands, by = c("cancer", "name"))
random_lncs_vs_cand1$HR = as.numeric(random_lncs_vs_cand1$HR)
random_lncs_vs_cand1$risk = ""
random_lncs_vs_cand1$risk[random_lncs_vs_cand1$HR > 1] = "Unfavourable"
random_lncs_vs_cand1$risk[random_lncs_vs_cand1$HR < 1] = "Favourable"

full = merge(full, allCands, by = c("cancer", "name"))
full$HR = as.numeric(full$HR)
full$risk = ""
full$risk[full$HR > 1] = "Unfavourable"
full$risk[full$HR < 1] = "Favourable"

write.csv(full, file="summary_cross_validation_clinical_random_lncRNAs_Nov1.csv", quote=F, row.names=F)

random_lncs_vs_cand1$HR = log2(random_lncs_vs_cand1$HR)
random_lncs_vs_cand1$sig = ""
random_lncs_vs_cand1$sig[random_lncs_vs_cand1$wp_lnc_clinical_fdr >= -log10(0.05)] = "Yes"
random_lncs_vs_cand1$sig[random_lncs_vs_cand1$wp_lnc_clinical_fdr < -log10(0.05)] = "No"
random_lncs_vs_cand1$sig = factor(random_lncs_vs_cand1$sig, levels=c("Yes", "No"))
random_lncs_vs_cand1$type = as.factor(random_lncs_vs_cand1$type)
random_lncs_vs_cand1$canc_type = as.factor(random_lncs_vs_cand1$canc_type)

#supp table 4
#write.csv(random_lncs_vs_cand1, file="/u/kisaev/Dec2020/SuppTable4_final_cands_wcindices.csv", quote=F, row.names=F)

#part a
mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

random_lncs_vs_cand1 = unique(random_lncs_vs_cand1[,c("HR", "diff_meds_lnc_clinical",
"name", "canc_type", "sig")])

pdf("/u/kisaev/Dec2020/figure2_e_lncRNA_cands_vs_clinical_variables.pdf", width=6, height=6)
g1 = ggplot(random_lncs_vs_cand1, aes(x=HR, y=diff_meds_lnc_clinical, label=name)) +
 geom_point(aes(fill=canc_type, colour=sig),
       pch=21, size=2)+theme_bw()+
geom_text_repel(data = subset(random_lncs_vs_cand1, name %in% known_lncs),
      min.segment.length = unit(0, 'lines'),
                     nudge_y = .2)+ labs(fill = "Cancer") +
geom_hline(yintercept=0, linetype="dashed", color = "black")+
geom_vline(xintercept=0, linetype="dashed", color = "black")+
labs(x="log2(HR)", y="med lncRNA c-index - med clinical c-index")+
 theme(legend.position="bottom", text = element_text(size=11)) +
 scale_fill_manual(values = mypal5[1:22]) +
 scale_colour_manual(values = c("black", "white"))+theme_bw()
 g1

dev.off()

#part b
#pdf("figure2_e_lncRNA_cands_vs_random_lncs_variables.pdf", width=6, height=6)
#ggplot(random_lncs_vs_cand1, aes(x=diff_meds_lnc_random, y=wp_lnc_random_fdr)) + geom_point(aes(colour=type))+
#geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
#geom_vline(xintercept=0, linetype="dashed", color = "black")+
#theme_bw() +
#labs(x="(median lncRNA candidate c-index) - (median random lncRNAs c-index)", y="-log10(adjusted p-val)")+
#scale_colour_manual(values = mypal5[1:22]) + theme(legend.position="bottom", text = element_text(size=12))
#dev.off()


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
