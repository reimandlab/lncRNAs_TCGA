library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(Hmisc)


#Data---------------------------------------------------

#cds mutations new june 12 
coding_drivers = fread("cds_drivers.txt")

#[] mutations in all CRMs, subset of these are FMREs

#new file June 12: 

load("encode_merge__patient_element_snv_list.rsav")
mutations_in_crms = (patient_element_snv_list)

#crm mutations new june 12 
#fmre mutations 
fmres = fread("fmre_drivers.txt")

#[4] mutations in all CDS, subset of these are CDS drivers
load("gc19_pc.cds__patient_element_snv_list.rsav")
cds_mutations = patient_element_snv_list

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion<- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(names(mutations_in_crms) %in% fmres$id)
mutations_in_crms = mutations_in_crms[z]

#2. Which patients have coding mutations 
z = which(names(cds_mutations) %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z]

#RESULTS_from_001_------------------------------------------------

results_pairs = fread("686_fmre_cds_pairs_fishers_analysis_June12th_KI.txt")

library(plyr)
library(dplyr)

#keep just gene name for CDS
clean_gene = function(gene){
	#gene = results_pairs$CDS_mut[1]
	r =  unlist(strsplit(gene, "::"))[3]
	return(r)
}

results_pairs$CDS_mut = unlist(llply(results_pairs$CDS_mut, clean_gene))

clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "::"))[2]
	return(r)
}

results_pairs$FMRE_mut = unlist(llply(results_pairs$FMRE_mut, clean_fmre))
#results_pairs$FMRE_mut = as.character(results_pairs$FMRE_mut)

#1. order by significance in driver analysis 
#using fdr_element column for this 

fmres$id = llply(fmres$id, clean_fmre)
fmres = as.data.table(fmres)
fmres = fmres[order(fdr_element)]

order= as.character(fmres$id)

#relevel order 
results_pairs$FMRE_mut <- factor(results_pairs$FMRE_mut, levels = order)

#relevel pcg order
coding_drivers$id = llply(coding_drivers$id, clean_gene)
coding_drivers = as.data.table(coding_drivers)
coding_drivers = coding_drivers[order(-fdr_element)]
order = as.character(coding_drivers$id)

#relevel order
results_pairs$CDS_mut <- factor(results_pairs$CDS_mut, levels = order)

#2. color is pvalue (fisher's pvalue) FDR -> -log10 FDR 
#FDR > 0.1 is NA 

results_pairs$fdr_plotting = -log10(results_pairs$fdr)
results_pairs$fdr_plotting[results_pairs$fdr >0.1] = NA

#now just need to plot heatmap 
results_pairs$CDS_mut = as.character(results_pairs$CDS_mut)
results_pairs = as.data.frame(results_pairs)

#shorten fmre name
clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "-"))[1]
	return(r)
}
results_pairs$FMRE_mut = unlist(llply(as.character(results_pairs$FMRE_mut), clean_fmre))

fmres$id = llply(fmres$id, clean_fmre)
fmres = as.data.table(fmres)
fmres = fmres[order(fdr_element)]

order= as.character(fmres$id)

#relevel order 
results_pairs$FMRE_mut <- factor(results_pairs$FMRE_mut, levels = order)

#Analysis part2---------------------------------------------------

#- can you tell me if TP53 and ZKSCAN3 are correlated transcriptionally? cancer type by cancer type; PCAWG initially. 

#RNA data 
pcg_rna = readRDS(file="all_rna_may8th.rds")

#subset to TP53 and ZKSCAN3
z = which(colnames(pcg_rna) %in% c("ENSG00000141510", "ENSG00000189298", "patient", "canc"))
pcg_rna = pcg_rna[,z]

colnames(patient_table) = c("patient", "cancer")
pcg_rna = merge(pcg_rna, patient_table, by = "patient")

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg_rna$cancer))

get_cor = function(canc){
	print(canc)
	canc_data = subset(pcg_rna, cancer == canc)
	colnames(canc_data)[c(2:3)] = c("TP53", "ZKSCAN3")
	canc_data[,2:3] = log1p(canc_data[,2:3])

	if(dim(canc_data)[1] >=4){

	#calculate correlation via linear model 
	lm_pval = summary(lm(canc_data$TP53 ~ canc_data$ZKSCAN3))$coefficients[2,4]
	lm_estimate = summary(lm(canc_data$TP53 ~ canc_data$ZKSCAN3))$coefficients[2,1]
	
	#calculate correlation via spearman cor 
	spearman_rho = rcorr(canc_data$TP53 , canc_data$ZKSCAN3, type=c("spearman"))$r[2]
	spearman_pval = rcorr(canc_data$TP53 , canc_data$ZKSCAN3, type=c("spearman"))$P[2]

	#plot 
	# Scatter plot with correlation coefficient
	#:::::::::::::::::::::::::::::::::::::::::::::::::
	#sp <- ggscatter(canc_data, x = "ZKSCAN3", y = "TP53",
   	#	add = "reg.line",  # Add regressin line
   	#	add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   	#	conf.int = TRUE # Add confidence interval
   	#)
		# Add correlation coefficient
	#print(sp + stat_cor(method = "spearman") + theme_bw() + ggtitle(paste(canc, "TP53 ~ ZKSCAN3")))
		#FPKM is log1p transformed 

	row = c(canc, lm_estimate, lm_pval, spearman_rho, spearman_pval)
	names(row) = c("cancer", "lm_estiamte", "lm_pval", "spearman_rho", "spearman_pval")
	return(row)
}
}

results = llply(cancers, get_cor, .progress="text")

results = Filter(Negate(is.null), results)
#combine into one dataframe
results = do.call(rbind.data.frame, results)

colnames(results) = c("cancer", "lm_estiamte", "lm_pval", "spearman_rho", "spearman_pval")
results$sig_cor = ""
results$sig_cor[as.numeric(as.character(results$spearman_pval)) <= 0.05] = "sig"
results$sig_cor[as.numeric(as.character(results$spearman_pval)) > 0.05] = "NOTsig"

results$cor_type = ""
results$cor_type[as.numeric(as.character(results$spearman_rho)) > 0] = "pos"
results$cor_type[as.numeric(as.character(results$spearman_rho)) < 0] = "neg"
results$fdr = p.adjust(as.numeric(as.character(results$spearman_pval)), method="fdr")
results = as.data.table(results)
results = results[order(fdr)]

cancers = (unique(results$cancer))

#order by FDR then re-plot 
get_cor_plot = function(canc){
	
	fdr = results$fdr[which(results$cancer == canc)]
	print(canc)
	canc_data = subset(pcg_rna, cancer == canc)
	colnames(canc_data)[c(2:3)] = c("TP53", "ZKSCAN3")
	canc_data[,2:3] = log1p(canc_data[,2:3])

	if(dim(canc_data)[1] >=4){
	#plot 
	# Scatter plot with correlation coefficient
	#:::::::::::::::::::::::::::::::::::::::::::::::::
	sp <- ggscatter(canc_data, x = "ZKSCAN3", y = "TP53",
   		add = "reg.line",  # Add regressin line
   	    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   		conf.int = TRUE # Add confidence interval
   	)
		# Add correlation coefficient
	sp + stat_cor(method = "spearman") + theme_bw() + ggtitle(paste(canc, "TP53 ~ ZKSCAN3", "fdr=", round(fdr, digits=4)))
		#FPKM is log1p transformed 
}
}


pdf("pcawg_cancers_tp53_zkscan3_correlations_june13.pdf")
llply(cancers, get_cor_plot, .progress="text")
dev.off()


table(results$sig_cor, results$cor_type)
        
#         neg pos
#  NOTsig   4  10
#  sig      1   5

#sig cancers = Bladder-TCC, CNS-Oligo , ColoRect-AdenoCA , Head-SCC, Lung-SCC, Ovary-AdenoCA 

#- chr6:278 is associated with ZKSCAN3. there should be three ovarian samples with PCAWG RNAdata and mutation in chr6:278. Is TP53 lower in these samples? 

ovary = subset(pcg_rna, cancer == "Ovary-AdenoCA")

#which patients have FMRE mut 
clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "::"))[2]
	return(r)
}
names(mutations_in_crms) = llply(names(mutations_in_crms), clean_fmre)

#shorten fmre name
clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "-"))[1]
	return(r)
}
names(mutations_in_crms) = unlist(llply(as.character(names(mutations_in_crms)), clean_fmre))

z = which(names(mutations_in_crms) == "chr6:27870028")

mut = mutations_in_crms[[z]]

z = which(ovary$patient %in% mut)

ovary$fmre = ""
ovary$fmre[z] = "FMRE"
ovary$fmre[-z] = "noFMRE"
ovary[,2] = log1p(ovary[,2])

pdf("TP53_exp_FMRE_vs_no_OVARY_PCAWG.pdf")
p <- ggboxplot(ovary, x = "fmre", y = "ENSG00000141510",
         color = "fmre", 
         palette = "jco", title = paste("Ovarian cancer - TP53 exp ~ FMRE", table(ovary$fmre)[1], "muts", table(ovary$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test") + theme_bw()

dev.off()


#- Of the ovarian Tp53 mutated samples, are there differences on zkscan3 expression relative to samples with WT tp53?

ovary = subset(pcg_rna, cancer == "Ovary-AdenoCA")

#which OV patients have TP53 mutation? 

#keep just gene name for CDS
clean_gene = function(gene){
	#gene = results_pairs$CDS_mut[1]
	r =  unlist(strsplit(gene, "::"))[3]
	return(r)
}

names(cds_mutations) = unlist(llply(names(cds_mutations), clean_gene))

z = which(names(cds_mutations) == "TP53")

mut = cds_mutations[[z]]

z = which(ovary$patient %in% mut)

ovary$tp53_mut = ""
ovary$tp53_mut[z] = "TP53_mut"
ovary$tp53_mut[-z] = "TP53_no_mut"
ovary[,2] = log1p(ovary[,2])
ovary[,3] = log1p(ovary[,3])

pdf("TP53_exp_FMRE_vs_no_OVARY_PCAWG.pdf")
p <- ggboxplot(ovary, x = "tp53_mut", y = "ENSG00000189298",
         color = "tp53_mut", 
         palette = "jco", title = paste("Ovarian cancer - ZKSCAN3 exp ~ TP53 mutation", table(ovary$tp53_mut)[1], "muts", table(ovary$tp53_mut)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test") + theme_bw()

dev.off()






























