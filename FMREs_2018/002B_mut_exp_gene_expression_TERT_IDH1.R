library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(Hmisc)
library(patchwork)
library(EnvStats)

#Data---------------------------------------------------

#cds mutations new july 2
coding_drivers = fread("july2_cds_drivers.txt")

#[3] mutations in all CRMs, subset of these are FMREs

#old file:

#load("encode_merge_oct2016_mutations__PANCANCER_in_elements.rsav")
#head(mutations_in_elements)
#mutations_in_crms = mutations_in_elements

#new file June 12: 

load("july2_encode_merge__patient_element_snv_list.rsav")
mutations_in_crms = (patient_element_snv_list)

#crm mutations new june 12 
#fmre mutations 
fmres = fread("july2_fmre_drivers.txt")

#[4] mutations in all CDS, subset of these are CDS drivers
load("july2_gc19_pc.cds__patient_element_snv_list.rsav")
cds_mutations = patient_element_snv_list

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(names(mutations_in_crms) %in% fmres$id)
mutations_in_crms = mutations_in_crms[z]

#2. Which patients have coding mutations 
z = which(names(cds_mutations) %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z]

#RESULTS_from_001_------------------------------------------------

results_pairs = fread("819_fmre_cds_pairs_fishers_analysis_with_lncRNAs_July12nd_KI.txt")

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

#subset to IDH1 and TERT
z = which(colnames(pcg_rna) %in% c("ENSG00000138413", "ENSG00000164362", "patient", "canc"))
pcg_rna = pcg_rna[,z]

colnames(patient_table) = c("patient", "cancer")
pcg_rna = merge(pcg_rna, patient_table, by = "patient")

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg_rna$cancer))

#- chr6:278 is associated with ZKSCAN3. there should be three ovarian samples with PCAWG RNAdata and mutation in chr6:278. Is TP53 lower in these samples? 

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
z = which(names(mutations_in_crms) == "chr5:1294859")

#------patients with ZKSCAN3 FMRE mutation -------------------------------------------------
mut_fmre = mutations_in_crms[[z]] 
#-------------------------------------------------------------------------------------------

#keep just gene name for CDS
clean_gene = function(gene){
	#gene = results_pairs$CDS_mut[1]
	r =  unlist(strsplit(gene, "::"))[3]
	return(r)
}

names(cds_mutations) = unlist(llply(names(cds_mutations), clean_gene))
z = which(names(cds_mutations) == "IDH1")
#------patients with ZKSCAN3 FMRE mutation -------------------------------------------------
mut_cds = cds_mutations[[z]]
#-------------------------------------------------------------------------------------------

check_zkscan_exp = function(canc){

	canc_exp = subset(pcg_rna, cancer == canc)
	z = which(canc_exp$patient %in% mut_fmre)
	if(!(length(z)==0)){
	canc_exp$fmre = ""
	canc_exp$fmre[z] = "FMRE"
	canc_exp$fmre[-z] = "noFMRE"
	
	z = which(canc_exp$patient %in% mut_cds)
	if(!(length(z)==0)){
	canc_exp$cds = ""
	canc_exp$cds[z] = "CDS_mut"
	canc_exp$cds[-z] = "noCDS_mut"

	z = which(colnames(canc_exp) == "ENSG00000138413")
	colnames(canc_exp)[z] = "IDH1"
	z = which(colnames(canc_exp) == "ENSG00000164362")
	colnames(canc_exp)[z] = "TERT"

	canc_exp$both = ""
	z = which((canc_exp$cds == "CDS_mut") & (canc_exp$fmre == "FMRE")) 
	if(!(length(z)==0)){
		canc_exp$both[z] = "both"
		canc_exp$both[-z] = "other"
	}

	canc_exp$none = ""
	z = which((canc_exp$cds == "noCDS_mut") & (canc_exp$fmre == "noFMRE")) 
	if(!(length(z)==0)){
		canc_exp$none[z] = "neither"
		canc_exp$none[-z] = "other"
	}
	
	canc_exp$mut_code = ""
	for(i in 1:nrow(canc_exp)){
		cds = canc_exp$cds[i]
		f = canc_exp$fmre[i]
		both = (cds == "CDS_mut") & (f=="FMRE")
		if(both){
			stat = "Both"
		}
		f_only = (cds == "noCDS_mut") & (f=="FMRE")
		if(f_only){
			stat = "FMRE"
		}
		cds_only = (cds == "CDS_mut") & (f=="noFMRE")
		if(cds_only){
			stat = "CDS"
		}
		none = (cds == "noCDS_mut") & (f=="noFMRE")
		if(none){
			stat = "None"
		}
		canc_exp$mut_code[i] = stat
	}


	canc_exp[,2:3] = log2(canc_exp[,2:3])
	ord = aggregate(canc_exp[, 3], list(canc_exp$mut_code), median)
	ord = as.data.table(ord)
	ord = ord[order(x)]

	#plot 1, y = TERT exp
	
	#-----------------------
	#A, x = IDH1 mut yes/no
	#-----------------------
	a1 <- ggboxplot(canc_exp, x = "mut_code", y = "TERT",
         color = "mut_code", order = ord$Group.1, 
         palette = "jco", title = paste(canc, "TERT exp ~ mut"),  
          add = "jitter") + stat_n_text() + ylab("log2 expression")
	# Change method
	a1 = a1 + stat_compare_means(method = "anova") + theme_bw()
	print(a1)
	
	#plot 2, y = IDH1 exp

	#-----------------------
	#A, x = TERT mut yes/no
	#-----------------------
	ord = aggregate(canc_exp[, 2], list(canc_exp$mut_code), median)
	ord = as.data.table(ord)
	ord = ord[order(x)]
	a1 <- ggboxplot(canc_exp, x = "mut_code", y = "IDH1",
         color = "mut_code", order = ord$Group.1, 
         palette = "jco", title = paste(canc, "IDH1 exp ~ mut"), 
          add = "jitter") + stat_n_text() + ylab("log2 expression")
	# Change method
	a2 = a1 + stat_compare_means(method = "anova") + theme_bw()
	print(a2)
}
}
}

pdf("CDS_FMRE_comutation_summary_all_cancers_TERT_IDH1.pdf")
llply(cancers, check_zkscan_exp, .progress="text")
dev.off()





























