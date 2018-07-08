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

#ALK and the other gene of interest for fusion
z = which(colnames(pcg_rna) %in% c("ENSG00000171094", "ENSG00000178568", "patient", "canc"))
pcg_rna = pcg_rna[,z]

colnames(patient_table) = c("patient", "cancer")
pcg_rna = merge(pcg_rna, patient_table, by = "patient")

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg_rna$cancer))


#------ALK---------------------------------------------------------------------------------

#ALK patient --> "DO40069"
canc_pat = pcg_rna$cancer[which(pcg_rna$patient =="DO40069")]

#first compare to all cancers 
pat_alk_exp = pcg_rna$ENSG00000171094[which(pcg_rna$patient =="DO40069")]
everyone_else = (pcg_rna$ENSG00000171094[which(!(pcg_rna$patient =="DO40069"))])
everyone_else_mean = mean(pcg_rna$ENSG00000171094[which(!(pcg_rna$patient =="DO40069"))])

#then compare to just that cancer
alk_canc = subset(pcg_rna, cancer == canc_pat)
alk_canc_mean = mean(alk_canc$ENSG00000171094[which(!(alk_canc$patient =="DO40069"))])

library(broom)

#do wilcoxon test
#all cancers
w = wilcox.test(pcg_rna$ENSG00000171094[which(pcg_rna$patient =="DO40069")], 
	pcg_rna$ENSG00000171094[which(!(pcg_rna$patient =="DO40069"))])

tidy(w)
#405.5 vs 0.26

#just thryoid
w = wilcox.test(alk_canc$ENSG00000171094[which(alk_canc$patient =="DO40069")], 
	alk_canc$ENSG00000171094[which(!(alk_canc$patient =="DO40069"))])

tidy(w)

#------ERBB4--------------------------------------------------------------------------------

#ERBB4 patient --> "DO4449"
canc_pat = pcg_rna$cancer[which(pcg_rna$patient =="DO4449")]

#first compare to all cancers 
pat_alk_exp = pcg_rna$ENSG00000178568[which(pcg_rna$patient =="DO4449")]
everyone_else = (pcg_rna$ENSG00000178568[which(!(pcg_rna$patient =="DO4449"))])
everyone_else_mean = mean(pcg_rna$ENSG00000178568[which(!(pcg_rna$patient =="DO4449"))])

#then compare to just that cancer
alk_canc = subset(pcg_rna, cancer == canc_pat)
alk_canc_mean = mean(alk_canc$ENSG00000178568[which(!(alk_canc$patient =="DO4449"))])

library(broom)

#do wilcoxon test
#all cancers
w = wilcox.test(pcg_rna$ENSG00000178568[which(pcg_rna$patient =="DO4449")], 
	pcg_rna$ENSG00000178568[which(!(pcg_rna$patient =="DO4449"))])

tidy(w)
#1.086 vs 0.94 not significant 

#just thryoid
w = wilcox.test(alk_canc$ENSG00000178568[which(alk_canc$patient =="DO4449")], 
	alk_canc$ENSG00000178568[which(!(alk_canc$patient =="DO4449"))])

tidy(w)

#1.086 vs 2.2





















