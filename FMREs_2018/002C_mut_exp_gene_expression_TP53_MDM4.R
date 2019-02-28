library(plyr)
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

#[1] FMREs, 41 in total
#load("dfr_encode_results_signf.rsav")
#head(dfr_encode_results_signf)
#fmres = dfr_encode_results_signf

# 84 unique patients have the FMRE "ENCODEmerge::chr17:57914548-57924685::NA::NA" 

#[2] coding drivers - need to filter for cancer_type==“PANCANCER” and 
#element_type==“gc19_pc.cds”, 47 in total

#old file:
#load("all_results_signf.rsav")
#head(all_results_signf)
#all_results_signf = as.data.table(all_results_signf)
#all_results_signf = filter(all_results_signf, cancer_type == "PANCANCER")
#all_results_signf = filter(all_results_signf, element_type == "gc19_pc.cds")
#coding_drivers = all_results_signf

#cds mutations new june 12 
coding_drivers = fread("july2_cds_drivers.txt")

#[3] mutations in all CRMs, subset of these are FMREs

#old file:

#load("encode_merge_oct2016_mutations__PANCANCER_in_elements.rsav")
#head(mutations_in_elements)
#mutations_in_crms = mutations_in_elements

#new file June 12: 

load("july12_encode_merge__patient_element_snv_list.rsav")
mutations_in_crms = (patient_element_snv_list)

#crm mutations new june 12 
#fmre mutations 
fmres = fread("july12_fmre_drivers.txt")

#[4] mutations in all CDS, subset of these are CDS drivers
load("july12_gc19_pc.cds__patient_element_snv_list.rsav")
cds_mutations = patient_element_snv_list

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[] cds_lnc_drivers

#[] july12_gc19_pc.cds__coords.rsav
load("july12_gc19_pc.cds__coords.rsav")
coords = element_coords
coords = as.data.table(coords)

#[] july12_lncrna.ncrna__coords.rsav
load("july12_lncrna.ncrna__coords.rsav")
lnc_coords = element_coords
lnc_coords = as.data.table(lnc_coords)

#[] july12_encode_merge__coords.rsav
load("july12_encode_merge__coords.rsav")
element_coords = element_coords
element_coords = as.data.table(element_coords)

#[] july12_lncrna.ncrna__patient_element_snv_list.rsav
load("july12_lncrna.ncrna__patient_element_snv_list.rsav")
lncrna_mutations = patient_element_snv_list

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

#tp53 muts 
tp53_muts = readRDS("tp53_pcawg_muts.rds")

#patients with zkscan FMRE --> 6:27870028
z = which(names(mutations_in_crms) == "ENCODEmerge::chr6:27870028-27871319::NA::NA")
fmre_pats = mutations_in_crms[[z]]

#which TP53 mutations do these patients have?
z = which(tp53_muts$Donor_ID %in% fmre_pats)
tp53_fmre_muts = tp53_muts[z,]
pats = unique(tp53_fmre_muts$Donor_ID)
pats[13]
#[1] "DO28763" <--- supposingly has TP53 mutation but not in CDS file 
tp53_fmre_muts[which(tp53_fmre_muts$Donor_ID == pats[13])]

#add 5 lncRNAs to mutations in cds
lncs = c("lncrna.ncrna::gencode::RN7SK::ENSG00000202198.1",
"lncrna.ncrna::gencode::NEAT1::ENSG00000245532.4",
"lncrna.ncrna::gencode::MALAT1::ENSG00000251562.3",
"lncrna.ncrna::gencode::RPPH1::ENSG00000259001.2",
"lncrna.ncrna::gencode::Z95704.4::ENSG00000248302.2")

z = which(names(lncrna_mutations) %in% lncs)
lncrna_mutations = lncrna_mutations[z]
#bind it to the cds list
mutations_in_cds = c(mutations_in_cds, lncrna_mutations)

#Compare ratios------------------------------------------------

unique_cds = unique(names(mutations_in_cds)) #48 + 5 lncRNAs 
unique_fmre = unique(names(mutations_in_crms)) #30

mdm4_fmre = "ENCODEmerge::chr1:204475015-204476599::NA::NA"

#Chcekc do we expect overlaps?
all_lnc_cord = as.data.table(filter(lnc_coords, id %in% unique_cds))
colnames(all_lnc_cord)[1:3] = c("Chr", "Start", "End")
all_lnc_cord = makeGRangesFromDataFrame(all_lnc_cord)

all_cds_cord = as.data.table(filter(coords, id %in% unique_cds))
colnames(all_cds_cord)[1:3] = c("Chr", "Start", "End")
all_cds_cord = makeGRangesFromDataFrame(all_cds_cord)

all_fmre_cord = as.data.table(filter(element_coords, id %in% unique_fmre))
colnames(all_fmre_cord)[1:3] = c("Chr", "Start", "End")
all_fmre_cord = makeGRangesFromDataFrame(all_fmre_cord)

findOverlaps(all_fmre_cord, all_cds_cord)
findOverlaps(all_fmre_cord, all_lnc_cord) #so only overlaps between FMREs and lncRNAs 

#Analysis part2---------------------------------------------------

#- can you tell me if TP53 and ZKSCAN3 are correlated transcriptionally? cancer type by cancer type; PCAWG initially. 

#RNA data 
pcg_rna = readRDS(file="all_rna_may8th.rds")

#subset to TP53 and MDM4
z = which(colnames(pcg_rna) %in% c("ENSG00000141510", "ENSG00000198625", "patient", "canc"))
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
z = which(names(mutations_in_crms) == "chr1:204475015")

#------patients with MDM4 FMRE mutation ----------------------------------------------------

mut_fmre = mutations_in_crms[[z]] 

#-------------------------------------------------------------------------------------------

#keep just gene name for CDS
clean_gene = function(gene){
	#gene = results_pairs$CDS_mut[1]
	r =  unlist(strsplit(gene, "::"))[3]
	return(r)
}

names(cds_mutations) = unlist(llply(names(cds_mutations), clean_gene))
z = which(names(cds_mutations) == "TP53")

#------patients with TP53 CDS mutation -----------------------------------------------------

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

	z = which(colnames(canc_exp) == "ENSG00000141510")
	colnames(canc_exp)[z] = "TP53"
	z = which(colnames(canc_exp) == "ENSG00000198625")
	colnames(canc_exp)[z] = "MDM4"

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

	canc_exp$mut_code = factor(canc_exp$mut_code, levels =c("None", "FMRE", "CDS", "Both"))

	#remove any groups made up of only 1 patient 
	t = as.data.table(table(canc_exp$mut_code))
	t = as.data.table(filter(t, N > 1))
	canc_exp = as.data.table(filter(canc_exp, mut_code %in% t$V1))

	#plot 1, y = TERT exp
	
	#-----------------------
	#A, x = IDH1 mut yes/no
	#-----------------------
	a1 <- ggboxplot(canc_exp, x = "mut_code", y = "MDM4",
         color = "mut_code", 
         palette = "jco", title = paste(canc, "MDM4 exp ~ mut"),  
          add = "jitter") + stat_n_text() + ylab("log2 expression (FPKM-UQ)")
	# Change method
	a1 = a1 + stat_compare_means(method = "anova") + theme_bw()
	print(a1)
	
	#plot 2, y = IDH1 exp

	#-----------------------
	#A, x = TERT mut yes/no
	#-----------------------

	a1 <- ggboxplot(canc_exp, x = "mut_code", y = "TP53",
         color = "mut_code", 
         palette = "jco", title = paste(canc, "TP53 exp ~ mut"), 
          add = "jitter") + stat_n_text() + ylab("log2 expression")
	# Change method
	a2 = a1 + stat_compare_means(method = "anova") + theme_bw()
	print(a2)
}
}
}

pdf("CDS_FMRE_comutation_summary_all_cancers_MDM4_TO53_feb28KI.pdf")
llply(cancers, check_zkscan_exp, .progress="text")
dev.off()





























