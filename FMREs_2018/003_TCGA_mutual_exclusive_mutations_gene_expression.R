library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(Hmisc)

#---TCGA directory 


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


###TCGA DATA---------------------------------------------------------------------------------------------------------

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)


#Data#-------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]

#split by cancer type 

z = which(str_detect(colnames(pcg), "ENSG"))	
pcg = as.data.frame(pcg)
pcg[,z] <- log1p(pcg[,z])

#2. Get lncRNA - median within each tissue type
tissues <- unique(pcg$type)

#RNA data -- pcawg
#pcg_rna = readRDS(file="all_rna_may8th.rds")

#subset to TP53 and ZKSCAN3
z = which(colnames(pcg) %in% c("ENSG00000141510", "ENSG00000189298", "patient", "type"))
pcg = pcg[,z]

#colnames(patient_table) = c("patient", "cancer")
#pcg_rna = merge(pcg_rna, patient_table, by = "patient")

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg$type))

get_cor = function(canc){
	print(canc)
	canc_data = subset(pcg, type == canc)
	colnames(canc_data)[c(2:3)] = c("ZKSCAN3", "TP53")
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
    #		conf.int = TRUE # Add confidence interval
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

#results = Filter(Negate(is.null), results)
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

table(results$sig_cor, results$cor_type)

results = as.data.table(results)
results = results[order(fdr)]

cancers = (unique(results$cancer))

#order by FDR then re-plot 
get_cor_plot = function(canc){
	
	fdr = results$fdr[which(results$cancer == canc)]
	print(canc)
	canc_data = subset(pcg, type == canc)
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
	sp + stat_cor(method = "spearman") + theme_bw() + ggtitle(paste(canc, "TP53 ~ ZKSCAN3", "fdr=", fdr))
		#FPKM is log1p transformed 
}
}


pdf("TCGA_cancers_tp53_zkscan3_correlations_june13.pdf")
llply(cancers, get_cor_plot, .progress="text")
dev.off()





































