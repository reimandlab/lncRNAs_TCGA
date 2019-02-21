#top5_cancers_extraction_script0.R

#Karina Isaev
#July 19th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script0 - extracts cancer types and plot distributions 

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library("colorout")
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

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

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
clin = fread("pcawg_specimen_histology_August2016_v9.tsv", data.table=F)

conversion <- fread("pcawgConversion.tsv", data.table=F)

more_clin = readRDS("Jan26_PCAWG_clinical")
more_clin = fread("donor.all_projects.tsv")

merge_cols = colnames(more_clin)[which(colnames(more_clin) %in% colnames(clin))]

clin = merge(clin, more_clin, by = merge_cols)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)


#---------------------------------------------------------
#Processing
#---------------------------------------------------------

#1. get cancer specific data 
clin = as.data.table(clin)
z = which(duplicated(clin$icgc_donor_id))
if(!(length(z)==0)){
clin = clin[-z,]}
cancers = unique(clin$histology_abbreviation)

#2. summarize clin variables in simple charts per cancer type 
get_charts = function(canc){
	print(canc)
	dat = as.data.table(filter(clin, histology_abbreviation == canc))
	#summarize each variables 
	cols = c("project_code", "histology_abbreviation", "histology_tier1", "histology_tier2", "histology_tier3", "histology_tier4",
		"tumour_histological_code", "tumour_histological_type", "tumour_stage", "tumour_grade", "percentage_cellularity", "level_of_cellularity", 
		"specimen_donor_treatment_type", "project_code", "donor_sex", "donor_vital_status", "donor_diagnosis_icd10", "first_therapy_type", 
		"first_therapy_response", "donor_age_at_diagnosis", "tobacco_smoking_history_indicator", "alcohol_history")
	
	get_summary = function(col){
		print(col)
		z = which(colnames(dat)==col)
		if(!(length(z))==0){

		if(!(is.numeric(unlist(dat[,..z])))){
			s = (table(dat[,..z]))
			s = as.data.table(s)
			colnames(s) = c(col, "freq")
			s = s[order(freq)]
			s$cancer = canc
		}
		
		if(is.numeric(unlist(dat[,..z]))){
			s = (summary(dat[,..z]))
			s = as.data.table(s)
			s$V1 = NULL
			s$V2 = NULL
			colnames(s)[1] = col
			s$N = NULL
			s$cancer = canc
		}
	
		return(s)
		}
	}

	cols_summaries = llply(cols, get_summary)
	cols_summaries = Filter(Negate(is.null), cols_summaries)

	library(gridExtra)
	library(grid)
	library(grid)  	

	for(i in 1:length(cols_summaries)){

	file = paste(canc, i, "new_sum.pdf", sep="_")	
	pdf(file)
	print(grid.table((cols_summaries[[i]])))
	dev.off()
	
	}
	print("done cancer")
}


llply(cancers, get_charts, .progress="text")



#just liver 
z = which(str_detect(conversion$dcc_project_code, "LIAD|LIHM|LICA|LIHC|LINC|LICA|LIRI"))
livers = conversion[z,]

pdf("liver_pcawg_countries.pdf", width=5, height=5)
t = as.data.table(table(livers$dcc_project_code))
colnames(t) = c("country-project","freq")
t = t[order(freq)]
grid.table(t)
dev.off()



z = which(str_detect(new_pcawg$project_code, "LIAD|LIHM|LICA|LIHC|LINC|LICA|LIRI"))
livers = new_pcawg[z,]
livers = as.data.table(filter(livers, study_donor_involved_in == "PCAWG"))










