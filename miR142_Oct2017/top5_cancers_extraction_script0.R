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
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

#---------------------------------------------------------
#Processing
#---------------------------------------------------------

# [1] Match RNA-Seq column names to conversion file in 
# order to get icgc donor IDs

#conversion$normal_rna_seq_submitter_sample_id = 63
#length(which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))) = 161
#length(which(conversion$tumor_rna_seq_submitter_sample_id %in% colnames(rna))) = 819
#length(which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))) = 1267

submitter <- conversion$tumor_rna_seq_submitter_sample_id[which(conversion$tumor_rna_seq_submitter_sample_id %in% colnames(rna))]
aliqot <- conversion$tumor_rna_seq_aliquot_id[which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))]

length(which(aliqot %in% submitter))
#[1] 819

#so aliqot has all the same IDs from submitter and so can just use aliqot ids 

###"TUMOUR SAMPLES"
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z]

lymph = clin$icgc_donor_id[clin$histology_tier2=="Lymphoid"]
z <- which(tum_pats %in% lymph)
tum_pats = tum_pats[z]

mts = c("DO27833" , "DO52682" , "DO52685" , "DO27851" , "DO52672" , "DO52692", 
 "DO27813" , "DO221124", "DO52717" , "DO52675" , "DO52647" , "DO27855", 
"DO27847" , "DO222308" ,"DO222305" ,"DO27835" , "DO52679")

##Convert RNA-Seq IDs to ICGC donor IDs
for(i in 1:ncol(rna)){
	z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna)[i])
	if(!(length(z)==0)){
		colnames(rna)[i] <- conversion$icgc_donor_id[z]
	}
}

#keep only lymphg patients 
z <- which(colnames(rna) %in% tum_pats)
rna = rna[,c(1,z)]

#Change gene names to match ensembl IDs 
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1:2], 1, ensg)
z <- which(check==TRUE)
rna <- rna[z,]

#2. Check how many IDs match the lncRNA IDs
#none match while trancript number is present 
#remove ie, ENSG00000201285.1 --> ENSG00000201285
#in both rna file and lncs file

extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract3) ; 

#Remove duplicate genes 
z <- which(duplicated(rna[,1]))
genes <- rna[z,1]
z <- which(rna[,1] %in% genes)
rna <- rna[-z,]

#change to HUGO gene names 
ucsc = ucsc[,c(6, 8)]
colnames(ucsc)[1] = "feature"

rna = merge(rna, ucsc, by = "feature")
colnames(rna)[182] = "Hugo"
saveRDS(rna, "180lymphomaPatients_RNA-Seqfile.rds")



