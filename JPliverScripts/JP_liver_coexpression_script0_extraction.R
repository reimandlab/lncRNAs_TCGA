###JP_liver_coexpression_script0_extraction.R

#Author: Karina_Isaev
#Date_started: July 11th 2017
#dir: Thesis/

#Description:
#Processing RNA-Seq file from GTEx
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data

###Preamble###############################################
options(stringsAsFactors=F)

###Libraries##############################################
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

###Data#####################################################

#RNA-Seq file 
rna <- readRDS("gtex_expression_file.rds")

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#Sample ID Conversion file
conv <- ""

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)

###Processing#################################################

#1. Want to only look at ENSG genes in rna file
#split feature column and extract third component, write function and apply

#seperate first by "::"
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1], 1, ensg)
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

rna[,1] <- apply(rna[,1], 1, extract3) ; 
lncs[,1] <- apply(lncs[,1], 1, extract3) ; 

#UCSC file 
lncs <- as.data.frame(lncs)
z <- which(duplicated(ucsc$hg19.ensGene.name2))
ucsc <- ucsc[-z,] #60,234
#keep antisense, lincRNAs from this file 
z <- which(ucsc$hg19.ensemblSource.sourc %in% c("lincRNA", "antisense", "sense_intronic", "sense_overlapping"))
ucsc_lncs <- ucsc[z,]
dim(ucsc_lncs) #13,506
z <- which(ucsc$hg19.ensemblSource.sourc %in% "protein_coding")
ucsc_prots <- ucsc[z,]
dim(ucsc_prots) #19,002

#3. Before subsetting expression file to lncs and PCGs, make sure all samples 
#in this martix also have clinical data 

#using conv file 

rna <- as.data.frame(rna)

conv <- as.data.frame(conv)
z <- which(colnames(rna) %in% conv[,83])
rna <- rna[,c(1,z)] #1267 patients with RNA and Clinical Data 

#convert IDs to ICGC_donor_ids

jp_liver_samples <- fread("300_JP_ID.txt", data.table=F)

z <- which(conv[,83] %in% colnames(rna))
convert_ids <- function(col){
	z <- which(conv[,83] %in% col[[1]])
	id <- conv[z,4]
	return(id)
}
l <- as.data.frame(colnames(rna)[2:1268])
colnames(rna)[2:1268] <- apply(l, 1, convert_ids)

#now check that samples have clinical data - they all do

z <- which(colnames(rna) %in% jp_liver_samples[,1])
jp_rna <- rna[,c(1,z)]

#now check that samples have clinical data - they all do
#4. remove duplicates 
rna2 <- jp_rna[! jp_rna$feature %in% unique(jp_rna[duplicated(jp_rna$feature), "feature"]), ]
lncs2 <- lncs[! lncs$CAT_geneName %in% unique(lncs[duplicated(lncs$CAT_geneName), "CAT_geneName"]), ]

#5. Now seperate rna file into lncs and non-lncs 
z <- which(rna2[,1] %in% ucsc_lncs$hg19.ensGene.name2)
lncs_expression <- rna2[z,] #13479
z <- which(rna2[,1] %in% ucsc_prots$hg19.ensGene.name2)
pcg_expression <- rna2[z,] #18039 

#6. subset clin file to patients in expression matrix
#z <- which(clin$icgc_donor_id %in% colnames(rna2))
#clin <- clin[z,]
#using JP ids 

saveRDS(lncs_expression, "liver_jp_lncRNA_expression_6028.rds")
saveRDS(pcg_expression, "liver_jp_pcg_expression.rds")
#saveRDS(clin, "liver_jp_clinical.rds")

#converting to find JP liver ids 

submitter_donor_id