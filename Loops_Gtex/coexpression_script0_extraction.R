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
rna <- as.data.frame(rna)

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

#seperate first by "."
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "\\..*"))
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 
lncs[,1] <- apply(lncs[,1], 1, extract) ; 

#2. Check how many IDs match the lncRNA IDs

#UCSC file 
lncs <- as.data.frame(lncs)
z <- which(duplicated(ucsc$hg19.ensGene.name2))
ucsc <- ucsc[-z,] #60,234
#keep antisense, lincRNAs from this file 
z <- which(ucsc$hg19.ensemblSource.sourc %in% c("lincRNA", "antisense"))
ucsc_lncs <- ucsc[z,]
dim(ucsc_lncs) #12,563
z <- which(ucsc$hg19.ensemblSource.sourc %in% "protein_coding")
ucsc_prots <- ucsc[z,]
dim(ucsc_prots) #19,002

#3. remove duplicates 
rna2 <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]
lncs2 <- lncs[! lncs$CAT_geneName %in% unique(lncs[duplicated(lncs$CAT_geneName), "CAT_geneName"]), ]

#4. Now seperate rna file into lncs and non-lncs 
z <- which(rna2[,1] %in% ucsc_lncs$hg19.ensGene.name2)
lncs_expression <- rna2[z,] #12187
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