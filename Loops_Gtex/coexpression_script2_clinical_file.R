###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 13th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data


#SCRIPT - PLOTTING VIOLIN PLOTS FOR LNCRNAS IN LOOPS WITH MEDIAN EXPRESSION GREATER THAN 2RPKM 


#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

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

#These are gene expression files based on the original RNA-Seq
#PCAWG file
lnc <- readRDS("gtex_lncRNA_expression_12187.rds")
all <- readRDS("gtex_pcg_expression_17680.rds")

#reference genes
ref_genes_lncs <- lnc[,1:2]
ref_genes_all <- all[,1:2]

#change feature column with gene names so that they are the rownames
rownames(lnc) <- lnc[,1] ; lnc <- lnc[,-c(1:2)] #13479 lncRNAs as defined by ensembl "lincRNA", "antisense", "sense_intronic", "sense_overlapping"
rownames(all) <- all[,1] ; all <- all[,-c(1:2)] #18039 PCGs as defined by ensembl "protein_coding"

#loops 
loops <- fread("processed_loop_fileKI.txt", data.table=F)

#ucsc gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#---------------------------------------------------------
#Processing lncRNAs and PCGs 
#---------------------------------------------------------

loops$gene_type <- ""

for(i in 1:nrow(loops)){
	id <- loops$Gene[i]
	z <- which(ucsc[,6] %in% id)
	loops$gene_type[i] <- ucsc[z,7]
}

table(loops$Side, loops$gene_type)

#Loop lncRNAs genes 
z <- which(rownames(lnc) %in% loops$Gene)
#******************************************************************
loop_lnc <- lnc[z,] #2021 lncRNAs in loops with gene expression data 
#******************************************************************

#Loop PCGs genes 
z <- which(rownames(all) %in% loops$Gene)
#******************************************************************
loop_all <- all[z,] #5264 PCGs in loops with gene expression data 
#******************************************************************


#------------------------------------------------------------------
#Processing clinical file - inlcude samples that have expression
#------------------------------------------------------------------
z <- which(clin[,1] %in% colnames(lnc))
clin <- clin[z,]

tissues <- as.data.frame(clin[,c(1,6)])
c <- as.data.frame(table(tissues$SMTS))
c[,1] <- as.character(c[,1])
c[,2] <- as.numeric(c[,2])

#order
or <- order(c[,2])
c <- c[or, ]

pdf("GTEx_tissue_distribution_inExpressionFIle_July13th.pdf", width=10)
g <- ggbarplot(c, "Var1", "Freq", xlab="Tissue",fill="Var1")
ggpar(g, x.text.angle=60, yticks.by =100, legend = "none")
dev.off()


#PCAWG Cancers
cancers <- c("Biliary", "Bladder", "Bone","Breast" , "Cervix" , "CNS", "Colon", "Esophagus",
"Head"  , "Kidney" ,   "Liver" , "Lung" , "Lymph" , "Myeloid" , "Ovary" , "Pancreas",
"Prostate" ,    "Skin" , "Stomach" ,  "Thyroid"  , "Uterus", "Cervix Uteri", "Brain")	

z <- which(c[,1] %in% cancers)
cancers_keep <- c[z,1]

