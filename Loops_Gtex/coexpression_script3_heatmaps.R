###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 14th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data


#SCRIPT - PLOTTING HEATMAP OF GENE MEDIAN EXPRESSION BY TISSUE TYPE
#SEE IF CAN CAPTURE TISSUE SPECIFICTY AND ALSO IF WE ARE LOSING GENES
#THAT ARE VERY TISSUE SPECIFIC IN EXPRESSION AND THUS HAVE LOW MEDIAN
#EXPRESSION OVERALL



#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
if (Sys.getenv("TERM") == "xterm-256color") library("colorout")

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

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(lnc))
clin <- clin[z,]

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


#---------------------------------------------------------
#Make matrix for heat-map. Columns will be tissue type, 
#rows will be gene median values 
#---------------------------------------------------------

##LNCRNAs
heatmap_matrix <- as.data.frame(matrix(ncol=length(unique(clin[,6])),nrow=length(unique(rownames(loop_lnc)))))
colnames(heatmap_matrix) <- unique(clin[,6])
rownames(heatmap_matrix) <- unique(rownames(loop_lnc))

tissues <- unique(clin[,6])
genes <- unique(rownames(loop_lnc))

for(i in 1:length(tissues)){
	tis <- tissues[i]
	pats_tis <- which(clin[,6] %in% tis)
	pats_tis <- clin[pats_tis,1]
	pats_exp <- which(colnames(loop_lnc) %in% pats_tis)
	pats_exp <- loop_lnc[,pats_exp]

	for(y in 1:length(genes)){
		k <- which(rownames(pats_exp) %in% genes[y])
		median_k <- median(as.numeric(pats_exp[k,]))
		heatmap_matrix[y, i] <- median_k
	}

}

heatmap_matrix <- floor(heatmap_matrix)

#---------------------------------------------------------
#Set up heatmap - first floor values
#---------------------------------------------------------

#1. Which genes have medians of 0 or 1 in all cancers
sums <- apply(heatmap_matrix, 1, sum)
zeroes <- which(sums==0)
#remove
heatmap_matrix <- heatmap_matrix[-zeroes,]

#2. Which genes have super high medians relative to others?
sums <- apply(heatmap_matrix, 1, sum)
z <- which(sums >20)
#remove for now 
heatmap_matrix <- heatmap_matrix[-z,]

#3. Which ones have medians less than 1?
meds <- apply(heatmap_matrix, 1, median)
z <- which(meds < 1)
heatmap_matrix <- heatmap_matrix[-z,]

#52 lncRNAs 
pdf("52_lncRNAs_heatmap.pdf", pointsize=5)
heatmap.2(as.matrix(heatmap_matrix),key=FALSE, dendrogram= "column", trace="none",col= mypal[1:5])
dev.off()




