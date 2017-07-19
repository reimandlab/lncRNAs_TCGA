###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 18th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Using cutoffs established previosuly for minimum median expression 
#and using genes that are in looops with genes expressed on both sides
#conduct co-expression analysis using logged linear regression 
#to identify PCGs co-expressed with each lncRNA at a signficant value

#next:
#for significant results 
#plot coefficient of regression versus adjusted pvalue volcano plot 

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

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

#lncRNA Expression (lncs in loops)
lnc_exp <- read.table("1931_lncs_GTEX_inLOOPS.txt", header=T, row.names=1)

#PCGs Expression (PCGs in loops)
pcg_exp <- read.table("5025_PCGs_GTEX_inLOOPS.txt", header=T, row.names=1)

#List of Tier1 lncRNAs
lncs <- fread("157_tier1_lncsGTEX.txt", data.table=F)

#subset expression file to these lncs
z <- which(rownames(lnc_exp) %in% lncs[,1])


#List of Tier1 PCGs
pcgs <- fread("2640_tier1_pcgsGTEX.txt", data.table=F)

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(lnc_exp))
clin <- clin[z,]
#keep only sample ID and tissue type for now 

#Loops
loops <- fread("4449_loops_wExp_data_GTEx.txt", data.table=F)

#Ensg to Hugo Conversion
conversion <- fread("esng_to_hugo_all_genes_used_inGTEX_analysis.txt", data.table=F)
colnames(conversion)[1] <- "Gene"

loops <- merge(loops, conversion, by="Gene")


##---------------------------------------------------------
#Co-Expression Analysis
#---------------------------------------------------------

#PARAMETERS REQUIRED

#1. lncRNA Gene Expression Matrix, rows = genes, columns = samples
#2. PCG Gene Expression Matrix, rows = genes, columns = samples
#3. List of each type of genes that are in loops with at least one gene on each side with median 
	#expression greater than the cutoff 
	#(2RPKM for lncRNAs & 5RPKM for PCGS) 
#4. Sample ID - Tissue Type two columns matrix 
















































