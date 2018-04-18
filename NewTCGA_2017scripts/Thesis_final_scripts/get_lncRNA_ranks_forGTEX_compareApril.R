###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)

rna = as.data.frame(rna)
dim(rna)
dim(pcg)
dim(norm)
dim(met)

gtex_genes = fread("genes_used_GTExApril17.txt", header=F, data.table=F)

all_genes = merge(rna, pcg, by="patient")
z = which(colnames(all_genes) %in% gtex_genes$V2)
all_genes = all_genes[,c(1,z)]
rownames(all_genes) = all_genes$patient
all_genes$patient = NULL

###---------------------------------------------------------------
###Add scores 
###---------------------------------------------------------------

#1. log1p
all_genes = log1p(all_genes)

#2. Get lncRNA - median within each tissue type
tissues <- unique(rna$Cancer)

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	z = which(rna$Cancer == tissue)
	pats = rna$patient[z]
	z = which(rownames(all_genes) %in% pats)
	tis = all_genes[z,]
	tis = as.data.frame(tis)
	tis$tis = tissue
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 

tcga_genes = readRDS("all_genes_used_in_TCGA_april17.rds")

getScores <- function(row){
	score=""
	expression <- data.frame(exp=as.numeric(row[1:(length(row)-1)]), gene=names(row)[1:(length(row)-1)])
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lncrnas
	lncs = tcga_genes$gene[tcga_genes$type == "lncRNA"]
	z <- which(expression$gene %in% lncs)
	expression <- expression[z, ]
	return(expression)
}

addScores <- function(dataframe){
	patients <- apply(dataframe, 1, getScores) #list of dataframes, need to coerce together
	names <- rownames(dataframe)
	patients <- rbindlist(patients)
	patients = as.data.frame(patients)
	patients$patient <- rep(names, each=length(unique(patients$gene))) 
	patients <- as.data.frame(patients)
	patients$canc <- dataframe$tissue[1]
	patients$data <- "TCGA"
	patients$tis = dataframe$tis[1]
	return(patients)
}	

scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_tissues_scored <-  rbindlist(scored)

#make boxplot of variation of lncRNA score for each patient within a cancer type 

for(i in 1:length(unique(all_tissues_scored$tis))){
	dat = subset(all_tissues_scored, tis == unique(all_tissues_scored$tis)[i])
	dat = dat[200:10000,]
	print(ggscatter(dat, x= "exp", "score", title= unique(all_tissues_scored$tis)[i]))
}

dev.off()

saveRDS(all_tissues_scored, file="all_lncRNAs_exp_scores_inTCGA_all_tissues_April17.rds")

