###coexpression_script0_extraction.R

#Author: Karina_Isaev
#Date_started: July 11th 2017
#dir: Thesis/GTEx_data/data/raw

#Description:
#Processing RNA-Seq file from GTEx
#Reformating loops file
#Note: these values are RPKM 
#There are 4,859 unique tissue samples with RNA-Seq data

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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)


###Data#####################################################

#RNA-Seq file 
rna <- readRDS("gtex_expression_file.rds")
rna <- as.data.frame(rna)

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#------------------------------------------------------------------
#Processing clinical file - inlcude samples that have expression
#------------------------------------------------------------------

tissues <- as.data.frame(clin[,c(1,6)])
c <- as.data.frame(table(tissues$SMTS))
c[,1] <- as.character(c[,1])
c[,2] <- as.numeric(c[,2])

#PCAWG Cancers top 5 
cancers <- c("Breast" , "Kidney" ,  "Liver" , "Ovary" , "Pancreas")	

z <- which(c[,1] %in% cancers)
cancers_keep <- c[z,1]

#Filter clinical file
z <- which(clin[,6] %in% cancers_keep)
clin <- clin[z,]

#keep only patients in gene expression file that are part of the tissues
#we want to study 
z <- which(colnames(rna) %in% clin[,1])
rna <- rna[,c(1,2,z)]

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(rna))
clin <- clin[z,]

####-------------------
### 633 samples TOTAL
####-------------------

##PCAWG top 5 high cancers 
high_lncs <- fread("high_lncsmed4.5top5cancersPCAWG.txt", sep=";")

##PCAWG sig lncRNAs < 0.05 pvalue from high lncs
sig <- fread("50sig_lncRNACancerAssociations.txt", sep=";")
sig$canc <- lapply(sig$canc, function(x) unlist(strsplit(x, " "))[1])
fdr_sig <- sig[fdr <0.1]

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

###Processing#################################################head()

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

#2. remove duplicates 
rna <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]

#3. subset to lncRNAs based on FANTOM annotation 
lnc_rna <- rna[which(rna[,1] %in% lncs$CAT_geneID),] #5835 

rownames(lnc_rna) <- lnc_rna[,2]
lnc_rna <- lnc_rna[,-c(1,2)]
lnc_rna <- t(lnc_rna)
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$tissue <- ""
lnc_rna$patient <- ""
lnc_rna$patient <- rownames(lnc_rna)

for(i in 1:nrow(lnc_rna)){
	z <- which(clin$SAMPID %in% lnc_rna$patient[i])
	lnc_rna$tissue[i] <- clin[z,6]
}

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
lnc_rna[,1:5835] <- log1p(lnc_rna[,1:5835])

#2. Get lncRNA - median within each tissue type
tissues <- unique(lnc_rna$tissue)

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- lnc_rna[lnc_rna$tissue==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA-RNAseq data 
#output: lncrna-Median ranked by percentile 

get_gene_medians <- function(dataframe){
	tis <- dataframe$tissue[1]
	meds <- apply(dataframe[,1:(dim(dataframe)[2]-2)], 2, median)
	bg_by_bin = by(order(meds), ceiling(seq_along(1:length(meds))/(length(meds)/100)), function(x) meds[x])
	#need to get canc-gene-percentile 
	data2 <- data.frame(Gene="",Percent = "", Tissue="")     
	for(i in 1:length(bg_by_bin)){
		percentile <- i
		genes <- names(bg_by_bin[[i]])
		data <- data.frame(Gene=genes,Percent = percentile, Tissue=tis)     
		data2 <- rbind(data2, data)
	}
	data2 <- data2[-1,]
	return(data2)
}

binned_lncs <- lapply(tissues_data, get_gene_medians)

#Function 3
#input: dataframe with lncrna-percentile-tissue
#output: print dataframe so can compare with PCAWG

print_ranked <- function(dataframe){
	filename <- paste(dataframe$Tissue[1], "ranked_by_bins_lncRNAs_GTEX.txt", sep="_")
	write.table(dataframe, file=filename, quote=F, row.names=F, sep=";")
}	
lapply(binned_lncs, print_ranked)































