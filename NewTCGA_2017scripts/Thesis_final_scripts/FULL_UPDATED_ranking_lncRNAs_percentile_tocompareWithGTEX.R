#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 28th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

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
library(ggthemes)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#List of canddidates and cox results
#allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
#List of canddidates and cox results
allCands <- readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
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
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   

#Combined into one dataframe because need to get ranks 
all <- merge(lnc_rna, pcg_rna, by = c("canc", "time", "status", "sex", "patient"))

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
all[,6:ncol(all)] <- log1p(all[,6:ncol(all)])

#2. Get lncRNA - median within each tissue type
tissues <- unique(all$canc)
tissues <- tissues[c(7,9,12,13)]

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- all[all$canc==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 
getScores <- function(row){
	score=""
	expression <- data.frame(exp=as.numeric(row[6:length(row)]), gene=names(row)[6:length(row)])
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lnc candidates - we just want their score 
	z <- which(expression$gene %in% allCands$gene)
	expression <- expression[z, ]
	return(expression)
}


addScores <- function(dataframe){
	patients <- apply(dataframe, 1, getScores) #list of dataframes, need to coerce together
	names <- rownames(dataframe)
	patients <- rbindlist(patients)
	patients$patient <- rep(names, each=25) #25 lncRNA candidates 
	patients <- as.data.frame(patients)
	patients$canc <- dataframe$canc[1]
	patients$data <- "TCGA"
	patients$canc <- lapply(patients$canc, function(x) unlist(strsplit(x, " "))[1])
	return(patients)
}	

scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_cancers_scored <-  rbindlist(scored)
all_cancers_scored <- as.data.frame(all_cancers_scored)

#one for each tissue/cancer type
#each gene is scored within each patient 
#can now make violin plot showing distributon of scores for each candidate lncRNA 
#just need to subset to genes interested in plotting 

saveRDS(all_cancers_scored, file="TCGA_four_cancers_scored_byindexMar26.rds")

#save list of genes in total used to also compare with GTEX 
write.table(colnames(all), file="all_genes_used_inRankingAnalysisTCGA_Mar26.txt", quote=F, row.names=F, sep=";")





























