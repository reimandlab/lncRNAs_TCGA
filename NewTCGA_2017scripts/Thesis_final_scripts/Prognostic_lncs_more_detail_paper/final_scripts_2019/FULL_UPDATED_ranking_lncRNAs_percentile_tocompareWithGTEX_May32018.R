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
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
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
#allCands <- readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")
#cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")

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

#Combined into one dataframe because need to get ranks 
#all <- merge(rna, pcg, by = c("patient", "Cancer"))
#rownames(all) = all$patient
#all = all[,1:25170]

#genes <- fread("all_genes_used_inRankingAnalysisPCAWG_Mar26.txt", sep=";")
#rownames(all) = all$patient
#all$patient = NULL
#z = which(colnames(all) %in% genes$x)
#all = all[,c(1, z)]

tcga_genes = fread("all_genes_used_inRankingAnalysisTCGA_May4th.txt")
tcga_genes$type = ""
tcga_genes$type[tcga_genes$x %in% fantom$CAT_geneID] = "lncRNA"
tcga_genes$type[(tcga_genes$type == "")] = "pcg"

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
z = which(str_detect(colnames(all), "ENSG"))	
all[,z] <- log1p(all[,z])

#2. Get lncRNA - median within each tissue type
#tissues <- tissues[c(7,9,12,13)]

t = as.data.table(table(all$type))
t = filter(t, N >=50)
all = subset(all, type %in% t$V1)
tissues <- unique(all$Cancer)
print(table(all$type))

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- all[all$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 

addScores <- function(dtt){
	names <- as.list(dtt$patient)
	getScores <- function(patient){
	z = which((dtt$patient) %in% patient)
	row = dtt[z,]
	z = which(str_detect(names(row), "ENSG"))

	score=""
	expression <- data.frame(exp=as.numeric(row[z]), gene=names(row)[z])
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lncrnas
	lncs = tcga_genes$x[tcga_genes$type == "lncRNA"]
	z <- which(expression$gene %in% lncs)
	expression <- expression[z, ]
	expression = as.data.frame(expression)
	expression$patient = ""
	expression$patient = patient 
	return(expression)
		}
	patients <- llply(names, getScores, .progress="text") #list of dataframes, need to coerce together
	#names <- rownames(dataframe)
	patients1 <- rbindlist(patients)
	patients1 <- as.data.frame(patients1)
	patients1$tis <- dtt$Cancer[1]
	patients1$data <- "TCGA"
	#patients1$tis = dataframe$tis[1]
	return(patients1)
}	

scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_tissues_scored <-  rbindlist(scored)

all_cancers_scored <- as.data.frame(all_tissues_scored)

#one for each tissue/cancer type
#each gene is scored within each patient 
#can now make violin plot showing distributon of scores for each candidate lncRNA 
#just need to subset to genes interested in plotting 
all_cancers_scored$exp = as.numeric(all_cancers_scored$exp)

saveRDS(all_cancers_scored, file="TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")

#save list of genes in total used to also compare with GTEX 
write.table(colnames(all), file="all_genes_used_inRankingAnalysisTCGA_May4th.txt", quote=F, row.names=F, sep=";")





























