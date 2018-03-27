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
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Multiplot function
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

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
cancers <- c("Liver" , "Ovary", "Pancreas", "Kidney")	

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
#high_lncs <- fread("high_lncsmed4top5cancersPCAWG.txt", sep=";")

##PCAWG sig lncRNAs < 0.05 pvalue from high lncs
#sig <- fread("42sig_lncRNACancerAssociations.txt", sep=";")
#sig$canc <- lapply(sig$canc, function(x) unlist(strsplit(x, " "))[1])
#fdr_sig <- sig[fdr <0.1]

allCands <- readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
#allCands = filter(allCands, gene %in% c("NEAT1", "RP11-622A1.2", "GS1-251I9.4", "ZNF503-AS2", "AC009336.24"))
#allCands = allCands[-6,]

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")

#seperate first by "."
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "\\..*"))
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 
lncs[,1] <- apply(lncs[,1], 1, extract) ; 

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

#PCAWG scores
pcawg_scores = readRDS("PCAWG_four_cancers_ranked_Genes_Mar26.rds")

#TCGA scores 
tcga_scores = readRDS("TCGA_four_cancers_scored_byindexMar26.rds")
z = which(tcga_scores$gene %in% pcawg_scores$gene)
tcga_scores = tcga_scores[z,]
all_cancers_scored = rbind(pcawg_scores, tcga_scores)

###Processing#################################################head()

#1. Want to only look at ENSG genes in rna file
#split feature column and extract third component, write function and apply

#2. remove duplicates 
rna <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]

#ALL GENES USED IN PCAWG
genes <- fread("all_genes_used_inRankingAnalysisPCAWG_Mar26.txt", sep=";")

z <- which(rna[,1] %in% genes$x)
rna <- rna[z,] #25125 both pcgs and lncrnas 

rownames(rna) <- rna[,1]
rna <- rna[,-c(1,2)]
rna <- t(rna)
rna <- as.data.frame(rna)
rna$tissue <- ""
rna$patient <- ""
rna$patient <- rownames(rna)

for(i in 1:nrow(rna)){
	z <- which(clin$SAMPID %in% rna$patient[i])
	rna$tissue[i] <- clin[z,6]
}

z <- which(rna$tissue %in% c(""))

#------------------------------------------------------------------
#Within each tissue type, score genes 
#------------------------------------------------------------------

#1. log1p
rna[,1:(ncol(rna)-2)] <- log1p(rna[,1:(ncol(rna)-2)])

#2. Get lncRNA - median within each tissue type
tissues <- unique(rna$tissue)

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- rna[rna$tissue==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 
getScores <- function(row){
	score=""
	expression <- data.frame(exp=as.numeric(row[1:(length(row)-2)]), gene=names(row)[1:(length(row)-2)])
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
	patients$patient <- rep(names, each=20) #20 lncRNA candidates 
	patients <- as.data.frame(patients)
	patients$canc <- dataframe$tissue[1]
	patients$data <- "GTEX"
	return(patients)
}	

scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_tissues_scored <-  rbindlist(scored)

#bind all tissues scored from GTEX with all cancers scored from pcawg
allScoredLncsgtexANDpcawg <- rbind(all_cancers_scored, all_tissues_scored)
allScoredLncsgtexANDpcawg <- as.data.frame(allScoredLncsgtexANDpcawg)
allScoredLncsgtexANDpcawg$canc <- as.character(allScoredLncsgtexANDpcawg$canc)

allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Ovarian"] = "Ovary"
allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Pancreatic"] = "Pancreas"

allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Kidney"] = "kirc"
allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Liver"] = "lihc"
allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Ovary"] = "ov"
allScoredLncsgtexANDpcawg$canc[allScoredLncsgtexANDpcawg$canc == "Pancreas"] = "paad"

#----------------------------------------------------------------------------------------------------------------
#PLOTS 
#----------------------------------------------------------------------------------------------------------------
#PLOT1 - for each cancer type, show the distribution of ranks between PCAWG, TCGA and GTEX with wilcoxon p-value
#----------------------------------------------------------------------------------------------------------------

tissues <- unique(allScoredLncsgtexANDpcawg$canc) #keep only candidates pancreas, kidney, liver and ovary 
allCands = subset(allCands, allCands$gene %in% allScoredLncsgtexANDpcawg$gene)


pdf("plot1_Version3_pcawgVSgtex_6candidatesWithinEachCancer_Mar27.pdf", pointsize=8)
for(i in 1:length(unique(allCands$gene))){
	allScored <- allScoredLncsgtexANDpcawg[allScoredLncsgtexANDpcawg$gene == allCands$gene[i],]
	cancer = allCands$Cancer[i]
	allScored = subset(allScored, canc == cancer)

	allScored$HR = ""
	for(k in 1:length(unique(allScored$gene))){
		lnc = unique(allScored$gene)[k]
		z = which((allScored$data == "PCAWG") & (allScored$gene == lnc))
		allScored$HR[z] = allCands$PCAWG_HR[which(allCands$gene == lnc)]
		z = which((allScored$data == "TCGA") & (allScored$gene == lnc))
		allScored$HR[z] = allCands$TCGA_HR[which(allCands$gene == lnc)]
	}
	allScored$HR[allScored$HR >= 1] = "Hazardous"
	allScored$HR[allScored$HR < 1] = "Protective"
	allScored$HR[allScored$data == "GTEX"] = "Normal"
	
	for(k in 1:nrow(allScored)){
		lnc = allCands$name[which(allCands$gene == allScored$gene[k])]
		allScored$gene[k] = lnc
	}

	my_comparisons <- list( c("PCAWG", "TCGA"), c("TCGA", "GTEX"), c("PCAWG", "GTEX") )
	f <- ggboxplot(allScored, x="data", y="score", fill="HR", palette=mypal, short.panel.labs=FALSE)
	f = f + stat_compare_means(comparisons = my_comparisons, label.y = c(1.05, 1.12, 1.16), label = "p.signif") 
	
	f <- ggpar(f, xlab="Candidate lncRNAs", main= paste(unique(allScored$gene), "Candidate Gene in", allScored$canc[1]) ,ylab="Score",
	 x.text.angle=65, font.tickslab=c(10, "plain", "black"), legend="right", ylim=c(0,1.6))
	print(f)
}
dev.off()

#----------------------------------------------------------------------------------------------------------------
#PLOT2 - add HR whether it's >1 or <1 in addition to whether the lncRNA is expressed more in cancer or normal
#----------------------------------------------------------------------------------------------------------------










