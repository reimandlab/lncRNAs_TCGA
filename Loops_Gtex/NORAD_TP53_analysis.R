#norad_analysisGTEX.R

###GTEX_42_candidatelncs_limma_diffCoexpression.R

#Author: Karina_Isaev
#Date_started: September 6ths 2017
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
library(limma)

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
high_lncs <- fread("high_lncsmed4top5cancersPCAWG.txt", sep=";")

##PCAWG sig lncRNAs < 0.05 pvalue from high lncs
sig <- fread("42sig_lncRNACancerAssociations.txt", sep=";")
sig$canc <- lapply(sig$canc, function(x) unlist(strsplit(x, " "))[1])
fdr_sig <- sig[fdr <0.1]

allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")

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

#ALL GENES USED IN PCAWG
genes <- fread("all_genes_used_inRankingAnalysis.txt", sep=";")

z <- which(rna[,2] %in% genes$x)
rna <- rna[z,] #25125 both pcgs and lncrnas 

rownames(rna) <- rna[,2]
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

#need to seperate into lncrnas and pcgs 
lncsPcawgList <- fread("lncsPcawgList.txt", sep=";")
pcgsPcawgList <- fread("pcgsPcawgList.txt", sep=";")

lnc_rna <- rna[,c(which(colnames(rna) %in% lncsPcawgList$x), 25126)]
pcg_rna <- rna[,c(which(colnames(rna) %in% pcgsPcawgList$x), 25126)]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------

#List of canddidates and cox results
allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
allCands$canc <- unlist(lapply(allCands$canc, function(x) unlist(strsplit(x, ' '))[1]))

#[1] - for each lncRNA, divide expression data by patient high or low lncRNA expression
#for each row of allCands

#FUNCTION1 - divide expression into high and low for each lnc within a cancer type 
#Apply this function to every row of allCands
getExpression <- function(row){
	lnc <- row[[1]]
	canc <- row[[5]]
	#First divide patients into high and low based on median expression of lnc
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5488:5489)]
	lncData <- lncData[lncData$tissue==canc,]
	med <- median(as.numeric(lncData[,1]))
	lncData$exp[lncData[,1] >= med] <- 1
	lncData$exp[lncData[,1] < med] <- 0 
	return(lncData)
}

divided <- apply(allCands, 1, getExpression)
print("pass")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  
getPCGS <- function(df){
	all <- merge(df, pcg_rna, by="patient")
	all <- all[,-(ncol(all))]
	colnames(all)[3] <- "canc"
	#col 2 is the lncRNA expression data 
	#cols 1, 3, 4 are patient data 

	#Remove PCGs with median E < 5 FPKM 	
	#get medians of all PCGs
	meds <- apply(all[,5:ncol(all)], 2, median)

	#names of pcgs with median <5 
	low <- names(meds[which(meds <5)]) 
	all <- all[,-(which(colnames(all) %in% low))] 
	#log all gene expression values
	all[,c(2,5:(ncol(all)))] <- log1p(all[,c(2,5:(ncol(all)))])

	return(all)
}

dividedWpcgs <- llply(divided, getPCGS, .progress = "text")
print("pass2")


##ALL SAMPLES#-----------------------------------------------------------------------

#norad data 
d <- dividedWpcgs[[1]]

#TP53 - NORAD 
d <- d[,c(1:4, 5887)]

#Run linear regression and save coefficients and pvalues 

column = 5

linear_regression <- function(column, d){
	gene <- colnames(d)[column]
	cancer <- d$canc
	#(1) across all samples
	lm0 <- lm(d[,column] ~ 1)
	lm1 <- lm(d[,column] ~ 1 + d[,2])
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2]
	coef_p <- summary(lm1)$coefficients[2,4]
	if((anov_p <= 0.05) & (coef_p <= 0.05)) {
		#plot regression
		fit <- lm(d[,column] ~ 1 + d[,2])
		print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  		geom_point() +
  		stat_smooth(method = "lm", col = "red") +
  		labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)), x=colnames(d)[2], y=gene))

		return(c(gene, coef, coef_p, anov_p, "allPatients"))
	}
}

#Just NORAD-TP53 Co-Expression 
pdf("NoradTP53Coexpression97GTEXovary.pdf", width=9, height=8, pointsize=8)
linear_regression(column, d)
dev.off()

#Plot boxplots differences between high and low 
d$exp[d$exp==0] <- "Low"
d$exp[d$exp==1] <- "High"

pdf("NoradTP53Boxplots_GTEX_sept15.pdf", pointsize=8)
p <- ggboxplot(d, x="exp", y="TP53", palette=mypal, fill="exp")
p <- p + stat_compare_means()
ggpar(p, legend = "right", legend.title = "LINC00657 Expression",
	xlab="LINC00657", ylab="TP53 log1p(RPKM) Expression")
dev.off()

























