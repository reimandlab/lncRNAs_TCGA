#JOB1_ALL_PATIENTS_COEXPRESSION.R

#Karina Isaev
#August 10th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run co-expression analysis with cancer type and sex? as confounders
#to get list of co-expressed PCGs with NEAT1 to run through m:Explorer  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
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
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(parallel)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Functions-------------------------------------------------

#Plotting linear regression line 
ggplotRegression <- function (fit) {
require(ggplot2)
ggplot(fit$model, aes_string(x = names(fit$model)[3], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[8]], 5),
                     " P =",signif(summary(fit)$coef[8,4], 5)))
}


#Run linear regression and obtain p-values and generate plot
linear_regression_wplot <- function(column){
	gene <- colnames(lnc_data)[column]
	cancer <- lnc_data$canc
	lm0 <- lm(lnc_data[,column] ~ 1)
	lm1 <- lm(lnc_data[,column] ~ 1 + JustLnc_dat)
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2] 
	coef_p <- summary(lm1)$coefficients[2,4]
	#Make plot 
	print(ggplot(lm1$model, aes_string(x = names(lm1$model)[3], y = names(lm1$model)[1])) + 
  	geom_point() +
  	theme_hc() + 
  	stat_smooth(method = "lm", col = "red") +
  	labs(title = paste(gene, "Adj R2 = ",round(signif(summary(lm1)$adj.r.squared, 5),digits=4),
                     "Intercept =",round(signif(lm1$coef[[1]],5), digits=4),
                     " Slope =",round(signif(lm1$coef[[2]], 5), digits=4),
                     " P =",round(signif(summary(lm1)$coef[2,4], 5),digits=4), x= "Neat1 Expression", y=paste(gene, "Expression"))))
	}



#Data-------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
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
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#Add stage data 
addStage <- function(row){
	pat <- row[[length(row)]]
	stage <- clin$tumour_stage[which(clin$icgc_donor_id %in% pat)[1]]
	return(stage)
}

lnc_rna$stage <- apply(lnc_rna, 1, addStage)

#Add time data 
addStage <- function(row){
	pat <- row[[length(row)]]
	stage <- clin$tumour_stage[which(clin$icgc_donor_id %in% pat)[1]]
	return(stage)
}

lnc_rna$stage <- apply(lnc_rna, 1, addStage)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)
pcg_rna$stage <- apply(pcg_rna, 1, addStage)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin_pats <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin_pats$icgc_donor_id %in% rownames(lnc_rna))
clin_pats <- clin_pats[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin_pats$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin_pats$icgc_donor_id),] #485 patients remain 

#make sure there aren't any lncRNAs in protein coding file 
table(ucsc[,7][ucsc[,8] %in% colnames(pcg_rna)])
#Remove 
z <- which(colnames(pcg_rna) %in% fantom[,2])
pcg_rna <- pcg_rna[,-z]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------

#List of canddidates and cox results
allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")

#[1] - for each lncRNA, divide expression data by patient high or low lncRNA expression
#for each row of allCands

#FUNCTION1 - divide expression into high and low for each lnc within a cancer type 
#Apply this function to every row of allCands
getExpression <- function(row){
	lnc <- row[[1]]
	canc <- row[[5]]
	#First divide patients into high and low based on median expression of lnc
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5608:5610)]
	lncData <- lncData[lncData$canc==canc,]
	med <- median(as.numeric(lncData[,1]))
	lncData$exp[lncData[,1] >= med] <- 1
	lncData$exp[lncData[,1] < med] <- 0 
	return(lncData)
}

divided <- apply(allCands, 1, getExpression)
print("pass")

#In the above list of dataframes matches order of genes in allCands file
lnc1 <- divided[[1]] ; lnc1$lnc <- colnames(lnc1)[1] ; colnames(lnc1)[5] <- paste("lnc1", "tag", sep="")
lnc3 <- divided[[2]] ; lnc3$lnc <- colnames(lnc3)[1] ; colnames(lnc3)[5] <- paste("lnc3", "tag", sep="")
colnames(lnc1)[1] <-"lncExp" ; colnames(lnc3)[1] <- "lncExp"
both <- as.data.table(merge(lnc1, lnc3, by=c("canc", "patient")))

#Get survival analysis 
res.cut <- surv_cutpoint(both, time = "time", event = "event",
   variables = c("DEPDC1", "WHSC1", "CRIM1"))



