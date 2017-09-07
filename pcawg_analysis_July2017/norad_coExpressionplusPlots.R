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
	stage <- clin$tumour_stage[which(clin$icgc_donor_id %in% pat)]
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
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5608:5609)]
	lncData <- lncData[lncData$canc==canc,]
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

#FUNCTION3 - run linear regression within each dataframe comparing lncRNA and each PCG
#(1) across all samples
#(2) just in high and just in low to compare gene lists obtained between high and low groups

#> which(colnames(d) %in% "PUM1")
#[1] 4549
#> which(colnames(d) %in% "PUM2")
#[1] 4550


##ALL SAMPLES#-----------------------------------------------------------------------

#norad data 
d <- dividedWpcgs[[1]]

#Run linear regression and save coefficients and pvalues 
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

#(1) across all samples 
pdf("norad_coexpressedPCGSsept6th.pdf", pointsize=8)
d <- as.data.frame(d)
results <- llply(5:ncol(d), linear_regression, d=d, .progress = "text")
#remove blank entries - those that weren't significantly associated 
results <- Filter(Negate(is.null), results) 
#convert to df
df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
colnames(df) <- c("PCG", "lm_coefficient", "lm_coef_pval", "lm_anov_pval", "Patients")
df$fdr_coef <- p.adjust(df$lm_coef_pval, method="fdr")
df$fdr_anov <- p.adjust(df$lm_anov_pval, method="fdr")
df$pass <- ""
for(i in 1:nrow(df)){
	coef <- df[i, 6]
	anov <- df[i,7]
	if((coef <=0.05) & (anov <= 0.05)){
			f <- 1
		}
	if(!((coef <=0.05) & (anov <= 0.05))){
			f <- 0
		}
		df$pass[i] <- f
	}
df$canc <- d$canc[1]
df$lnc <- colnames(d)[2]
all <- df
dev.off()

#Rank genes and figure out what genes do 
#Rank first by adjusted pvalue 
all <- as.data.table(all)
all <- all[order(as.numeric(fdr_coef))]

#Then by absolute coefficient slope/correlation value 
all <- all[order(-(abs(as.numeric(lm_coefficient))))]
#pathway enrichment looking at top 100 coexpressed genes based on the above ranking
genes <- all[1:100,1]
allnoradpaths <- gprofiler(genes, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr", include_graph=T)

##HIGH NORAD SAMPLES#------------------------------------------------------------------

#norad data 
d <- dividedWpcgs[[1]]

linear_regression_highOnly <- function(column, d){
	gene <- colnames(d)[column]
	cancer <- d$canc
	dataframe2 <- d[d$exp ==1, ]
	#(2) across only high expressing lncRNA samples 
	lm0 <- lm(dataframe2[,column] ~ 1)
	lm1 <- lm(dataframe2[,column] ~ 1 + dataframe2[,2])
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2]
	coef_p <- summary(lm1)$coefficients[2,4]
	if((anov_p <= 0.05) & (coef_p <= 0.05)) {
		#plot regression
		fit <- lm(dataframe2[,column] ~ 1 + dataframe2[,2])
		print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  		geom_point() +
  		stat_smooth(method = "lm", col = "red") +
  		labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)), x=colnames(d)[2], y=gene))

		return(c(gene, coef, coef_p, anov_p, "HIGHPatients"))
	}
}

#(2) across high lncRNA samples 
pdf("norad_HIGHgroup_coexpressedPCGSsept6th.pdf", pointsize=8)
d <- as.data.frame(d)
results <- llply(5:ncol(d), linear_regression_highOnly, d=d, .progress = "text")
#remove blank entries - those that weren't significantly associated 
results <- Filter(Negate(is.null), results) 
#convert to df
df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
colnames(df) <- c("PCG", "lm_coefficient", "lm_coef_pval", "lm_anov_pval", "Patients")
df$fdr_coef <- p.adjust(df$lm_coef_pval, method="fdr")
df$fdr_anov <- p.adjust(df$lm_anov_pval, method="fdr")
df$pass <- ""
for(i in 1:nrow(df)){
	coef <- df[i, 6]
	anov <- df[i,7]
	if((coef <=0.05) & (anov <= 0.05)){
			f <- 1
		}
	if(!((coef <=0.05) & (anov <= 0.05))){
			f <- 0
		}
		df$pass[i] <- f
	}
df$canc <- d$canc[1]
df$lnc <- colnames(d)[2]
high <- df
dev.off()

#Rank genes and figure out what genes do 
#Rank first by adjusted pvalue 
high <- as.data.table(high)
high <- high[order(as.numeric(fdr_coef))]

#Then by absolute coefficient slope/correlation value 
high <- high[order(-(abs(as.numeric(lm_coefficient))))]
genes <- high[1:100,1]
highnoradpaths <- gprofiler(genes, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr", include_graph=T)

##LOW NORAD SAMPLES#------------------------------------------------------------------

#norad data 
d <- dividedWpcgs[[1]]

linear_regression_lowOnly <-  function(column, d){
	gene <- colnames(d)[column]
	cancer <- d$canc
	dataframe3 <- d[d$exp ==0, ]
	#(3) across only low expressing lncRNA samples 
	lm0 <- lm(dataframe3[,column] ~ 1)
	lm1 <- lm(dataframe3[,column] ~ 1 + dataframe3[,2])
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2]
	coef_p <- summary(lm1)$coefficients[2,4]
	if((anov_p <= 0.05) & (coef_p <= 0.05)) {
		#plot regression
		fit <- lm(dataframe3[,column] ~ 1 + dataframe3[,2])
		print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  		geom_point() +
  		stat_smooth(method = "lm", col = "red") +
  		labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)), x=colnames(d)[2], y=gene))

		return(c(gene, coef, coef_p, anov_p, "LOWPatients"))
	}
}

#(2) across high lncRNA samples 
pdf("norad_LOWgroup_coexpressedPCGSsept6th.pdf", pointsize=8)
d <- as.data.frame(d)
results <- llply(5:ncol(d), linear_regression_lowOnly, d=d, .progress = "text")
#remove blank entries - those that weren't significantly associated 
results <- Filter(Negate(is.null), results) 
#convert to df
df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
colnames(df) <- c("PCG", "lm_coefficient", "lm_coef_pval", "lm_anov_pval", "Patients")
df$fdr_coef <- p.adjust(df$lm_coef_pval, method="fdr")
df$fdr_anov <- p.adjust(df$lm_anov_pval, method="fdr")
df$pass <- ""
for(i in 1:nrow(df)){
	coef <- df[i, 6]
	anov <- df[i,7]
	if((coef <=0.05) & (anov <= 0.05)){
			f <- 1
		}
	if(!((coef <=0.05) & (anov <= 0.05))){
			f <- 0
		}
		df$pass[i] <- f
	}
df$canc <- d$canc[1]
df$lnc <- colnames(d)[2]
low <- df
dev.off()

#Rank genes and figure out what genes do 
#Rank first by adjusted pvalue 
low <- as.data.table(low)
low <- low[order(as.numeric(fdr_coef))]

#Then by absolute coefficient slope/correlation value 
low <- low[order(-(abs(as.numeric(lm_coefficient))))]
genes <- low[1:100,1]
lownoradpaths <- gprofiler(genes, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr", include_graph=T)
