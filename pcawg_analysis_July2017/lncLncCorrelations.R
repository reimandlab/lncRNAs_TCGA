#42candidate_lncs_limma_diffCoexpression .R

#Karina Isaev
#September 5th, 2017

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
library(limma)
library(corrplot)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data#-------------------------------------------------

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

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

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

#FUNCTION2 - divide by cancer type get correlations of candidate lncRNAs in that cancer type
values <- unique(allCands$canc)
correlateLNCS <- function(value){
	lncs <- filter(allCands, canc==value)
	lncsData <- lnc_rna[which(lnc_rna$canc %in% value),]
	lncsData <- lncsData[, which(colnames(lncsData) %in% lncs$gene)]
	lncsData <- log1p(lncsData)
	cor.mtest <- function(mat, conf.level = 0.95){
  	mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
            p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
            lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
            uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
        }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
	}
	res1 <- cor.mtest(lncsData,0.95)
	M <- cor(lncsData)
	## specialized the insignificant value according to the significant level
	print(corrplot(M, method="circle", type="upper", order="hclust", p.mat = res1[[1]], sig.level=0.05, insig = "pch"))
}

pdf("NewcandidateCancersWithCorrelationsofLNCRNAs_Sept13th.pdf", pointsize=8, width=8, height=9)
llply(values, correlateLNCS)
dev.off()

#FUNCTION3 - calculate linear regression between top 
ovarian <- values[[1]]

lncs <- filter(allCands, canc==ovarian)
lncsData <- lnc_rna[which(lnc_rna$canc %in% ovarian),]
lncsData <- lncsData[, which(colnames(lncsData) %in% lncs$gene)]

#Run linear regression and save coefficients and pvalues 
linear_regression <- function(column, d){
	if(!(column == 8)){
	gene <- colnames(d)[column]
	#(1) across all samples
	lm0 <- lm(d[,column] ~ 1)
	lm1 <- lm(d[,column] ~ 1 + d[,8]) #8 because LINC00657 is in column 8 
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2]
	coef_p <- summary(lm1)$coefficients[2,4]
	#if((anov_p <= 0.05) & (coef_p <= 0.05)) {
		#plot regression
		fit <- lm(d[,column] ~ 1 + d[,8])
		print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  		geom_point() +
  		stat_smooth(method = "lm", col = "red") +
  		labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)), x=colnames(d)[8], y=gene))

		return(c(gene, coef, coef_p, anov_p, "allPatients"))
	#}
	}
}

#(1) across all samples 
d <- lncsData
pdf("logged_ovary_LNCLNC_coexpressedPCGSsept13th.pdf", pointsize=8)
d <- as.data.frame(d)
d <- log1p(d)
results <- llply(1:ncol(d), linear_regression, d=d, .progress = "text")

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
df$lnc <- colnames(d)[8]
dev.off()









































