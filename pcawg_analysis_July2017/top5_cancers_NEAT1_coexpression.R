#top5_cancers_NEAT1_coexpression.R

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
library(ggthemes)

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
	gene <- colnames(neat1)[column]
	cancer <- neat1$canc
	lm0 <- lm(neat1[,column] ~ cancer)
	lm1 <- lm(neat1[,column] ~ cancer + neat1_dat)
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[8] 
	coef_p <- summary(lm1)$coefficients[8,4]
	#Make plot 
	print(ggplot(lm1$model, aes_string(x = names(lm1$model)[3], y = names(lm1$model)[1])) + 
  	geom_point() +
  	theme_hc() + 
  	stat_smooth(method = "lm", col = "red") +
  	labs(title = paste(gene, "Adj R2 = ",round(signif(summary(lm1)$adj.r.squared, 5),digits=4),
                     "Intercept =",round(signif(lm1$coef[[1]],5), digits=4),
                     " Slope =",round(signif(lm1$coef[[8]], 5), digits=4),
                     " P =",round(signif(summary(lm1)$coef[8,4], 5),digits=4), x= "Neat1 Expression", y=paste(gene, "Expression"))))
	}

#Run linear regression and save coefficients and pvalues 
linear_regression <- function(column){
	gene <- colnames(neat1)[column]
	cancer <- neat1$canc
	lm0 <- lm(neat1[,column] ~ cancer)
	lm1 <- lm(neat1[,column] ~ cancer + neat1_dat)
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[8] 
	coef_p <- summary(lm1)$coefficients[8,4]
	if((anov_p <= 0.05) & (coef_p <= 0.05)) {
		return(c(gene, coef, coef_p, anov_p))
	}
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

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#---------------------------------------------------------
#Pre-Processing - set up NEAT1/PCG matrix for LM
#---------------------------------------------------------
neat1 <- lnc_rna[,c(which(colnames(lnc_rna) == "NEAT1"),5608,5609)]
neat1 <- merge(neat1, pcg_rna, by="patient")
neat1 <- neat1[,-20065]
colnames(neat1)[3] <- "canc"

#1. Remove PCGs with median E < 5 FPKM 
#get medians of all PCGs
meds <- apply(neat1[,4:20062], 2, median)
#names of pcgs with median <5 
low <- names(meds[which(meds <5)]) 
neat1 <- neat1[,-(which(colnames(neat1) %in% low))] #6028 pcgs remain
#log all gene expression values
neat1[,c(2,4:6032)] <- log1p(neat1[,c(2,4:6032)])

#---------------------------------------------------------
#Analysis - Linear Regression Using Logged Values and 
#cancer type as confounder  
#---------------------------------------------------------

neat1_dat <- neat1[,2]

#Main Analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pdf("all_coexpressions.pdf")  
lapply(4:6032, linear_regression_wplot)
dev.off()

#Obtain regression results 
results <- lapply(4:6032, linear_regression)
#remove blank entries - those that weren't significantly associated 
results <- Filter(Negate(is.null), results) 
#convert to df
df <- data.frame(matrix(unlist(results), nrow=3132, byrow=T))
colnames(df) <- c("PCG", "lm_coefficient", "lm_coef_pval", "lm_anov_pval")
df$fdr_coef <- p.adjust(df$lm_coef_pval, method="fdr")
df$fdr_anov <- p.adjust(df$lm_anov_pval, method="fdr")
df$pass <- ""
for(i in 1:nrow(df)){
	coef <- df[i, 5]
	anov <- df[i,6]
	if((coef <=0.05) & (anov <= 0.05)){
		f <- 1
	}
	if(!((coef <=0.05) & (anov <= 0.05))){
		f <- 0
	}
	df$pass[i] <- f
}

#all pass fdr test... 
#but median coefficient of linear regression is 0.077454

#Summarize coefficients/pvalues with volcano plot 
#+++++++++++++++++++++++++++++++++++++++++++
df$logp <- -log10(df$fdr_coef)
df$lm_coefficient <- as.numeric(df$lm_coefficient)

g <- ggscatter(df, x="lm_coefficient", y="logp", palette=mypal, size="logp", color=mypal[3])
g <- g+ labs(title = "3,132 Protein Coding Genes Co-Expressed with Neat1") + 
     theme(plot.title = element_text(hjust = 0.5)) + 
     geom_vline(aes(xintercept=0.1), linetype="dashed") +  geom_vline(aes(xintercept=-0.1), linetype="dashed")

ggpar(g, xlab = "Coefficient of Regression", ylab = "-log10(pval)", legend = "none")
dev.off()


#subset to only those with absolute coefficient size >0.1
df <- subset(df, abs(as.numeric(df$lm_coefficient)) >0.1 ) #900 remain 
df$cor <- ""
for(i in 1:nrow(df)){
	coef <- df[i,2]
	if(coef > 0){
		f <- "pos"
	}
	if(coef < 0){
		f <- "neg"
	}
	df$cor[i] <- f
}

#neg pos                                                                                                                       
#486 414    

write.table(df, file="3132_pcgs_coexpressedWneat1.txt", quote=F, row.names=F)














