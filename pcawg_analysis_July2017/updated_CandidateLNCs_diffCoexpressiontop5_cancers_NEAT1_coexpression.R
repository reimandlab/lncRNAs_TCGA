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

#Run linear regression and save coefficients and pvalues 
linear_regression <- function(column){
	gene <- colnames(lnc_data)[column]
	cancer <- lnc_data$canc
	lm0 <- lm(lnc_data[,column] ~ 1)
	lm1 <- lm(lnc_data[,column] ~ 1 + JustLnc_dat)
	anov_p <- anova(lm0, lm1)[2,6]
	coef <- lm1$coefficients[2]
	coef_p <- summary(lm1)$coefficients[2,4]
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
lnc <- allCands$gene[1]
canc <- allCands$canc[1]

#get lnc data 
lnc_data <- lnc_rna[,c(which(colnames(lnc_rna) == lnc),5608,5609)]
#subset data by cancer type that that lnc is a candidate in 
lnc_data <- lnc_data[lnc_data$canc == canc, ]

lnc_data <- merge(lnc_data, pcg_rna, by="patient")
lnc_data <- lnc_data[,-19957]
colnames(lnc_data)[3] <- "canc"

#1. Remove PCGs with median E < 5 FPKM 
#get medians of all PCGs
meds <- apply(lnc_data[,4:19956], 2, median)
#names of pcgs with median <5 
low <- names(meds[which(meds <5)]) 
lnc_data <- lnc_data[,-(which(colnames(lnc_data) %in% low))] #6785 pcgs remain
#log all gene expression values
lnc_data[,c(2,4:(ncol(lnc_data)))] <- log1p(lnc_data[,c(2,4:(ncol(lnc_data)))])

#---------------------------------------------------------
#Analysis - Linear Regression Using Logged Values and 
#cancer type as confounder  
#---------------------------------------------------------

JustLnc_dat <- lnc_data[,2]

#Main Analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#PLOTTING LINEAR REGRESSIOn

#pdf("all_coexpressions.pdf")  
#lapply(4:6032, linear_regression_wplot)
#dev.off()

#Obtain regression results 

#Divide lnc



results <- lapply(4:(ncol(lnc_data)), linear_regression)






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

pdf("3312_coexpressedSignificantly_genes.pdf", pointsize=6, height=8, width=9)
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














