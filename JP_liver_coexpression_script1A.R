#---------------------------------------------------------
#JP_liver_coexpression_script1.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: June 26th 2017
#dir: Thesis/pcawg_liver_JP

#------------
#Description:
#------------

#Have two RNA-Seq files, one for lncRNAs (all from UCSC including low confidence) 
#and one for PCGs. Conduct here LM and NGB regression to identify 
#signficantly co-expressed lncRNA-PCGs in liver cancer patients 
#using an array of confounders 
#------------
#confounders:
#------------
#[1]. PCG CNA
#[2]. lncRNA CNA
#[3]. Clinical features 
#[4]. is methylation available?

#SCRIPT1 - PROCESS GENES, VISUALIZE EXPRESSION OF GENES AND 
#SET UP DATA FRAME FOR REGRESSION ANALYSIS 

#---------------------------------------------------------
#Preamble
#---------------------------------------------------------

options(stringsAsFactors=F)

#---------------------------------------------------------
#Libraries
#---------------------------------------------------------
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

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#These are gene expression files based on the original RNA-Seq
#PCAWG file
lnc <- readRDS("liver_jp_lncRNA_expression_6028.rds")
all <- readRDS("liver_jp_pcg_expression.rds")

#change feature column with gene names so that they are the rownames
rownames(lnc) <- lnc[,1] ; lnc <- lnc[,-1] #13479 lncRNAs as defined by ensembl "lincRNA", "antisense", "sense_intronic", "sense_overlapping"
rownames(all) <- all[,1] ; all <- all[,-1] #18039 PCGs as defined by ensembl "protein_coding"

#---------------------------------------------------------
#Processing lncRNAs and PCGs 
#---------------------------------------------------------

#1. REMOVE ALL THOSE THAT ARE 0 ACROSS THE BOARD BUT SAVE THEM
##lncs
sums <- apply(lnc, 1, sum)
z <- which(sums ==0) #845 have 0 expression in all patients 
zero_e <- rownames(lnc)[z] ; write.table(zero_e, file="jp_liver_0expressing_lncs.txt", quote=F, row.names=F)
#remove from expression file
lnc <- lnc[-z,]
##pcgs
sums <- apply(all, 1, sum)
z <- which(sums ==0) #845 have 0 expression in all patients 
zero_e <- rownames(all)[z] ; write.table(zero_e, file="jp_liver_0expressing_pcgs.txt", quote=F, row.names=F)
#remove from expression file
all <- all[-z,]

#2. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 1
meds <- apply(lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
#12,634 lncRNAs
z <- which(meds$meds < 1) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds2 <- meds[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan1.pdf", pointsize=4, width=8.5)
g <- gghistogram(meds2, x="meds", title="723 lncRNAs with Median Expression Greater Than 1FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g <- ggpar(g, xticks.by = 2)
#dev.off()

##greater than 2s
z <- which(meds2$meds < 2) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds3 <- meds2[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g2 <- gghistogram(meds3, x="meds", title="317 lncRNAs with Median Expression Greater Than 2FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g2 <- ggpar(g2, xticks.by = 2)
#dev.off()

##greater than 3s
z <- which(meds3$meds < 3) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds4 <- meds3[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g3 <- gghistogram(meds4, x="meds", title="169 lncRNAs with Median Expression Greater Than 3FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g3 <- ggpar(g3, xticks.by = 2)
#dev.off()

##greater than 4s
z <- which(meds4$meds < 4) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds5 <- meds4[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g4 <- gghistogram(meds5, x="meds", title="105 lncRNAs with Median Expression Greater Than 4FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g4 <- ggpar(g4, xticks.by = 2)
#dev.off()

##greater than 5s
z <- which(meds5$meds < 5) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds6 <- meds5[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g5 <- gghistogram(meds6, x="meds", title="70 lncRNAs with Median Expression Greater Than 5FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g5 <- ggpar(g5, xticks.by = 2)
#dev.off()


pdf("arranged_median_summaries.pdf", pointsize=3, width=10, height=12)
grid.arrange(g, g2, g3, g4, g5, ncol = 1, nrow=5)
dev.off()