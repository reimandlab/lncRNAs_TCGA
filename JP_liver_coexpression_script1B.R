#---------------------------------------------------------
#JP_liver_coexpression_script1B.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: July 10th 2017
#dir: Thesis/pcawg_liver_JP

#------------
#Description:
#------------

#In this script, using JP n=68 samples
#using 70 Tier1 lncRNAs
#make one wide figure showing all violin plots 
#from left to right with increasing medians
#colour them by whether they are lincRNA or antisense 

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

mypal = pal_npg("nrc", alpha = 0.7)(10)

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

#ucsc gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

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
z <- which(meds >=5)
lncs_tier1 <- lnc[z,] #70
meds <- apply(lncs_tier1, 1, median)
#order rows by increasing median
o <-order(meds)
lncs_tier1 <- lncs_tier1[o,]
lncs_tier2 <- lnc[-z,] #12,564

#3. ORDERED VIOLIN PLOTS FOR TOP 70 LNCS
all_genes <- as.data.frame(matrix(nrow=4760, ncol=3)) 
colnames(all_genes) <- c("gene", "patient", "geneE")
for(i in 1:nrow(lncs_tier1)){
	gene <- rownames(lncs_tier1)[i]
	if(i ==1){
		all_genes[1:68,1] <- gene
		all_genes[1:68,2] <- colnames(lncs_tier1)
		all_genes[1:68,3] <- as.numeric(lncs_tier1[i,1:68])
	}
	if(!(i==1)){
		z <- which(is.na(all_genes[,1]))[1]
		all_genes[z:(z+67),1] <- gene
		all_genes[z:(z+67),2] <- colnames(lncs_tier1)
		all_genes[z:(z+67),3] <- as.numeric(lncs_tier1[i,1:68])
	}
}

all_genes$Type <- ""
for(i in 1:nrow(all_genes)){
	id <- all_genes$gene[i]
	z <- which(ucsc[,6] %in% id)
	all_genes$Type[i] <- ucsc[z,7]
}


#plot 
pdf("not_logged_70lncRNAs_tier1.pdf", pointsize=3, width=14, height=9)
g <- ggviolin(all_genes, x = "gene", y = "geneE", color = "Type", draw_quantiles = 0.5, xlab="Gene", palette=mypal, ylab="FPKM")
ggpar(g,
 font.tickslab = c(8,"bold", "black"),
 xtickslab.rt = 45, ytickslab.rt = 45)
dev.off()

#see some outliers - improve by logging 
all_genes$geneE <- log1p(all_genes$geneE)
pdf("70lncRNAs_tier1.pdf", pointsize=3, width=14, height=9)
g <- ggviolin(all_genes, x = "gene", y = "geneE", color = "Type", draw_quantiles = 0.5, xlab="Gene", palette=mypal, ylab="log1p(FPKM)")
ggpar(g,
 font.tickslab = c(8,"bold", "black"),
 xtickslab.rt = 45, ytickslab.rt = 45)
dev.off()
#---------------------------------------------------------
#Save final litsts of genes in each tier to use further
#---------------------------------------------------------
tier1_lncs <- rownames(lncs_tier1) ; write.table(tier1_lncs, file="tier1_lncs.txt", quote=F, row.names=F)
tier2_lncs <- rownames(lncs_tier2) ; write.table(tier2_lncs, file="tier2_lncs.txt", quote=F, row.names=F)

