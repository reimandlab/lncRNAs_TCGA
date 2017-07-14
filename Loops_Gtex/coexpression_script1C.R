###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 13th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data


#SCRIPT - PLOTTING VIOLIN PLOTS FOR LNCRNAS IN LOOPS WITH MEDIAN EXPRESSION GREATER THAN 2RPKM 


#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Libraries#------------------------------------------------
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
lnc <- readRDS("gtex_lncRNA_expression_12187.rds")
all <- readRDS("gtex_pcg_expression_17680.rds")

#reference genes
ref_genes_lncs <- lnc[,1:2]
ref_genes_all <- all[,1:2]

#change feature column with gene names so that they are the rownames
rownames(lnc) <- lnc[,1] ; lnc <- lnc[,-c(1:2)] #13479 lncRNAs as defined by ensembl "lincRNA", "antisense", "sense_intronic", "sense_overlapping"
rownames(all) <- all[,1] ; all <- all[,-c(1:2)] #18039 PCGs as defined by ensembl "protein_coding"

#loops 
loops <- fread("processed_loop_fileKI.txt", data.table=F)

#ucsc gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#---------------------------------------------------------
#Processing lncRNAs and PCGs 
#---------------------------------------------------------

loops$gene_type <- ""

for(i in 1:nrow(loops)){
	id <- loops$Gene[i]
	z <- which(ucsc[,6] %in% id)
	loops$gene_type[i] <- ucsc[z,7]
}

table(loops$Side, loops$gene_type)

#Loop lncRNAs genes 
z <- which(rownames(lnc) %in% loops$Gene)
#******************************************************************
loop_lnc <- lnc[z,] #2021 lncRNAs in loops with gene expression data 
#******************************************************************

#Loop PCGs genes 
z <- which(rownames(all) %in% loops$Gene)
#******************************************************************
loop_all <- all[z,] #5264 PCGs in loops with gene expression data 
#******************************************************************



#2. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 2
meds <- apply(loop_lnc, 1, median)
z <- which(meds >=2)
lncs_tier1 <- loop_lnc[z,] #73
meds <- apply(lncs_tier1, 1, median)
#order rows by increasing median
o <-order(meds)
lncs_tier1 <- lncs_tier1[o,]
lncs_tier2 <- loop_lnc[-z,] #1,948

#3. ORDERED VIOLIN PLOTS FOR TOP 70 LNCS
all_genes <- as.data.frame(matrix(nrow=624515, ncol=4)) 
colnames(all_genes) <- c("gene", "patient", "geneE", "Type")
for(i in 1:nrow(lncs_tier1)){
	gene <- rownames(lncs_tier1)[i]
	gene_type <- which(ucsc[,6] %in% gene)
	if(i ==1){
		all_genes[1:8555,1] <- gene
		all_genes[1:8555,2] <- colnames(lncs_tier1)
		all_genes[1:8555,3] <- as.numeric(lncs_tier1[i,1:8555])
		all_genes[1:8555,4] <- ucsc[gene_type,7]
	}
	if(!(i==1)){
		z <- which(is.na(all_genes[,1]))[1]
		all_genes[z:(z+8554),1] <- gene
		all_genes[z:(z+8554),2] <- colnames(lncs_tier1)
		all_genes[z:(z+8554),3] <- as.numeric(lncs_tier1[i,1:8555])
		all_genes[z:(z+8554),4] <- ucsc[gene_type,7]
	}
}

#plot 
pdf("not_logged_73lncRNAs_tier1_GTex.pdf", pointsize=3, width=14, height=9)
g <- ggviolin(all_genes, x = "gene", y = "geneE", color = "Type", draw_quantiles = 0.5, xlab="Gene", palette=mypal, ylab="RPKM")
ggpar(g,
 font.tickslab = c(8,"bold", "black"),
 xtickslab.rt = 45, ytickslab.rt = 45)
dev.off()

#see some outliers - improve by logging 
all_genes$geneE <- log1p(all_genes$geneE)
pdf("73lncRNAs_tier1_GTex.pdf", pointsize=3, width=14, height=9)
g <- ggviolin(all_genes, x = "gene", y = "geneE", color = "Type", draw_quantiles = 0.5, xlab="Gene", palette=mypal, ylab="log1p(RPKM)")
ggpar(g,
 font.tickslab = c(8,"bold", "black"),
 xtickslab.rt = 45, ytickslab.rt = 45)
dev.off()