###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 13th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data


#SCRIPT - PLOTTING NUMBER OF LOOPS REMAINING AFTER CHANGING MEDIAN CUTOFF



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


#Results Storage -------------------------------------------------------
results <- as.data.frame(matrix(ncol=3)) ; colnames(results) <- c("median","Count", "Type")


#------------------------------------------------------------------------
#Based on different median RPKM cutoffs how many pairs remain - MEDIAN 1
#------------------------------------------------------------------------
#1. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 1
#lncs*********************************************
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 1) ; 
meds2 <- meds[-z,]
lncs_med_1 <- meds2$gene
#pcgs*********************************************
meds <- apply(loop_all, 1, median)
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 1) ; 
meds2 <- meds[-z,]
pcgs_med_1 <- meds2$gene
#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% lncs_med_1)
id_loops_2 <- which(loops$Gene %in% pcgs_med_1)
loops_expression <- loops[c(id_loops, id_loops_2),]
dim(loops)
dim(loops_expression)
#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#2,000 loops remain with at least two genes with gene expression data (both genes have median RPKM >=1)
#*********************************************
median1 <- length(unique(loops_expression$ID))
#---------------------------------------------------------
#How many lnc-lnc, lnc-pcg loops remain? 
#---------------------------------------------------------
loop_types <- table(loops_expression$ID, loops_expression$gene_type)
loop_types <- as.matrix(loop_types)
loop_type_list <- c()
for(i in 1:nrow(loop_types)){
	counts_lnc <- sum(as.numeric(loop_types[i,1:2]))
	counts_pcg <- sum(as.numeric(loop_types[i,3]))
	if(counts_lnc > 0 & counts_pcg >0 ){
		loop_type_list <- c(loop_type_list, "lnc-pcg")
	}
	if(counts_lnc > 1 & counts_pcg == 0 ){
		loop_type_list <- c(loop_type_list, "lnc-lnc")
	}
	if(counts_lnc == 0 & counts_pcg >1 ){
		loop_type_list <- c(loop_type_list, "pcg-pcg")
	}
}
loop_types <- cbind(loop_types, loop_type_list)
table(loop_types[,4])
#row indicating how many total loops, lnc-lnc, lnc-pcg loops remain

row <- c("1RPKM", table(loop_types[,4])[2],  "lnc-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("1RPKM", table(loop_types[,4])[3], "PCG-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("1RPKM", table(loop_types[,4])[1], "lnc-lnc")
names(row) <- colnames(results)
results <- rbind(results, row)

#*********************************************

#1 loop with both genes lncRNAs above 1 median L9037
#Genes in loop:
loops$Gene[which(loops$ID %in% "L9037")]
# "ENSG00000245532" "ENSG00000251562"



#------------------------------------------------------------------------
#Based on different median RPKM cutoffs how many pairs remain - MEDIAN 2
#------------------------------------------------------------------------
#1. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 2
#lncs*********************************************
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 2) ; 
meds2 <- meds[-z,]
lncs_med_1 <- meds2$gene
#pcgs*********************************************
meds <- apply(loop_all, 1, median)
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 2) ; 
meds2 <- meds[-z,]
pcgs_med_1 <- meds2$gene
#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% lncs_med_1)
id_loops_2 <- which(loops$Gene %in% pcgs_med_1)
loops_expression <- loops[c(id_loops, id_loops_2),]
dim(loops)
dim(loops_expression)
#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#1601 loops remain with at least two genes with gene expression data (both genes have median RPKM >=1)
#*********************************************
median2 <- length(unique(loops_expression$ID))
#---------------------------------------------------------
#How many lnc-lnc, lnc-pcg loops remain? 
#---------------------------------------------------------
loop_types <- table(loops_expression$ID, loops_expression$gene_type)
loop_types <- as.matrix(loop_types)
loop_type_list <- c()
for(i in 1:nrow(loop_types)){
	counts_lnc <- sum(as.numeric(loop_types[i,1:2]))
	counts_pcg <- sum(as.numeric(loop_types[i,3]))
	if(counts_lnc > 0 & counts_pcg >0 ){
		loop_type_list <- c(loop_type_list, "lnc-pcg")
	}
	if(counts_lnc > 1 & counts_pcg == 0 ){
		loop_type_list <- c(loop_type_list, "lnc-lnc")
	}
	if(counts_lnc == 0 & counts_pcg >1 ){
		loop_type_list <- c(loop_type_list, "pcg-pcg")
	}
}
loop_types <- cbind(loop_types, loop_type_list)
table(loop_types[,4])
#row indicating how many total loops, lnc-lnc, lnc-pcg loops remain
row <- c("2RPKM", table(loop_types[,4])[2],  "lnc-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("2RPKM", table(loop_types[,4])[3], "PCG-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("2RPKM", table(loop_types[,4])[1], "lnc-lnc")
names(row) <- colnames(results)
results <- rbind(results, row)
#*********************************************





#------------------------------------------------------------------------
#Based on different median RPKM cutoffs how many pairs remain - MEDIAN 3
#------------------------------------------------------------------------
#1. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 3
#lncs*********************************************
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 3) ; 
meds2 <- meds[-z,]
lncs_med_1 <- meds2$gene
#pcgs*********************************************
meds <- apply(loop_all, 1, median)
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 3) ; 
meds2 <- meds[-z,]
pcgs_med_1 <- meds2$gene
#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% lncs_med_1)
id_loops_2 <- which(loops$Gene %in% pcgs_med_1)
loops_expression <- loops[c(id_loops, id_loops_2),]
dim(loops)
dim(loops_expression)
#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#1,277 loops remain with at least two genes with gene expression data (both genes have median RPKM >=1)
#*********************************************
median3 <- length(unique(loops_expression$ID))
#---------------------------------------------------------
#How many lnc-lnc, lnc-pcg loops remain? 
#---------------------------------------------------------
loop_types <- table(loops_expression$ID, loops_expression$gene_type)
loop_types <- as.matrix(loop_types)
loop_type_list <- c()
for(i in 1:nrow(loop_types)){
	counts_lnc <- sum(as.numeric(loop_types[i,1:2]))
	counts_pcg <- sum(as.numeric(loop_types[i,3]))
	if(counts_lnc > 0 & counts_pcg >0 ){
		loop_type_list <- c(loop_type_list, "lnc-pcg")
	}
	if(counts_lnc > 1 & counts_pcg == 0 ){
		loop_type_list <- c(loop_type_list, "lnc-lnc")
	}
	if(counts_lnc == 0 & counts_pcg >1 ){
		loop_type_list <- c(loop_type_list, "pcg-pcg")
	}
}
loop_types <- cbind(loop_types, loop_type_list)
table(loop_types[,4])
#row indicating how many total loops, lnc-lnc, lnc-pcg loops remain
row <- c("3RPKM", table(loop_types[,4])[2],  "lnc-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("3RPKM", table(loop_types[,4])[3], "PCG-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("3RPKM", table(loop_types[,4])[1], "lnc-lnc")
names(row) <- colnames(results)
results <- rbind(results, row)
#*********************************************





#------------------------------------------------------------------------
#Based on different median RPKM cutoffs how many pairs remain - MEDIAN 4
#------------------------------------------------------------------------
#1. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 4
#lncs*********************************************
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 4) ; 
meds2 <- meds[-z,]
lncs_med_1 <- meds2$gene
#pcgs*********************************************
meds <- apply(loop_all, 1, median)
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 4) ; 
meds2 <- meds[-z,]
pcgs_med_1 <- meds2$gene
#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% lncs_med_1)
id_loops_2 <- which(loops$Gene %in% pcgs_med_1)
loops_expression <- loops[c(id_loops, id_loops_2),]
dim(loops)
dim(loops_expression)
#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#1,024 loops remain with at least two genes with gene expression data (both genes have median RPKM >=1)
#*********************************************
median4 <- length(unique(loops_expression$ID))
#---------------------------------------------------------
#How many lnc-lnc, lnc-pcg loops remain? 
#---------------------------------------------------------
loop_types <- table(loops_expression$ID, loops_expression$gene_type)
loop_types <- as.matrix(loop_types)
loop_type_list <- c()
for(i in 1:nrow(loop_types)){
	counts_lnc <- sum(as.numeric(loop_types[i,1:2]))
	counts_pcg <- sum(as.numeric(loop_types[i,3]))
	if(counts_lnc > 0 & counts_pcg >0 ){
		loop_type_list <- c(loop_type_list, "lnc-pcg")
	}
	if(counts_lnc > 1 & counts_pcg == 0 ){
		loop_type_list <- c(loop_type_list, "lnc-lnc")
	}
	if(counts_lnc == 0 & counts_pcg >1 ){
		loop_type_list <- c(loop_type_list, "pcg-pcg")
	}
}
loop_types <- cbind(loop_types, loop_type_list)
table(loop_types[,4])
#row indicating how many total loops, lnc-lnc, lnc-pcg loops remain
row <- c("4RPKM", table(loop_types[,4])[2],  "lnc-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("4RPKM", table(loop_types[,4])[3], "PCG-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("4RPKM", table(loop_types[,4])[1], "lnc-lnc")
names(row) <- colnames(results)
results <- rbind(results, row)
#*********************************************





#------------------------------------------------------------------------
#Based on different median RPKM cutoffs how many pairs remain - MEDIAN 5
#------------------------------------------------------------------------
#1. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 5
#lncs*********************************************
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 5) ; 
meds2 <- meds[-z,]
lncs_med_1 <- meds2$gene
#pcgs*********************************************
meds <- apply(loop_all, 1, median)
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 5) ; 
meds2 <- meds[-z,]
pcgs_med_1 <- meds2$gene
#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% lncs_med_1)
id_loops_2 <- which(loops$Gene %in% pcgs_med_1)
loops_expression <- loops[c(id_loops, id_loops_2),]
dim(loops)
dim(loops_expression)
#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#866 loops remain with at least two genes with gene expression data (both genes have median RPKM >=1)
#*********************************************
median5 <- length(unique(loops_expression$ID))
#---------------------------------------------------------
#How many lnc-lnc, lnc-pcg loops remain? 
#---------------------------------------------------------
loop_types <- table(loops_expression$ID, loops_expression$gene_type)
loop_types <- as.matrix(loop_types)
loop_type_list <- c()
for(i in 1:nrow(loop_types)){
	counts_lnc <- sum(as.numeric(loop_types[i,1:2]))
	counts_pcg <- sum(as.numeric(loop_types[i,3]))
	if(counts_lnc > 0 & counts_pcg >0 ){
		loop_type_list <- c(loop_type_list, "lnc-pcg")
	}
	if(counts_lnc > 1 & counts_pcg == 0 ){
		loop_type_list <- c(loop_type_list, "lnc-lnc")
	}
	if(counts_lnc == 0 & counts_pcg >1 ){
		loop_type_list <- c(loop_type_list, "pcg-pcg")
	}
}
loop_types <- cbind(loop_types, loop_type_list)
table(loop_types[,4])
#row indicating how many total loops, lnc-lnc, lnc-pcg loops remain
row <- c("5RPKM", table(loop_types[,4])[2],  "lnc-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("5RPKM", table(loop_types[,4])[3], "PCG-PCG")
names(row) <- colnames(results)
results <- rbind(results, row)

row <- c("5RPKM", table(loop_types[,4])[1], "lnc-lnc")
names(row) <- colnames(results)
results <- rbind(results, row)
#*********************************************

results <- results[-1,]
results$Count <- as.numeric(results$Count)


#------------------------------------------------------------------------
#PLOT
#------------------------------------------------------------------------

pdf("GTEx_loops_median_cutoffs_barplot_KI_july14th.pdf", pointsize=4, height=12, width=13)
ggbarplot(results, "median", "Count",
  fill = "Type", color = "Type", palette = "Paired",
 xlab="Median RPKM Cutoff", ylab="Number of Loops That Remain",
 label = FALSE) + 

ggtitle("Loop Distribution as Median Cutoff Requirement Changes - GTEx") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()












