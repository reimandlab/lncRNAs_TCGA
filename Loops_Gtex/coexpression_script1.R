###coexpression_script1.R

#Author: Karina_Isaev
#Date_started: July 12th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data

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

#---------------------------------------------------------
#How many pairs of genes remain? 
#---------------------------------------------------------
id_loops <- which(loops$Gene %in% rownames(loop_lnc))
id_loops_2 <- which(loops$Gene %in% rownames(loop_all))

loops_expression <- loops[c(id_loops, id_loops_2),]

dim(loops)
dim(loops_expression)

#how many pairs
s <- as.data.frame(table(loops_expression$ID))
z <- (which(s$Freq>=2)) ; loops_keep <- as.character(s$Var1[z])
z <- which(loops_expression$ID %in% loops_keep)
loops_expression <- loops_expression[z,]
#4,632 loops remain with at least two genes with gene expression data 

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

#---------------------------------------------------------
#Loop Gene Distribution - lncRNAs
#---------------------------------------------------------

#1. REMOVE ALL THOSE THAT ARE 0 ACROSS THE BOARD BUT SAVE THEM
##lncs
sums <- apply(loop_lnc, 1, sum)
z <- which(sums ==0) #12 have 0 expression in all patients 
zero_e <- rownames(loop_lnc)[z] ; write.table(zero_e, file="GTEx_0expressing_lncs.txt", quote=F, row.names=F)
#remove from expression file
loop_lnc <- loop_lnc[-z,]

#2. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 1
meds <- apply(loop_lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 1) ; 
#1858/2009 have median Expression less than 1
#Don't plot them 
meds2 <- meds[-z,]
#151 remain - make histogram 
g <- gghistogram(meds2, x="meds", title="151 lncRNAs with Median Expression Greater Than 1RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g <- ggpar(g, xticks.by = 2)

##greater than 2s
z <- which(meds2$meds < 2) ; 
meds3 <- meds2[-z,]
#73 remain - make histogram 
g2 <- gghistogram(meds3, x="meds", title="73 lncRNAs with Median Expression Greater Than 2RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g2 <- ggpar(g2, xticks.by = 2)

##greater than 3s
z <- which(meds3$meds < 3) ; 
meds4 <- meds3[-z,]
#38 remain - make histogram 
g3 <- gghistogram(meds4, x="meds", title="38 lncRNAs with Median Expression Greater Than 3RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g3 <- ggpar(g3, xticks.by = 2)

##greater than 4s
z <- which(meds4$meds < 4) ; 
meds5 <- meds4[-z,]
#28 remain - make histogram 
g4 <- gghistogram(meds5, x="meds", title="28 lncRNAs with Median Expression Greater Than 4RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g4 <- ggpar(g4, xticks.by = 2)

##greater than 5s
z <- which(meds5$meds < 5) ; 
meds6 <- meds5[-z,]
#20 remain - make histogram 
g5 <- gghistogram(meds6, x="meds", title="20 lncRNAs with Median Expression Greater Than 5RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g5 <- ggpar(g5, xticks.by = 2)


pdf("GTEx_lncRNAs_arranged_median_summaries.pdf", pointsize=3, width=15, height=12)
grid.arrange(g, g2, g3, g4, g5, ncol = 1, nrow=5)
dev.off()

#---------------------------------------------------------
#Loop Gene Distribution - PCGs
#---------------------------------------------------------

#1. REMOVE ALL THOSE THAT ARE 0 ACROSS THE BOARD BUT SAVE THEM
##pcgs
sums <- apply(loop_all, 1, sum)
z <- which(sums ==0) #13 have 0 expression in all patients 
zero_e <- rownames(loop_all)[z] ; write.table(zero_e, file="GTEx_0expressing_pcgs.txt", quote=F, row.names=F)
#remove from expression file
loop_all <- loop_all[-z,]

#2. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 1
meds <- apply(loop_all, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
z <- which(meds$meds < 1) ; 
#1817/5251 have median Expression less than 1
#Don't plot them 
meds2 <- meds[-z,]
g <- gghistogram(meds2, x="meds", title="3,434 PCGs with Median Expression Greater Than 1RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of PCGs")
g <- ggpar(g, xticks.by = 50)

##greater than 2s
z <- which(meds2$meds < 2) ; 
meds3 <- meds2[-z,]
g2 <- gghistogram(meds3, x="meds", title="2,994 PCGs with Median Expression Greater Than 2RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of PCGs")
g2 <- ggpar(g2, xticks.by = 50)

##greater than 3s
z <- which(meds3$meds < 3) ; 
meds4 <- meds3[-z,]
g3 <- gghistogram(meds4, x="meds", title="2,610 PCGs with Median Expression Greater Than 3RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of PCGs")
g3 <- ggpar(g3, xticks.by = 50)

##greater than 4s
z <- which(meds4$meds < 4) ; 
meds5 <- meds4[-z,]
#28 remain - make histogram 
g4 <- gghistogram(meds5, x="meds", title="2,258 PCGs with Median Expression Greater Than 4RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of PCGs")
g4 <- ggpar(g4, xticks.by = 50)

##greater than 5s
z <- which(meds5$meds < 5) ; 
meds6 <- meds5[-z,]
g5 <- gghistogram(meds6, x="meds", title="1,978 PCGs with Median Expression Greater Than 5RPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of PCGs")
g5 <- ggpar(g5, xticks.by = 50)


pdf("GTEx_PCGs_arranged_median_summaries.pdf", pointsize=3, width=15, height=15)
grid.arrange(g, g2, g3, g4, g5, ncol = 1, nrow=5)
dev.off()

