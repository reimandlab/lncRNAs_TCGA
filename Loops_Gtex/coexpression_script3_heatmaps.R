###coexpression_script1B.R

#Author: Karina_Isaev
#Date_started: July 14th 2017
#dir: Thesis/GTEx_data/data/processed

#Description:
#Visulizing distribution of gene expression of gene in loops vs those that aren't 
#For lncRNAs and PCGs
#Note: these values are RPKM 
#There are 8,557 unique tissue samples with RNA-Seq data


#SCRIPT - PLOTTING HEATMAP OF GENE MEDIAN EXPRESSION BY TISSUE TYPE
#SEE IF CAN CAPTURE TISSUE SPECIFICTY AND ALSO IF WE ARE LOSING GENES
#THAT ARE VERY TISSUE SPECIFIC IN EXPRESSION AND THUS HAVE LOW MEDIAN
#EXPRESSION OVERALL



#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

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
ensg_hugo <- rbind(ref_genes_all, ref_genes_lncs)
write.table(ensg_hugo, file="esng_to_hugo_all_genes_used_inGTEX_analysis.txt", quote=F, row.names=F)

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

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(lnc))
clin <- clin[z,]

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

#Keep Loop IDs for genes that have expression data
a <- which(loops$Gene %in% rownames(loop_lnc))
b <- which(loops$Gene %in% rownames(loop_all))

#4996 unique loop IDs remain
loops_wexp <- loops[c(a,b),] #7285 genes remain but how many unique loops with at least 2 genes in them? one on each side?
loops_test <- table(loops_wexp$ID, loops_wexp$Side)
rm <- c()

for(i in 1:nrow(loops_test)){
	l <- rownames(loops_test)[i]
	x <- loops_test[i,1]
	y <- loops_test[i,2]
	check <- (x>=1 & y>=1)
	if(!(check)){
		rm <- c(rm, l)
	}
}

#547 to remove
z <- which(loops_wexp$ID %in% rm)
loops_wexp <- loops_wexp[-z,]

#4,449 unique loops with at least one gene on each side that have gene expression data
table(loops_wexp$gene_type)
#write loop file
write.table(loops_wexp, file="4449_loops_wExp_data_GTEx.txt", quote=F, row.names=F)

#---------------------------------------------------------
#Subset gene expression matrix to include lncRNAs/PCGs 
#in those 4,449 loops that have at least one gene on each 
#side with GE data
#---------------------------------------------------------

a <- which(rownames(loop_lnc) %in% loops_wexp$Gene)
#1931/2021 lncRNAs remain
loop_lnc <- loop_lnc[a,]
write.table(loop_lnc, file="1931_lncs_GTEX_inLOOPS.txt", quote=F)


b <- which(rownames(loop_all) %in% loops_wexp$Gene)
#5025/5264 PCGs remain 
loop_all <- loop_all[b,]
write.table(loop_all, file="5025_PCGs_GTEX_inLOOPS.txt", quote=F)

#---------------------------------------------------------
#Make matrix for heat-map. Columns will be tissue type, 
#rows will be gene median values 
#---------------------------------------------------------

##LNCRNAs
heatmap_matrix <- as.data.frame(matrix(ncol=length(unique(clin[,6])),nrow=length(unique(rownames(loop_lnc)))))
colnames(heatmap_matrix) <- unique(clin[,6])
rownames(heatmap_matrix) <- unique(rownames(loop_lnc))

tissues <- unique(clin[,6])
genes <- unique(rownames(loop_lnc))

for(i in 1:length(tissues)){
	tis <- tissues[i]
	pats_tis <- which(clin[,6] %in% tis)
	pats_tis <- clin[pats_tis,1]
	pats_exp <- which(colnames(loop_lnc) %in% pats_tis)
	pats_exp <- loop_lnc[,pats_exp]

	for(y in 1:length(genes)){
		k <- which(rownames(pats_exp) %in% genes[y])
		median_k <- median(as.numeric(pats_exp[k,]))
		heatmap_matrix[y, i] <- median_k
	}

}

heatmap_matrix <- floor(heatmap_matrix)

#---------------------------------------------------------
#Set up heatmap - first floor values
#---------------------------------------------------------

#Cutoff median RPKM - 1

#Which lncRNAs would have been excluded because overall they have median <1?
meds <- apply(loop_lnc, 1, median)
z <- which(meds < 1) #1777 lncRNAs
lncs_med_lessT1 <- names(meds[z])
z <- which(rownames(heatmap_matrix) %in% lncs_med_lessT1)
heatmap_matrix_low <- heatmap_matrix[z,]
heatmap_matrix_high <- heatmap_matrix[-z,] #154

#1. Which genes have medians of 0 or 1 in all cancers
sums <- apply(heatmap_matrix_low, 1, sum)
zeroes <- which(sums==0) #1477 have 0 medians across the board
#remove
heatmap_matrix_low <- heatmap_matrix_low[-zeroes,]
#300 left

#2. Which genes have super high medians relative to others?
sums <- apply(heatmap_matrix_low, 1, sum)
#highest sum of medians = Xist 

my_palette <- colorRampPalette(mypal[5:10])(n = 500)

pdf("300_low_lncRNAs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_low), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D2'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

#---------------------------------------------------------
#Cutoff median RPKM - 2

#Which lncRNAs would have been excluded because overall they have median <2?
meds <- apply(loop_lnc, 1, median)
z <- which(meds < 2) #1858 lncRNAs
lncs_med_lessT1 <- names(meds[z])
z <- which(rownames(heatmap_matrix) %in% lncs_med_lessT1)
heatmap_matrix_low <- heatmap_matrix[z,]
heatmap_matrix_high <- heatmap_matrix[-z,] #73

#1. Which genes have medians of 0 or 1 in all cancers
sums <- apply(heatmap_matrix_low, 1, sum)
zeroes <- which(sums==0) #1477 have 0 medians across the board
#remove
heatmap_matrix_low <- heatmap_matrix_low[-zeroes,]
#381 left

#2. Which genes have super high medians relative to others?
sums <- apply(heatmap_matrix_low, 1, sum)
#highest sum of medians = Xist 

my_palette <- colorRampPalette(mypal[5:10])(n = 500)

pdf("381_medianLessThan2_low_lncRNAs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_low), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D2'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

#to add them to the list of 73, check how many have medians >=2 in at least two tissues 
keep <- c()

for(i in 1:nrow(heatmap_matrix_low)){
	z <- length(which(heatmap_matrix_low[i,] >= 2))
	if(z>=3){
		keep <- c(keep, rownames(heatmap_matrix_low)[i])
	}
	}	

#117 left out of 381 if cutoff is median >=2 in at least two tissues
#84  left out of 381 if cutoff is median >=2 in at least three tissues (84+73 == 157 lncRNAs)

z <- which(rownames(heatmap_matrix_low) %in% keep)
heatmap_matrix_low <- heatmap_matrix_low[z,]

#these 84 have median >=2 in at least three tissues 
pdf("84_highConfidence_low_lncRNAs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_low), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D2'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

my_palette <- colorRampPalette(mypal[c(5:10, 1:2)])(n = 100)
#these 73 have median >=2 when looking at all samples together
pdf("73_highConfidence_high_lncRNAs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_high), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D2'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

#final list of 157 lncRNAS to use
all_lncs <- c(rownames(heatmap_matrix_low), rownames(heatmap_matrix_high))
all_lncs <- as.data.frame(all_lncs)
write.table(all_lncs, file="157_tier1_lncsGTEX.txt", quote=F, row.names=F)

#For final list of lncRNAs (that meet minimum expression requirements either Tier1 or Tier2)
#Check if any of them showed significant association with proliferation from Chang's CRISPRi screen
#Check their expression via Cage 
#Run them through the lncRNA-disease database

#PCA using low vs high expressing lncRNAs



#---------------------------------------------------------
#Make matrix for heat-map. Columns will be tissue type, 
#rows will be gene median values 
#PCGS****************************************************
#---------------------------------------------------------

##PCGs
heatmap_matrix <- as.data.frame(matrix(ncol=length(unique(clin[,6])),nrow=length(unique(rownames(loop_all)))))
colnames(heatmap_matrix) <- unique(clin[,6])
rownames(heatmap_matrix) <- unique(rownames(loop_all))

tissues <- unique(clin[,6])
genes <- unique(rownames(loop_all))

for(i in 1:length(tissues)){
	tis <- tissues[i]
	pats_tis <- which(clin[,6] %in% tis)
	pats_tis <- clin[pats_tis,1]
	pats_exp <- which(colnames(loop_all) %in% pats_tis)
	pats_exp <- loop_all[,pats_exp]

	for(y in 1:length(genes)){
		k <- which(rownames(pats_exp) %in% genes[y])
		median_k <- median(as.numeric(pats_exp[k,]))
		heatmap_matrix[y, i] <- median_k
	}

}
heatmap_matrix <- floor(heatmap_matrix)

#---------------------------------------------------------
#Cutoff median RPKM - 5 PCGS******************************

#Which PCGs would have been excluded because overall they have median <5?
meds <- apply(loop_all, 1, median)
z <- which(meds < 5) # 3106 PCGs
pcgs_med_lessT1 <- names(meds[z])
z <- which(rownames(heatmap_matrix) %in% pcgs_med_lessT1)
heatmap_matrix_low <- heatmap_matrix[z,]

heatmap_matrix_high <- heatmap_matrix[-z,] #1919 PCGs
#check how many outliers, ie super high expressed
medians <- apply(heatmap_matrix_high, 1, median)
z <- which(medians >500)
rm <- which(rownames(heatmap_matrix_high) %in% names(z))
heatmap_matrix_high <- heatmap_matrix_high[-rm,] #1914 PCGs left

#2. Which genes have medians of 0 or 1 in all cancers
sums <- apply(heatmap_matrix_low, 1, sum)
zeroes <- which(sums==0) #672 have 0 medians across the board
#remove
heatmap_matrix_low <- heatmap_matrix_low[-zeroes,]
#2434 left

#to add them to the list of 1914, check how many have medians >=5 in at least two tissues 
keep <- c()

for(i in 1:nrow(heatmap_matrix_low)){
	z <- length(which(heatmap_matrix_low[i,] >= 5))
	if(z>=3){
		keep <- c(keep, rownames(heatmap_matrix_low)[i])
	}
	}	

#726  left out of 2434 if cutoff is median >=5 in at least three tissues (726+1914 == 2640 PCGs)

z <- which(rownames(heatmap_matrix_low) %in% keep)
heatmap_matrix_low <- heatmap_matrix_low[z,]

mypal[5] <- "white"
my_palette <- colorRampPalette(mypal[c(5:10, 1:2)])(n = 4000)
#these 726 have median >=5 RPKM in at least three tissues 
pdf("726_highConfidence_low_PCGs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_low), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

my_palette <- colorRampPalette(mypal[c(5:10, 1:2)])(n = 500)
#these 1914 have median >=5 when looking at all samples together
pdf("1914_highConfidence_high_PCGs_heatmap.pdf", pointsize=2)
heatmap.2(as.matrix(heatmap_matrix_high), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D2'),
	distfun = function(x) dist(x,method = 'euclidean'))
dev.off()

#final list of 2640 PCGs to use
all_pcgs <- c(rownames(heatmap_matrix_low), rownames(heatmap_matrix_high))
all_pcgs <- as.data.frame(all_pcgs)
write.table(all_pcgs, file="2640_tier1_pcgsGTEX.txt", quote=F, row.names=F)


#--------------DONE----------------------------------------------------------------------------------------------------------




