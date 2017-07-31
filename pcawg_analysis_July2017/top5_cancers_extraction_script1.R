#top5_cancers_extraction_script1.R

#Karina Isaev
#July 20th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script1 - using the top 5 cancer types chosen
#PLOTS:

#1. Heatmap of median lncRNA expression in a cancer type 

#2. PCA to show that lncRNA expression in general
#can seperate cancer types 

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

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

#Cancers to use 
tum_types <- fread("top5_cancers_andHISTO_to_keepJuly20.txt", data.table=F)

#---------------------------------------------------------
#Processing
#---------------------------------------------------------

###"NORMAL SAMPLES"
z <- which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))
norm_pats <- conversion$icgc_donor_id[z]

###"TUMOUR SAMPLES" - PROCESSING RNA FILE 
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z]

for(i in 1:ncol(rna)){
	z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna)[i])
	if(!(length(z)==0)){
		colnames(rna)[i] <- conversion$icgc_donor_id[z]
	}
}

extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1:2], 1, ensg)
z <- which(check==TRUE)
rna <- rna[z,]

#2. Check how many IDs match the lncRNA IDs
#none match while trancript number is present 
#remove ie, ENSG00000201285.1 --> ENSG00000201285
#in both rna file and lncs file

extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract3) ; 

#Remove duplicate genes 
z <- which(duplicated(rna[,1]))
genes <- rna[z,1]
z <- which(rna[,1] %in% genes)
rna <- rna[-z,]

#Using UCSC keep only antisense, lincRNA and protein-coding genes
z <- which(rna[,1] %in% ucsc$hg19.ensGene.name2)
rna <- rna[z,]
 
#rows
rownames(rna) <- rna[,1]
rna <- rna[,-1]

##Keep only patient samples that also have clinical data
z <- which(colnames(rna) %in% clin$icgc_donor_id)
rna <- rna[,z]

#Divide RNA into lnc_RNA and pcg_RNA
z <- which(ucsc$hg19.ensemblSource.source == "protein_coding")
zz <- which(rownames(rna) %in% ucsc$hg19.ensGene.name2[z])

###################
###PCG_RNA#########
###################
pcg_rna <- rna[zz,]


###################
###LNC_RNA#########
###################
lnc_rna <- rna[-zz,]


###SUBSET CLINICAL AND EXPRESSION FILE TO ONLY THE TOP 5 CANCERS/HISTOS

##CLIN
clin_top5 <- subset(clin, clin$icgc_donor_id %in% tum_pats)
z <- which(duplicated(clin_top5$icgc_donor_id))
clin_top5 <- clin_top5[-z,]
clin_top5 <- subset(clin_top5, (clin_top5$histology_tier2 %in% tum_types$V1) & (clin_top5$histology_tier4 %in% tum_types$V2))
clin_top5$combined_tum_histo <- ""
clin_top5$combined_tum_histo <- paste(clin_top5[,15], clin_top5[,17])

#EXPRESSION - pcgs
z <- which(colnames(pcg_rna) %in% clin_top5$icgc_donor_id)
pcg_rna_top5 <- pcg_rna[,z] #20166

#EXPRESSION - lncs
z <- which(colnames(lnc_rna) %in% clin_top5$icgc_donor_id)
lnc_rna_top5 <- lnc_rna[,z] #12598
pca_lncs <- lnc_rna_top5

#---------------------------------------------------------
#lncRNA Median Expression Cancer Heatmap 
#---------------------------------------------------------

##LNCRNAs
heatmap_matrix <- as.data.frame(matrix(ncol=dim(tum_types)[1],nrow=length(unique(rownames(lnc_rna_top5)))))
colnames(heatmap_matrix) <-  paste(tum_types[,1], tum_types[,2])
rownames(heatmap_matrix) <- unique(rownames(lnc_rna_top5))

tissues <- paste(tum_types[,1], tum_types[,2])
genes <- unique(rownames(lnc_rna_top5))

for(i in 1:length(tissues)){
	tis <- tissues[i]
	pats_tis <- which(clin_top5$combined_tum_histo %in% tis)
	pats_tis <- clin_top5[pats_tis,9]
	pats_exp <- which(colnames(lnc_rna_top5) %in% pats_tis)
	pats_exp <- lnc_rna_top5[,pats_exp]

	for(y in 1:length(genes)){
		k <- which(rownames(pats_exp) %in% genes[y])
		median_k <- median(as.numeric(pats_exp[k,]))
		heatmap_matrix[y, i] <- median_k
	}
}

##Are there any lncs that median of 0 in all cancers?
sums <- apply(heatmap_matrix, 1, sum) #3658
s <- which(sums == 0)
z <- which(rownames(heatmap_matrix) %in% names(s))
heatmap_matrix <- heatmap_matrix[-z,] #8940 genes left 

floored_heatmap_matrix <- floor(heatmap_matrix)
##now how many genes have a median of 0 in all tissues?
sums <- apply(floored_heatmap_matrix, 1, sum) #7328
s <- which(sums == 0)
z <- which(rownames(floored_heatmap_matrix) %in% names(s))
floored_heatmap_matrix <- floored_heatmap_matrix[-z,] #1612 genes left 
#that have median of at least 1 FPKM in at least one cancer type 

##some outliers
maxs <- apply(floored_heatmap_matrix, 1, max)
summary(maxs)
   # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   # 1.0      1.0      2.0    127.5      3.0 147231.0 

#cap it at 100, only 8 lncRNAs with max median E greater 15
s <- which(maxs >100)
z <- which(rownames(floored_heatmap_matrix) %in% names(s))
floored_heatmap_matrix <- floored_heatmap_matrix[-z,] #1604 genes left 

mypal[5] <- "white"
my_palette <- colorRampPalette(mypal[c(5,8,3,2)])(n = 150)

#pdf("top5cancers_lncRNAs_heatmapJuuly20.pdf", pointsize=2)
heatmap.2(as.matrix(floored_heatmap_matrix), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D'),
	distfun = function(x) dist(x,method = 'euclidean'), srtCol=35, cexCol=0.9, key.title=NA, keysize=1, labRow = FALSE,
          margins = c(9, 2), main = list("Median Expression of 1,604 lncRNAs Across Cancers", cex = 1.2))

dev.off()

#what if divide the lncRNAs into high and low median groups to see heatmap more clearly 

#divide by medians

meds <- apply(floored_heatmap_matrix, 1, median)
#median of medians is 0
z <- which(meds==0)
s <- names(z)
z <- which(rownames(floored_heatmap_matrix) %in% s)
lncs1 <- floored_heatmap_matrix[z,]
#medians of median >=1
lncs2 <- floored_heatmap_matrix[-z,]

pdf("top5cancers_1073medianofmedians0_lncRNAs_heatmapJuly31.pdf", pointsize=2, width=12)
heatmap.2(as.matrix(lncs1), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D'),
	distfun = function(x) dist(x,method = 'euclidean'), srtCol=35, cexCol=0.9, key.title=NA, keysize=1, labRow = FALSE,
          margins = c(9, 2), main = list("Median Expression of 1,073 lncRNAs Across Cancers", cex = 1.2))
dev.off()

pdf("top5cancers_531medianofmedianshigherthan0_lncRNAs_heatmapJuly31.pdf", pointsize=2, width=12)
heatmap.2(as.matrix(lncs2), trace="none",col= my_palette, hclustfun = function(x) hclust(x,method = 'ward.D'),
	distfun = function(x) dist(x,method = 'euclidean'), srtCol=35, cexCol=0.9, key.title=NA, keysize=1, labRow = FALSE,
          margins = c(9, 2), main = list("Median Expression of 531 lncRNAs Across Cancers", cex = 1.2))
dev.off()



































