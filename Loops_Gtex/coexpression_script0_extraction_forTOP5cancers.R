###coexpression_script0_extraction.R

#Author: Karina_Isaev
#Date_started: July 11th 2017
#dir: Thesis/GTEx_data/data/raw

#Description:
#Processing RNA-Seq file from GTEx
#Reformating loops file
#Note: these values are RPKM 
#There are 4,859 unique tissue samples with RNA-Seq data

###Preamble###############################################
options(stringsAsFactors=F)

###Libraries##############################################
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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)


###Data#####################################################

#RNA-Seq file 
rna <- readRDS("gtex_expression_file.rds")
rna <- as.data.frame(rna)

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#------------------------------------------------------------------
#Processing clinical file - inlcude samples that have expression
#------------------------------------------------------------------

tissues <- as.data.frame(clin[,c(1,6)])
c <- as.data.frame(table(tissues$SMTS))
c[,1] <- as.character(c[,1])
c[,2] <- as.numeric(c[,2])

#PCAWG Cancers top 5 
cancers <- c("Breast" , "Kidney" ,  "Liver" , "Ovary" , "Pancreas")	

z <- which(c[,1] %in% cancers)
cancers_keep <- c[z,1]

#Filter clinical file
z <- which(clin[,6] %in% cancers_keep)
clin <- clin[z,]

#keep only patients in gene expression file that are part of the tissues
#we want to study 
z <- which(colnames(rna) %in% clin[,1])
rna <- rna[,c(1,2,z)]

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(rna))
clin <- clin[z,]

####-------------------
### 633 samples TOTAL
####-------------------

##PCAWG top 5 high cancers 
high_lncs <- fread("high_lncsmed4.5top5cancersPCAWG.txt", sep=";")

##PCAWG sig lncRNAs < 0.05 pvalue from high lncs
sig <- fread("50sig_lncRNACancerAssociations.txt", sep=";")
sig$canc <- lapply(sig$canc, function(x) unlist(strsplit(x, " "))[1])
fdr_sig <- sig[fdr <0.1]

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")
#subset to those in high_lncs
lncs <- subset(lncs, lncs$CAT_geneName %in% high_lncs$gene)

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

###Processing#################################################head()

#1. Want to only look at ENSG genes in rna file
#split feature column and extract third component, write function and apply

#seperate first by "."
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "\\..*"))
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 
lncs[,1] <- apply(lncs[,1], 1, extract) ; 

#2. remove duplicates 
rna <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]

#3. Subset to only high_lncs expression 
lnc_rna <- subset(rna, rna[,1] %in% lncs$CAT_geneID) #243/244 lncRNAs available 

rownames(lnc_rna) <- lnc_rna[,2]
lnc_rna <- lnc_rna[,-c(1,2)]
lnc_rna <- t(lnc_rna)
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$tissue <- ""
lnc_rna$patient <- ""
lnc_rna$patient <- rownames(lnc_rna)

for(i in 1:nrow(lnc_rna)){
	z <- which(clin$SAMPID %in% lnc_rna$patient[i])
	lnc_rna$tissue[i] <- clin[z,6]
}

#------------------------------------------------------------------
#Plot high lncRNAs per tissue ordered by increasing medians 
#Then compare to the same lncRNAs but in PCAWG 
#------------------------------------------------------------------

lnc_rna[,1:243] <- log1p(lnc_rna[,1:243])

cancers_list <- unique(high_lncs$canc)
#only care about pancreas, ovary, liver and clear cell 
cancers_list <- cancers_list[c(1,4,5,6)]

#Function 1
#input: list of cancers who had high expressing genes 
#output: list of cancer-genes associated with it because they had high expression
get_genes <- function(cancer){
	genes <- high_lncs$gene[high_lncs$canc==cancer]
	return(list(cancer, genes))
}

canc_genes_list <- lapply(cancers_list, get_genes)

#Function 2
#input: list of cancer-genes associated with it 
#output: list of dataframes with only samples associated with that tissue and rna expression of only those genes 

get_rna_data <- function(canc_genes){
	tissue <- unlist(strsplit(canc_genes[[1]], " "))[1]
	genes <- canc_genes[[2]]
	z2 <- which(lnc_rna$tissue == tissue)  
	z <- which(colnames(lnc_rna) %in% genes)
	df <- lnc_rna[z2,c(z, ncol(lnc_rna)-1)]
	return(df)
}

canc_gene_expression <- lapply(canc_genes_list, get_rna_data)

#Function 3
#input: list of datframes
#output: dataframe prepared for plotting violin plots side by side for each gene 

for(i in 1:length(canc_gene_expression)){
	dataframe <- canc_gene_expression[[i]]
	#nrow = #genes * #patients
	rows <- (nrow(dataframe) * (ncol(dataframe)-1))
	df3 <- as.data.frame(matrix(nrow=rows, ncol=4))
	colnames(df3) <- c("gene", "patient", "exp", "PCAWG")

	genes <- colnames(dataframe)[1:(ncol(dataframe)-1)]
	
	meds <- list()
	sig_genes <- sig$gene[which(sig$canc %in% dataframe$tissue[1])]
	fdr_genes <- fdr_sig$gene[which(fdr_sig$canc %in% dataframe$tissue[1])]

	for(j in 1:length(genes)){
		gene <- genes[j]
		z <- which(is.na(df3[["gene"]]))[1]       
		gene_exp <- which(colnames(dataframe) %in% gene)
		gene_exp <- dataframe[,gene_exp]
		h <- length(gene_exp)
		df3[z:(z+h-1),1] <- gene
		df3[z:(z+h-1),2] <- rownames(dataframe) 
		df3[z:(z+h-1),3] <- gene_exp
		sig_tag <- length(which(sig_genes %in% gene))
		fdr_tag <- length(which(fdr_genes %in% gene))
		if(sig_tag == 0){
			sig_tag <- "Not Sig"
		}
		if(sig_tag == 1){
			sig_tag <- "Sig"
		}
		if(fdr_tag==1){
			sig_tag <- "FDR Sig"
		}
		df3[z:(z+h-1),4] <- sig_tag
		med <- median(gene_exp)
		meds[[j]] <- c(gene, med)
	}
	#get list of ordered genes
	meds <- as.data.table(matrix(unlist(meds), ncol=2, byrow=T))	
	meds <- meds[order(V2)]
	order <- meds$V1
	tissue <- dataframe$tissue[1]
	pdf(paste(tissue, "high_lncs_expression.pdf", sep="_"), pointsize=9, width=16, height=12)
	g <- ggboxplot(df3, x= "gene", y="exp", add="boxplot", order=order, palette=mypal, fill="PCAWG")
	g <- ggpar(g, main= paste(tissue, "High Expressing lncRNAs in Cancer,", nrow(dataframe), "GTEx Samples") ,xlab = "lncRNA", ylab = "log1p(RPKM)", xtickslab.rt = 65, font.tickslab = c(7,"plain", "black"))
	g <- g + geom_hline(aes(yintercept=log1p(5)), colour="#990000")
	g <- g + geom_hline(aes(yintercept=log1p(10)), colour="#990000")
	print(g) 
	dev.off()
}

















