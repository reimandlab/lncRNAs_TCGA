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

#PCAWG Cancers
cancers <- c("Biliary", "Bladder", "Bone","Breast" , "Cervix" , "CNS", "Colon", "Esophagus",
"Head"  , "Kidney" ,   "Liver" , "Lung" , "Lymph" , "Myeloid" , "Ovary" , "Pancreas",
"Prostate" ,    "Skin" , "Stomach" ,  "Thyroid"  , "Uterus", "Cervix Uteri", "Brain")	

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
### 4,859 samples TOTAL
####-------------------

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)

#loop file
load("merged_loop_and_all_full_length_gene_June_13_UCSC_gene_file.rsav")
loops <- as.data.frame(loop_gene)

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

#2. Check how many IDs match the lncRNA IDs

#UCSC file 
lncs <- as.data.frame(lncs)
z <- which(duplicated(ucsc$hg19.ensGene.name2))
ucsc <- ucsc[-z,] #60,234
#keep antisense, lincRNAs from this file 
z <- which(ucsc$hg19.ensemblSource.sourc %in% c("lincRNA", "antisense"))
ucsc_lncs <- ucsc[z,]
dim(ucsc_lncs) #12,563
z <- which(ucsc$hg19.ensemblSource.sourc %in% "protein_coding")
ucsc_prots <- ucsc[z,]
dim(ucsc_prots) #19,002

#3. remove duplicates 
rna2 <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]
lncs2 <- lncs[! lncs$CAT_geneName %in% unique(lncs[duplicated(lncs$CAT_geneName), "CAT_geneName"]), ]

#4. Now seperate rna file into lncs and non-lncs 
z <- which(rna2[,1] %in% ucsc_lncs$hg19.ensGene.name2)
lncs_expression <- rna2[z,] #12,187 - these aren't yet filtered by being the FANTOM5 list but can later be tagged as such
z <- which(rna2[,1] %in% ucsc_prots$hg19.ensGene.name2)
pcg_expression <- rna2[z,] #17,680

#5. Extract gene that are in loops
#cols to keep from loops file

new_loops <- as.data.frame(matrix(ncol=3)) ; colnames(new_loops) <- c("ID", "Gene", "Side")

cols <- colnames(loops)[c(8,9,13,20)]
loops <- loops[cols]
for(i in 1:nrow(loops)){
	id <- loops[i,1]
	x <- loops[i,2]
	#adding x-side genes
	check <- length(unlist(strsplit(x, ",")))
	if(check >1){
		genes <- unlist(strsplit(x, ","))
		for(s in 1:check){
			row <- c(id, genes[s], "x")
			names(row) <- colnames(new_loops)
			new_loops <- rbind(new_loops, row)
		}
	}
	if(check == 1){
		row <- c(id, x, "x")
		names(row) <- colnames(new_loops)
		new_loops <- rbind(new_loops, row)
	}
	y <- loops[i,3]
	check <- length(unlist(strsplit(y, ",")))
	#adding y-side genes
	if(check >1){
		genes <- unlist(strsplit(y, ","))
		for(s in 1:check){
			row <- c(id, genes[s], "y")
			names(row) <- colnames(new_loops)
			new_loops <- rbind(new_loops, row)
		}
	}
	if(check == 1){
		row <- c(id, y, "y")
		names(row) <- colnames(new_loops)
		new_loops <- rbind(new_loops, row)
	}

}

#remove the first .na line
new_loops <- new_loops[-1,]

#save loop file in processed folder
write.table(new_loops, file="processed_loop_fileKI.txt", quote=F, row.names=F)

#Further Process the two files into whether they are genes in loops or not in loops
saveRDS(lncs_expression, "gtex_lncRNA_expression_12187.rds")
saveRDS(pcg_expression, "gtex_pcg_expression_17680.rds")


