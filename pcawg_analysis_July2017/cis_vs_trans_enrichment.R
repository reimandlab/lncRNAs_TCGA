#42candidate_lncs_limma_diffCoexpression .R

#Karina Isaev
#Based on the list of differentially expressed genes 
#for some reason they are expressed differently between two groups of patients split
#based on the expression of a given lncRNA 
#if that lncRNA is responsible for regulating those genes 
#in one group when lncRNA is highly expressed some genes also get more expressed 
#so lncRNA might be activating them or their upstream regulators 
#how many of those genes are cis versus trans 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
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
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(parallel)
library(limma)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data#-------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
#z <- which(duplicated(ucsc[,6]))
#ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "ID"

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

#Right now just looking at the top 6 in PCAWG and TCGA 
cancer_pairs = data.frame(cancer = c(rep(unique(allCands$canc)[1], 3), rep(unique(allCands$canc)[2], 2), unique(allCands$canc)[3]), genes = c("LINC00665", 
	"ZNF503-AS2", "GS1-251I9.4", "NEAT1", "ADORA2A-AS1", "AC006126.4"))

allCands = subset(allCands, (allCands$canc %in% cancer_pairs$canc & allCands$gene %in% cancer_pairs$genes))
allCands = allCands[-5,]

#---------------------------------------------------------
#Get coordinates of DE genes 
#---------------------------------------------------------

DE_genes = fread("GSI_DE_genes.txt", data.table=F)
DE_genes = merge(DE_genes, ucsc, by = "ID")
colnames(DE_genes)[8:14] = c("Transcript", "Chr", "Strand", "Tstart", "Tend", "ENSG", "Type")

DE_genes = subset(DE_genes, abs(DE_genes$logFC) >= 0.2)

DE_genes_bedfile = DE_genes[,c(9, 11,12,1,2,10,3,6,13)]
write.table(DE_genes_bedfile, file="GSI_ovary_DE_genes.bed", quote=F, col.names=F,row.names=F, sep="\t")

#---------------------------------------------------------
#Get coordinates of lncRNA (GSI in this case)
#---------------------------------------------------------

gsi = allCands[3,]
colnames(gsi)[1] = "ID"
gsi = merge(gsi, ucsc, by = "ID")
colnames(gsi)[7:13] = c("Transcript", "Chr", "Strand", "Tstart", "Tend", "ENSG", "Type")
gsi = gsi[,c(8, 10,11,1,3,9,12)]
write.table(gsi, file="GSI_coordinates.bed", quote=F, col.names=F,row.names=F, sep="\t")

#---------------------------------------------------------
#Hypergeometric test - is there an enrichment of DE genes
#on the same chromosome?
#---------------------------------------------------------

DE_genes = fread("GSI_DE_genes.txt", data.table=F)
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

colnames(ucsc)[8] = "ID"

DE_genes = merge(DE_genes, ucsc, by = "ID")
colnames(DE_genes)[8:14] = c("Transcript", "Chr", "Strand", "Tstart", "Tend", "ENSG", "Type")
DE_genes = subset(DE_genes, abs(DE_genes$logFC) >= 0.2)

#pop size : 6786 (total genes tested)
#sample size : 1112 (DE expressed genes)
#Number of items in the pop that are classified as successes : 248 (248 genes from chromosome 8)
#Number of items in the sample that are classified as successes : 98 

#P(Observed 98 or more) = 1-P(Observed less than 98)
#1.0-phyper(62-1, 1998, 5260-1998, 131)

1 - (phyper((98-1),248,(6786-248),1112))











































