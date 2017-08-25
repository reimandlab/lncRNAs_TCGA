#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 11th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
library(ggthemes)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

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
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

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
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 


#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------

#For each patient add survival status and days since last seen 
lnc_rna$status <- ""
lnc_rna$time <- ""
lnc_rna$sex <- ""

#lncs
for(i in 1:nrow(lnc_rna)){
  pat <- rownames(lnc_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  lnc_rna$status[i] <- clin$donor_vital_status[z]
  lnc_rna$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
lnc_rna[,1:(dim(lnc_rna)[2]-5)] <- log1p(lnc_rna[,1:(dim(lnc_rna)[2]-5)])

#2. Get lncRNA - median within each tissue type
tissues <- unique(lnc_rna$canc)

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- lnc_rna[lnc_rna$canc==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA-RNAseq data 
#output: lncrna-Median ranked by percentile 

get_gene_medians <- function(dataframe){
	tis <- dataframe$canc[1]
	meds <- apply(dataframe[,1:(dim(dataframe)[2]-5)], 2, median)
	bg_by_bin = by(order(meds), ceiling(seq_along(1:length(meds))/(length(meds)/100)), function(x) meds[x])
	#need to get canc-gene-percentile 
	data2 <- data.frame(Gene="",Percent = "", Tissue="")     
	for(i in 1:length(bg_by_bin)){
		percentile <- i
		genes <- names(bg_by_bin[[i]])
		data <- data.frame(Gene=genes,Percent = percentile, Tissue=tis)     
		data2 <- rbind(data2, data)
	}
	data2 <- data2[-1,]
	return(data2)
}

binned_lncs <- lapply(tissues_data, get_gene_medians)

#Function 3
#input: dataframe with lncrna-percentile-tissue
#output: print dataframe so can compare with PCAWG

print_ranked <- function(dataframe){
	filename <- paste(dataframe$Tissue[1], "ranked_by_bins_lncRNAsPCAWG.txt", sep="_")
	write.table(dataframe, file=filename, quote=F, row.names=F, sep=";")
}	
lapply(binned_lncs, print_ranked)































