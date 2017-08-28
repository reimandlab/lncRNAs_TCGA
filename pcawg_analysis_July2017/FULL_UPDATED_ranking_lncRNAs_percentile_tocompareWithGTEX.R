#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 28th, 2017

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

#List of canddidates and cox results
allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
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

#make sure there aren't any lncRNAs in protein coding file 
table(ucsc[,7][ucsc[,8] %in% colnames(pcg_rna)])
#Remove 
z <- which(colnames(pcg_rna) %in% fantom[,2])
pcg_rna <- pcg_rna[,-z]

#lncs
lncs <- lnc_rna[,-c(5608,5609)]

#Combined into one dataframe because need to get ranks 
all <- cbind(lncs, pcg_rna)

#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------

#For each patient add survival status and days since last seen 
all$status <- ""
all$time <- ""
all$sex <- ""

#lncs
for(i in 1:nrow(all)){
  pat <- rownames(all)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  all$status[i] <- clin$donor_vital_status[z]
  all$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        all$time[i] <- t
}

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
all[,1:(dim(all)[2]-5)] <- log1p(all[,1:(dim(all)[2]-5)])

#2. Get lncRNA - median within each tissue type
tissues <- unique(all$canc)
tissues <- tissues[c(1,2,4,5,6)]

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- all[all$canc==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 
getScores <- function(row){
	score=""
	expression <- data.frame(exp=as.numeric(row[1:(length(row)-5)]), gene=names(row)[1:(length(row)-5)])
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lnc candidates - we just want their score 
	z <- which(expression$gene %in% allCands$gene)
	expression <- expression[z, ]
	return(expression)
}


addScores <- function(dataframe){
	patients <- apply(dataframe, 1, getScores) #list of dataframes, need to coerce together
	names <- rownames(dataframe)
	patients <- rbindlist(patients)
	patients$patient <- rep(names, each=38) #38 lncRNA candidates 
	patients <- as.data.frame(patients)
	patients$canc <- dataframe$canc[1]
	patients$data <- "PCAWG"
	patients$canc <- lapply(patients$canc, function(x) unlist(strsplit(x, " "))[1])
	return(patients)
}	

scored <- lapply(tissues_data, addScores) #list of dataframes
all_cancers_scored <-  rbindlist(scored)
all_cancers_scored <- as.data.frame(all_cancers_scored)

#one for each tissue/cancer type
#each gene is scored within each patient 
#can now make violin plot showing distributon of scores for each candidate lncRNA 
#just need to subset to genes interested in plotting 

#write file so can use with GTEX 
saveRDS(all_cancers_scored, file="all_cancers_scored.rds")

#save list of genes in total used to also compare with GTEX 
write.table(colnames(all), file="all_genes_used_inRankingAnalysis.txt", quote=F, row.names=F, sep=";")





























