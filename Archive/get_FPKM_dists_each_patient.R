#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#List of canddidates and cox results
#allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
#List of canddidates and cox results
#allCands <- readRDS("all_candidates_combined_cancers_typesAnalysis_May3rd.rds")
#cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")

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

tcga_genes = fread("all_genes_used_inRankingAnalysisTCGA_May4th.txt")
tcga_genes$type = ""
tcga_genes$type[tcga_genes$x %in% fantom$CAT_geneID] = "lncRNA"
tcga_genes$type[is.na(tcga_genes$type)] = "pcg"

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

z = which(str_detect(colnames(rna), "ENSG"))	

#2. Get lncRNA - median within each tissue type
tissues <- unique(rna$Cancer)
#tissues <- tissues[c(7,9,12,13)]

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- rna[rna$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")


#Function 2 get distribution of FPKM for each lncRNA for each patient 
#colour by patient? 

getScores <- function(patient){
	z = which(names %in% patient) 
	row = dataframe[z,]
	z = which(str_detect(names(row), "ENSG"))

	expression <- data.frame(exp=as.numeric(row[z]), gene=names(row)[z])
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	
	#subset to just lncrnas
	lncs = tcga_genes$x[tcga_genes$type == "lncRNA"]
	z <- which(expression$gene %in% lncs)
	expression <- expression[z, ]
	expression = as.data.frame(expression)
	expression$patient = ""
	expression$patient = patient 
	return(expression)
}

addScores <- function(dataframe){
	names <- dataframe$patient
	patients <- llply(names, getScores, .progress="text") #list of dataframes, need to coerce together
	#names <- rownames(dataframe)
	patients1 <- rbindlist(patients)
	patients1 <- as.data.frame(patients1)
	patients1$canc <- dataframe$Cancer[1]
	patients1 = subset(patients1, exp < summary(patients1$exp)[5])
	patients1$exp = log1p(patients1$exp)
	#plot distribution of FPKM 
	ggdensity(patients1, x="exp", color="patient", legend="none", rug = TRUE, title=paste("log1p FPKM by patient", dataframe$Cancer[1]))
	return(patients1)
}	

tissues_data = tissues_data[1:3]
scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_tissues_scored <-  rbindlist(scored)
dev.off()








