###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
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

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
rna = rna[-which(rna$Cancer %in% canc_rm),]

#------FUNCTIONS----------------------------------------------------

#######
##[1]##-------------------------------------------------------------
#######

#get gene ID from name 
#input: GeneName
#output: GeneID 

get_ensg = function(lnc){
	z = which(fantom$CAT_geneName == lnc)
	return(fantom$CAT_geneID[z])
}

#######
##[2]##-------------------------------------------------------------
#######

#get lncRNA expression across cancers plot 
#input: lncRNA id
#output: boxplot of expression across cancer types

#test case: pca3 - ENSG00000225937

get_exp_plots = function(lnc){
	dat = rna[,which(colnames(rna) %in% c("type", lnc))]
	colnames(dat)[1] = "lncRNAExp"
	dat$lncRNAExp = log1p(dat$lncRN)
	dat = as.data.table(dat)
	dat = dat[order(lncRNAExp)]

	#get order of cancer types by median 
	sum = as.data.table(dat %>% dplyr::group_by(type) %>% 
		dplyr::summarize(median=median(lncRNAExp)))
	sum = sum[order(median)]

	dat$type = factor(dat$type, levels=unique(sum$type))
	#p = ggboxplot(dat, x="type", y = "lncRNAExp", title=lnc, error.plot = "linerange")+
	#rotate_x_text(45)
	
	#try density
	p = ggdensity(dat, x = "lncRNAExp", title=lnc, color="type")
	print(p)
	print("done plot")
}

#######
##[3]##-------------------------------------------------------------
#######

#check how many patients have expression value X FPKM 
#(because we want at least 15 patients with 100FPKM expression)

get_num_pats = function(lnc, canc, exp_cut){
	dat = rna[,which(colnames(rna) %in% c("Cancer", lnc))]
	dat = dat[which(dat$Cancer == canc),]
	colnames(dat)[1] = "lncRNAExp"
	num_pats = length(which(dat[,1] > exp_cut))
	if(num_pats >= 15){
		return("great success")
	}
	if(num_pats < 15){
		return("problem")
	}
}

















