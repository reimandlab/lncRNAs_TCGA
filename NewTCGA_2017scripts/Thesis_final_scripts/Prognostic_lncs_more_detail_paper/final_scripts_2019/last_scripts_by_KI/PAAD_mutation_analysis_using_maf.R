setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
#require("powerSurvEpi")
library(SIBER)
library(EnvStats)

#------FEATURES-----------------------------------------------------

#cands -- should be this file
#cands = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

cands = readRDS("lncRNAs_selected_by_EN_april14.rds") #1000 runs of cross-validations using new updated dataset (GBM=124, OV and LUAD)

#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set
#fantom 

#UCSC gene info
ucsc <- fread("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#fantom 
fantom <- fread("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
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

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

filtered_data_tagged = readRDS("22_cancer_types_with_lncRNA_candidates_labelled_high_low.rds")
paad_tagged = filtered_data_tagged[[17]]
paad_maf = fread("PAAD_maf_hg38_2019.txt")

#subset to gene we are interested in = RNF43
paad_maf = as.data.table(filter(paad_maf, Hugo_Symbol == "RNF43"))

#get patient ids in file
paad_maf$patient = sapply(paad_maf$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
paad_tagged$patient = rownames(paad_tagged)
z = which(paad_tagged$patient %in% paad_maf$patient)
paad_tagged$RNF43_muts = ""
paad_tagged$RNF43_muts[z] = "yes"
paad_tagged$RNF43_muts[-z] = "no"

lnc = which(colnames(paad_tagged) == "ENSG00000265148")

table(paad_tagged[,lnc], paad_tagged$RNF43_muts)







