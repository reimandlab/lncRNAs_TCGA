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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#--------This script ------------------------------------------------

#summarize results from co-expression analysis of PCGs
#how many per lcnRNA
#how many pathways per lncRNA
#how many cancer genes, RBPs, TFs...

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#--------------------------------------------------------------------
#RESULTS-------------------------------------------------------------
#--------------------------------------------------------------------

coexp = readRDS("all_results_for_each_cancer_from_coexpression_analysis_june27th_allCands.rds")

#can run FDR on all PCGs for all Cancer types 
#175 * 20,0000 = enormous FDR test 

#1. Get FDR by cancer type 
cancers = as.list(unique(coexp$canc))
canc_fdr = function(cancer){
  dat = subset(coexp, canc == cancer)
	dat$pvalue = as.numeric(dat$pvalue)
  dat$fdr = p.adjust(dat$pvalue, method = "fdr")
	z = which(is.na(dat$pvalue))
  if(!length(z)==0){
    dat = dat[-z,]
  }
  dat = as.data.table(dat)
	dat = dat[order(fdr)]
	return(dat)
}

library(dplyr)
library(plyr)
canc_dats = llply(cancers, canc_fdr) 
canc_dats = ldply(canc_dats, data.frame)

#plot scatter plot, FC versus p-value
#coexp$fdr = -log10(coexp$fdr)
canc_dats$mean_diff = as.numeric(canc_dats$mean_diff)
#z = which(canc_dats$mean_diff == 0)
#canc_dats = canc_dats[-z,]
canc_dats = as.data.table(canc_dats)
#coexp = filter(coexp, pvalue <= 0.05)

#ggscatter(coexp, x = "mean_diff", y = "fdr", size=0.5, 
#   color="fdr") + geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 0) +
#    geom_vline(xintercept = (log1p(4))-(log1p(2))) +
#    geom_vline(xintercept = (log1p(2))-(log1p(4)))
#dev.off()


#2. Summarize per lncRNA/cancer, how many PCGs upregulated in risk group
#and how many upregulated in non-risk group 
#(>0 --> more expressed in risk group, <0, more expressed in low risk group)

#keep sig fdr, add risk/non-risk tag
canc_dats = as.data.table(canc_dats)
canc_dats = as.data.table(filter(canc_dats, fdr <=0.01))
canc_dats$risk_type = ""
canc_dats$risk_type[canc_dats$mean_diff < (log1p(2)-log1p(4))] = "NonRisk"
canc_dats$risk_type[canc_dats$mean_diff >= (log1p(4)-log1p(2))] = "Risk"
canc_dats = canc_dats[which(!(canc_dats$risk_type =="")),]
canc_dats$combo = paste(canc_dats$lnc, canc_dats$canc, sep="_")

#PCG lncRNA results
pcg_lnc = readRDS("summary_pcg_analysis_wHRs_july17.rds")
pcg_lnc = pcg_lnc[order(-NumPCGs)]
pcg_lnc = as.data.table(filter(pcg_lnc, NumPCGs >=500))

#-------------------ANALYSIS--------------------------------------------
#Generate heatmap using those PCGs sig up/down regulated in lncRNA 
#risk or non-risk groups

combos = unique(pcg_lnc$combo)

gen_heatmap = function(lnc_canc_combo){
  #get which pcgs enriched in group
  z = which(canc_dats$combo %in% lnc_canc_combo)
  #cancer type
  #subset gene expression to those pcgs
  #label patients by either high/low lncRNA expression 

}














