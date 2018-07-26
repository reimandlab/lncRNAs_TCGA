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

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#------FUNCTIONS-----------------------------------------------------

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

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

coexp = readRDS("coexpression_results_processed_july24.rds")

#PCG lncRNA results
pcg_lnc = readRDS("summary_pcg_analysis_wHRs_jul2y24.rds") #all these have at least 1, 50-pcg signature 
pcg_lnc = pcg_lnc[order(-NumPCGs)]
pcg_lnc$HR = as.numeric(pcg_lnc$HR)
pcg_lnc$lnc_stat = ""
pcg_lnc$lnc_stat[which(pcg_lnc$HR < 0)] = "Favourable"
pcg_lnc$lnc_stat[which(pcg_lnc$HR > 0)] = "Unfavourable"

#-------------------ANALYSIS--------------------------------------------
#Generate heatmap using those PCGs sig up/down regulated in lncRNA 
#risk or non-risk groups

#For each cancer type get all required data 
#PCG and lncRNA expression
combos = unique(pcg_lnc$combo)
cancs = sapply(combos, function(x){unlist(strsplit(x, "_"))[2]})
cancs = unique(cancs)

##1-----------------all expression--------------------------------------

get_tissue_specific <- function(combo){
  canc = unlist(strsplit(combo, "_"))[2]
  lnc = unlist(strsplit(combo, "_"))[1]
  tis = all[all$Cancer==canc,]
  tis$combo = combo
  print(combo)
  return(tis)
}
tissues_data <- llply(combos, get_tissue_specific, .progress="text")

##2-----------------label patients by risk------------------------------

get_lnc_canc = function(dat){
  cancer = dat$Cancer[1]
  combo = dat$combo[1]
  lnc = unlist(strsplit(combo, "_"))[1]

  pcgs = colnames(pcg)[2:19351]
  #keep only pcgs that are selected to be in lncRNA signature 
  z = which(coexp$combo == combo)
  lnc_pcgs = unique(coexp$pcg[z])

  dat_keep = dat[,which(colnames(dat) %in% c("patient", lnc, lnc_pcgs))]
  rownames(dat_keep) = dat_keep$patient
  dat_keep$patient = NULL
  #figure out which patients are high risk and which patients low risk
  dat_keep$median <- ""
  median2 <- quantile(as.numeric(dat_keep[,1]), 0.5)

       if(median2 ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat_keep[,1] > 0)
        l2 = which(dat_keep[,1] ==0)
        dat_keep$median[l1] = 1
        dat_keep$median[l2] = 0
        }

      if(!(median2 ==0)){
        l1 = which(dat_keep[,1] >= median2)
        l2 = which(dat_keep[,1] < median2)
        dat_keep$median[l1] = 1
        dat_keep$median[l2] = 0
    }

    #which one is high risk --> need surivval data
    dat_keep$patient = rownames(dat_keep)
    
      dat_keep$median[dat_keep$median ==0] = "Low"
      dat_keep$median[dat_keep$median==1] = "High"

      #cox ph
      z = which((allCands$gene == lnc) & (allCands$cancer == cancer))

      HR = as.numeric(allCands$HR[z])
      
      if(HR <1){
        risk = "Low"
        dat_keep$risk = ""
        dat_keep$risk[dat_keep$median=="High"] ="noRISK"
        dat_keep$risk[dat_keep$median=="Low"] ="RISK"
      }
      if(HR >1){
        risk = "High"
        dat_keep$risk = ""
        dat_keep$risk[dat_keep$median=="High"] ="RISK"
        dat_keep$risk[dat_keep$median=="Low"] ="noRISK"
      }

      dat_keep$lnc = colnames(dat_keep)[1]
      dat_keep$canc = cancer
      colnames(dat_keep)[1] = "lncRNA"

  return(dat_keep)  

}

#all lncRNAs with status 
all_canc_lnc_data = llply(tissues_data, get_lnc_canc, .progress="text")

##3-----------------get correlation pairs-----------------------------------

cancer = cancs[1]
library(tidyverse)

get_summary = function(cancer){
  
  #collect data from all lncRNAs in cancer type 
  keep = c()
  for(i in 1:length(all_canc_lnc_data)){
    z = all_canc_lnc_data[[i]]$canc[1] == cancer
    if(z){
      keep = c(keep, i)
    }
  }

  #cancer data
  canc_dats = all_canc_lnc_data[keep]
  canc_dats = reshape::merge_all(canc_dats)

  print(cancer)
    pcgs = colnames(canc_dats)[which(str_detect(colnames(canc_dats), "ENSG"))]
    lncs = unique(canc_dats$lnc)
    genes = c(pcgs, lncs)
    canc_exp = subset(all, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    
    res2 = rcorr(as.matrix(canc_exp), type="spearman")
    res2 = flattenCorrMatrix(res2$r, res2$P)
    res2$fdr = p.adjust(res2$p, method="fdr")
    res2 = as.data.table(res2)
    res2 = res2[order(fdr)]

    z = which(res2$row %in% lncs)
    res2$rowgene[z] = "lncRNA"
    res2$rowgene[-z] = "mRNA"
    z = which(res2$column %in% lncs)
    res2$columngene[z] = "lncRNA"
    res2$columngene[-z] = "mRNA"
    
    #total pairs 
    tot_pairs = nrow(res2)
    res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05, abs(cor) >= 0.3))
    sig_pairs = nrow(res2)    

    #%
    perc = sig_pairs/tot_pairs

    row = c(as.character(cancer), tot_pairs, sig_pairs, perc)
    return(row)

}

canc_results = llply(cancs, get_summary, .progress = "text")
#remove null
canc_results = Filter(Negate(is.null), canc_results)

canc_results = do.call(rbind.data.frame, canc_results)
colnames(canc_results) = c("cancer", "total_pairs", "sig_pairs", "perc")












