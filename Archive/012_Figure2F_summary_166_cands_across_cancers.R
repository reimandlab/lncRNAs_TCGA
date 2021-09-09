#------------------------------------------------------------------------------
#This script summarizes the 166 lncRNA candidates - which cancer types they
#belong to and what their HR is relative to multivariate model p-value (FDR) 
#------------------------------------------------------------------------------

#source code
source("check_lnc_exp_cancers.R")
library(corrplot)

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

#HiC data 
load("hic_data.rsav")
#remove rownames
rownames(hic_data) = c(1:nrow(hic_data))

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

#------PLOT----------------------------------------------------------

#geom_tile_plot -
#x = cancer type
#y = lncRNA
#tile fill = fdr multivariate model p-value
#tile shape = geom_point circle --> Hazard Ratio 
#colour of circle red versus blue 

#1. Get type of cancer for easier to read format

canc_conv = readRDS("cancers_conv_july23.rds")
