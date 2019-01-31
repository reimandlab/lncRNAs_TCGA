library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#this script below prepares the RNA and clinical files for analysis 
source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

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

#Combined into one dataframe because need to get ranks 
com = colnames(pcg)[which(colnames(pcg) %in% colnames(rna))]
all <- merge(rna, pcg, by = com)

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
cands = merge(cands, ucsc, by = "HGNC.symbol")

#outlier gene exp 
load("expr_discr_cc.rsav")

change_cols = function(dat){
  cols = sapply(colnames(dat), function(x){gsub("\\.", "-", x)})
  colnames(dat) = cols
  return(dat)
}

ic_cancs = llply(expr_discr_cc, change_cols)

for(i in 1:length(ic_cancs)){
  dff = ic_cancs[[i]]
  dff = t(dff)
  dff = as.data.frame(dff)
  dff$canc = ""
  dff$canc = names(expr_discr_cc)[[i]]
  ic_cancs[[i]] = dff
}


###MAIN FUNCTION#####################################################

get_fishers = function(canc_dat){
  
  cancer = canc_dat$canc[1]
  canc_dat$canc = NULL
  
  #keep only ion channels 
  z = which(colnames(canc_dat) %in% cands$HGNC.symbol)
  canc_dat = canc_dat[,z]

  ics1 = colnames(canc_dat)

  pairs = combn(ics1, 2)
  pairs = t(pairs)
  pairs = as.data.table(pairs)
  colnames(pairs) = c("IC1", "IC2")

  results_pairs = as.data.frame(matrix(ncol=8)) ; colnames(results_pairs) = c("IC1", "IC2", "fishers_pval", "fishers_OR", "num_overlap", "IC1_high", "IC2_high", "cancer")

  #for each cds/fmre combo
  for(i in 1:nrow(pairs)){
    print(i)

    pair = pairs[i,]
    
    #get patients that have either of these mutations or both 
    z = which(colnames(canc_dat) %in% pair)
    pair_dat = canc_dat[,z]  
    cont_table = table(pair_dat)
    if((dim(cont_table)[2] ==2) & (dim(cont_table)[1] ==2)){

    f = fisher.test(cont_table, alt = "greater")
    row = c(pair[[1]], pair[[2]], f$p.value, f$estimate, cont_table[4], cont_table[2], cont_table[3], cancer)
    names(row) = colnames(results_pairs)
    results_pairs = rbind(results_pairs, row)

        }
      }

  results_pairs = as.data.table(results_pairs)
  results_pairs = results_pairs[-1,]
  results_pairs$fishers_pval = as.numeric(results_pairs$fishers_pval)
  results_pairs$fdr = p.adjust(results_pairs$fishers_pval, method="fdr")
  results_pairs = results_pairs[order(fdr)]
  return(results_pairs)
}


all_fishers = llply(ic_cancs, get_fishers, .progress="text")

































