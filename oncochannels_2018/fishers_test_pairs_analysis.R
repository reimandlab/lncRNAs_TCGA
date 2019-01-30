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
  ic_cancs[[i]]$canc = names(expr_discr_cc)[[i]]
}


unique_cds = unique(names(mutations_in_cds)) #48
unique_fmre = unique(names(mutations_in_crms)) #30

results_pairs = as.data.frame(matrix(ncol=6)) ; colnames(results_pairs) = c("CDS_mut", "FMRE_mut", "fishers_pval", "fishers_OR", "num_overlap", "canc_overlapping_pats")

#for each cds/fmre combo
for(i in 1:length(unique_cds)){
  for(y in 1:length(unique_fmre)){
    pair = c(unique_cds[i], unique_fmre[y])
    #get patients that have either of these mutations or both 
    #FMRE
    z1 = which(names(mutations_in_crms) %in% pair)
    pats_crms = as.data.frame(mutations_in_crms[z1])
    pats_crms = as.data.frame(pats_crms[!duplicated(pats_crms), ])
    pats_crms$mut = "FMRE"
    colnames(pats_crms)[1] = "patient"
    #CDS
    z2 = which(names(mutations_in_cds) %in% pair)
    pats_cds = as.data.frame(mutations_in_cds[z2])
    pats_cds = as.data.frame(pats_cds[!duplicated(pats_cds), ])
    pats_cds$mut = "CDS"
    colnames(pats_cds)[1] = "patient"
    #combine 
    patients_wmuts = rbind(pats_crms, pats_cds)
    patients_wmuts = patients_wmuts[!duplicated(patients_wmuts), ]
    
    #set up contigency table
    #number of pats with both muts, #fmre only, #cds only, #no muts
    #both 
    both = as.data.table(table(patients_wmuts$patient))
    both = filter(both, N ==2)
    both_pats = both$V1
    both = dim(both)[1]
    #overlap
    num_overlap = both
    
    #test only if have at least 1 patient in common 
    if(!(length(both_pats)==0)){
      canc_both_pats = paste(unique(patient_table$V2[patient_table$V1 %in% both_pats]), collapse="_")
      
      #fmre only
      fmre = unique(as.character(patients_wmuts$patient[patients_wmuts$mut == "FMRE"]))
      fmre = fmre[-(which(fmre %in% both_pats))]
      FMRE_yes = length(unique(fmre))
      
      #cds only
      cds = unique(patients_wmuts$patient[patients_wmuts$mut == "CDS"])
      cds = cds[-which(cds %in% both_pats)]
      CDS_yes = length(unique(cds))
      
      #none 
      none = length(which(!(patient_table$V1 %in% unique(patients_wmuts$patient))))
      
      #contigency table 
      cont_table = matrix(c(both, CDS_yes,FMRE_yes , none), nrow=2, byrow=T)
      colnames(cont_table) = c("FMRE_yes", "FMRE_no")
      rownames(cont_table) = c("CDS_yes", "CDS__no")
      #p <-tableGrob(cont_table)
      #print(grid.arrange(top=paste(pair[1], pair[2]), p))
      f = fisher.test(cont_table, alt = "greater")
      row = c(pair, f$p.value, f$estimate, num_overlap, canc_both_pats)
      names(row) = colnames(results_pairs)
      results_pairs = rbind(results_pairs, row)
    }
  }
}





