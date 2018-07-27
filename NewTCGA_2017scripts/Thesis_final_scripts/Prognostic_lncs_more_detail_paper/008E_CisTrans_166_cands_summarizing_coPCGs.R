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
library(GenomicRanges)

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
#all <- merge(rna, pcg, by = c("patient", "Cancer"))
#all = all[,1:25170]

#------FUNCTIONS-----------------------------------------------------


#--------This script ------------------------------------------------

#intersect lncRNA coordinates with PCGs to obtain dataset of 
#lncRNA - cis PCGs and lncRNA - trans PCGs 

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#--------------------------------------------------------------------
#ANALYSIS------------------------------------------------------------
#--------------------------------------------------------------------

library(plyr)
library(dplyr)

#1. Get coordinates of all lncRNA candidates 
# & make into genomic ranges object 

#2. Get coordinates of all PCGs 
# & make into genomic ranges objects 

#convert to granges 
pcgs = makeGRangesFromDataFrame(pcgs)

countOverlaps(lncs_cords, pcgs)
findOverlaps(lncs_cords, pcgs)
subsetByOverlaps(lncs_cords, pcgs)

lncs = as.list(unique(allCands$gene))

get_cis = function(lnc){

  #lnc ----------------------------------------------------------------------
  z = which(ucsc$hg19.ensGene.name2 == lnc)
  lncs_cords = ucsc[z,]
  lncs_cords = lncs_cords[,c("hg19.ensGene.chrom", "hg19.ensGene.txStart", 
  "hg19.ensGene.txEnd", "hg19.ensGene.strand", "hg19.ensemblToGeneName.value",
  "hg19.ensemblSource.source")]
  lncs_cords = as.data.table(lncs_cords)
  colnames(lncs_cords) = c("chr", "start", "end", "strand", "name", "type")
  
  #convert to granges 
  lncs_cords_gr = makeGRangesFromDataFrame(lncs_cords)
  values(lncs_cords_gr) <- DataFrame(name = lncs_cords$name)


  #pcsg ---------------------------------------------------------------------
  z = which(ucsc$hg19.ensemblSource.source == "protein_coding")
  pcgs = ucsc[z,]
  pcgs = pcgs[,c("hg19.ensGene.chrom", "hg19.ensGene.txStart", 
  "hg19.ensGene.txEnd", "hg19.ensGene.strand", "hg19.ensemblToGeneName.value",
  "hg19.ensemblSource.source")]
  pcgs = as.data.table(pcgs)
  colnames(pcgs) = c("chr", "start", "end", "strand", "name", "type")
  z = which(str_detect(pcgs$chr, "_"))
  pcgs = pcgs[-z]

  #convert to granges 
  pcgs = makeGRangesFromDataFrame(pcgs)





}































