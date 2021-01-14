library(corrplot)

set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
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

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

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


#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
allCands$combo=paste(allCands$gene, allCands$cancer, sep="_")
allCands$gene_name = allCands$gene_symbol

#--------This script ------------------------------------------------

#intersect lncRNA coordinates with PCGs to obtain dataset of
#lncRNA - cis PCGs and lncRNA - trans PCGs

#--------------------------------------------------------------------
#ANALYSIS------------------------------------------------------------
#--------------------------------------------------------------------

library(plyr)
library(dplyr)

lncs = (unique(allCands$gene))

#GET COORDINATES AND SAVE AS BED FILES THEN USE BEDTOOLS

#1. Get coordinates of all lncRNA candidates
# & make into genomic ranges object

  z = which(hg38$ensgene %in% lncs)
  lncs_cords = hg38[z,]
  lncs_cords = lncs_cords[,c("chr", "start",
  "end", "strand", "symbol",
  "biotype")]
  colnames(lncs_cords) = c("chr", "start", "end", "strand", "name", "type")
  #missing lncRNA
  lncs_cords = as.data.table(lncs_cords)
  lncs_cords = lncs_cords[,c("chr", "start", "end", "name", "strand", "type")]
  write.table(lncs_cords, file="lncs_cords_142_cands.bed", quote=F, row.names=F, col.names=F, sep="\t")


#2. Get coordinates of all PCGs
# & make into genomic ranges objects

  z = which(hg38$ensgene %in% colnames(pcg))
  pcgs_cords = hg38[z,]
  pcgs_cords = pcgs_cords[,c("chr", "start",
  "end", "strand", "symbol",
  "biotype")]
  pcgs_cords = as.data.table(pcgs_cords)
  colnames(pcgs_cords) = c("chr", "start", "end", "strand", "name", "type")
  #z = which(str_detect(pcgs_cords$chr, "_"))
  #pcgs_cords = pcgs_cords[-z]
  pcgs_cords = pcgs_cords[,c("chr", "start", "end", "name", "strand", "type")]
  write.table(pcgs_cords, file="pcgs_cords_142_cands.bed", quote=F, row.names=F, col.names=F, sep="\t")


lncs_cords$strand[lncs_cords$strand == "1"] = "+"
lncs_cords$strand[lncs_cords$strand == "-1"] = "-"
pcgs_cords$strand[pcgs_cords$strand == "1"] = "+"
pcgs_cords$strand[pcgs_cords$strand == "-1"] = "-"

lncs_cords_gr = makeGRangesFromDataFrame(lncs_cords)
pcgs_cords_gr = makeGRangesFromDataFrame(pcgs_cords)

hits <- findOverlaps(lncs_cords_gr, pcgs_cords_gr, ignore.strand=TRUE, maxgap=10000)
hits_overlap = cbind(as.data.table(lncs_cords[queryHits(hits),]), as.data.table(pcgs_cords)[subjectHits(hits),])
print(head(hits_overlap))

colnames(hits_overlap) = c("lnc_chr", "lnc_start", "lnc_end", "lnc", "lnc_strand", "lnc_type", "pcg_chr",
  "pcg_start", "pcg_end", "pcg", "pcg_strand", "protein_type")
hits_overlap$pair = paste(hits_overlap$lnc, hits_overlap$pcg, sep="_")
#keep only unique pairs
hits_overlap = as.data.table(hits_overlap)
hits_overlap$distance = abs(as.numeric(hits_overlap$pcg_start) - as.numeric(hits_overlap$lnc_end))

z = which(duplicated(hits_overlap$pair))
hits_overlap = hits_overlap[-z,]

#95/179 lncRNAs have at least 1 nearby PCG (152 unique PCGs)
saveRDS(hits_overlap, file="lncRNA_cands_wPCGs_that_are_in_cis_10kb_nov16.rds")

#all the lncRNA candidates that are not in the above ^ file are "trans lncRNAs"
