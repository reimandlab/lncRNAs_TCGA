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
#z <- which(duplicated(ucsc[,8]))
#ucsc <- ucsc[-z,]

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
#allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_Aug8.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands = subset(val_cands, as.numeric(pval) < 0.05)

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

lncs = (unique(allCands$gene))

#GET COORDINATES AND SAVE AS BED FILES THEN USE BEDTOOLS 

#1. Get coordinates of all lncRNA candidates 
# & make into genomic ranges object 
  
  z = which(ucsc$hg19.ensGene.name2 %in% lncs)
  lncs_cords = ucsc[z,]
  lncs_cords = lncs_cords[,c("hg19.ensGene.chrom", "hg19.ensGene.txStart", 
  "hg19.ensGene.txEnd", "hg19.ensGene.strand", "hg19.ensemblToGeneName.value",
  "hg19.ensemblSource.source")]
  colnames(lncs_cords) = c("chr", "start", "end", "strand", "name", "type")
  #missing lncRNA 
  lncs_cords = as.data.table(lncs_cords)
  lncs_cords = lncs_cords[,c("chr", "start", "end", "name", "strand", "type")]
  write.table(lncs_cords, file="lncs_cords_166_cands_aug8.bed", quote=F, row.names=F, col.names=F, sep="\t")


#2. Get coordinates of all PCGs 
# & make into genomic ranges objects 

  z = which(ucsc$hg19.ensGene.name2 %in% colnames(pcg))
  pcgs = ucsc[z,]
  pcgs_cords = pcgs[,c("hg19.ensGene.chrom", "hg19.ensGene.txStart", 
  "hg19.ensGene.txEnd", "hg19.ensGene.strand", "hg19.ensemblToGeneName.value",
  "hg19.ensemblSource.source")]
  pcgs_cords = as.data.table(pcgs_cords)
  colnames(pcgs_cords) = c("chr", "start", "end", "strand", "name", "type")
  z = which(str_detect(pcgs_cords$chr, "_"))
  pcgs_cords = pcgs_cords[-z]
  pcgs_cords = pcgs_cords[,c("chr", "start", "end", "name", "strand", "type")]
  write.table(pcgs_cords, file="pcgs_cords_allPCGs_aug8.bed", quote=F, row.names=F, col.names=F, sep="\t")


##--------------------------------------------------------------------
###in bedtools 
##--------------------------------------------------------------------

module load bedtools 

#1. first need to sort lncRNA coordinates 
sort -k1,1 -k2,2n lncs_cords_166_cands_aug8.bed > lncs_cords_166_cands_aug8.sorted.bed
sort -k1,1 -k2,2n pcgs_cords_allPCGs_aug8.bed > pcgs_cords_allPCGs_aug8.sorted.bed

#2. merge transcripts into one 
bedtools merge -i lncs_cords_166_cands_aug8.sorted.bed -c 4 -o collapse -delim "|" > lncs_cords_166_cands_aug8.sorted.merged.bed
bedtools merge -i pcgs_cords_allPCGs_aug8.sorted.bed -c 4 -o collapse -delim "|" > pcgs_cords_allPCGs_aug8.sorted.merged.bed

#3. get file of lncRNA-cis interactions 
bedtools window -a lncs_cords_166_cands_aug8.sorted.bed -b pcgs_cords_allPCGs_aug8.sorted.bed -w 5000 > lncs_candidates_pcgs_intersected.bed

#get file of lncRNA-trans interactions 

##--------------------------------------------------------------------
###in R 
##--------------------------------------------------------------------

cis_ints = read.table("lncs_candidates_pcgs_intersected.bed")
colnames(cis_ints) = c("lnc_chr", "lnc_start", "lnc_end", "lnc", "lnc_strand", "lnc_type", "pcg_chr", 
  "pcg_start", "pcg_end", "pcg", "pcg_strand", "protein_type")
cis_ints$pair = paste(cis_ints$lnc, cis_ints$pcg, sep="_")
#keep only unique pairs 
cis_ints = as.data.table(cis_ints)
z = which(duplicated(cis_ints$pair))
cis_ints = cis_ints[-z,]

#80/166 lncRNAs have at least 1 nearby PCG
saveRDS(cis_ints, file="lncRNA_cands_wPCGs_that_are_in_cis_aug8.rds")

#all the lncRNA candidates that are not in the above ^ file are "trans lncRNAs"












