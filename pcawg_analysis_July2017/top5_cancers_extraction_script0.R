#top5_cancers_extraction_script0.R

#Karina Isaev
#July 19th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script0 - extracts cancer types and plot distributions 

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library("colorout")
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

#---------------------------------------------------------
#Processing
#---------------------------------------------------------

# [1] Match RNA-Seq column names to conversion file in 
# order to get icgc donor IDs

#conversion$normal_rna_seq_submitter_sample_id = 63
#length(which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))) = 161
#length(which(conversion$tumor_rna_seq_submitter_sample_id %in% colnames(rna))) = 819
#length(which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))) = 1267

submitter <- conversion$tumor_rna_seq_submitter_sample_id[which(conversion$tumor_rna_seq_submitter_sample_id %in% colnames(rna))]
aliqot <- conversion$tumor_rna_seq_aliquot_id[which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))]

length(which(aliqot %in% submitter))
#[1] 819

#so aliqot has all the same IDs from submitter and so can just use aliqot ids 

###"NORMAL SAMPLES"
z <- which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))
norm_pats <- conversion$icgc_donor_id[z]

###"TUMOUR SAMPLES"
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z]

###161 tum-matched normal samples -- but what cancer types are they from?
tum_match_norm <- tum_pats[which(tum_pats %in% norm_pats)]

z <- which(clin$icgc_donor_id %in% tum_match_norm) #324
tum_match_norm_clin <- clin[z,]
z <- which(duplicated(tum_match_norm_clin$icgc_donor_id))
tum_match_norm_clin <- tum_match_norm_clin[-z,]
#all 161 unique IDs have clinical information

#liver and kidney cancers have >50 tumour matched normal samples

###How many of total 1267 samples have clinical data?
z <- which(clin$icgc_donor_id %in% tum_pats) #2532
tum_clin <- clin[z,]
z <- which(duplicated(tum_clin$icgc_donor_id))
tum_clin <- tum_clin[-z,]
#all 1267 tumour RNA-seq samples also have matching clinical information 

#what is the distribution of samples/cancer type

###TUMOUR-MATCHED-NORMAL
tum_types_match <- as.data.table(table(tum_match_norm_clin$histology_tier2))
tum_types_match <- tum_types_match[order(N)]

pdf("161_tum_matched_sampels_cancer_distribution.pdf")
p <- ggbarplot(tum_types_match, "V1", "N", color="N", xlab="Tumour Histology") +
ggtitle("Distribution of Matched Tumour Normal RNA-Seq Samples")
ggpar(p,
 font.tickslab = c(8,"bold", "black"),
 xtickslab.rt = 65, legend = "right")
dev.off()

#---------------------------------------------------------
#So what are the top 5 cancers (most samples with clinical)
#and RNA-Seq data?
#---------------------------------------------------------

###JUST CANCER SITE
tum_types <- as.data.table(table(tum_clin$histology_tier2))
tum_types <- tum_types[order(N)]
tum_types <- tum_types[-1,]

pdf("1267_cancer_mainSITE_organ_distribution.pdf")
p <- ggbarplot(tum_types, "V1", "N", color="N", xlab="Tumour Organ Site") +
ggtitle("Distribution of Tumour RNA-Seq Samples by Cancer Site")
ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="right")
dev.off()

###CANCER TYPES WITH AT LEAST 50 PATIENTS 
tum_types <- tum_types[tum_types$N>=50]
pdf("50patients_min_cancers_dist.pdf")
p <- ggbarplot(tum_types, "V1", "N", color="N", xlab="Tumour Organ Site") +
ggtitle("Distribution of Tumour RNA-Seq Samples by Cancer Site")
ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="right")
dev.off()

###IN THESE CANCER TYPES, HOW MANY HISTOLOGICAL 
#SUBTYPES ARE THERE?
canc <- tum_types[,1] ; canc <- as.vector(unlist(canc))
z <- which(tum_clin$histology_tier2 %in% canc)
top50_tum_clin <- tum_clin[z,]

tum_types <- as.data.table(table(top50_tum_clin$histology_tier2, top50_tum_clin$histology_tier4))
tum_types <- tum_types[order(N)]
z <- which(tum_types$N ==0)
tum_types <- tum_types[-z,]

cols <- colorRampPalette(mypal)(30) 

pdf("50patients_min_cancers_dist_wHISTOsubtype.pdf")
p <- ggbarplot(tum_types, "V1", "N", fill="V2", xlab="Tumour Organ Site", palette=cols) +
ggtitle("Distribution of Tumour RNA-Seq Samples by Histological Subtype")
ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="right", font.legend = c(5, "plain", "black"), legend.title="Histological Subtype")
dev.off()


###MOST PATIENTS + GOOD DISTRIBUTION AMONG SUBTYPES
#--> KIDNEY, OVARY, LIVER, PANCREAS, BREAST 

canc <- c("Kidney", "Ovary", "Liver", "Breast", "Pancreas") ; canc <- as.vector(unlist(canc))
z <- which(tum_clin$histology_tier2 %in% canc)
top50_tum_clin <- tum_clin[z,]

tum_types <- as.data.table(table(top50_tum_clin$histology_tier2, top50_tum_clin$histology_tier4))
tum_types <- tum_types[order(N)]
z <- which(tum_types$N ==0)
tum_types <- tum_types[-z,]

cols <- colorRampPalette(mypal)(17) 

pdf("50patients_min_cancers_dist_wHISTOsubtype_TOP5.pdf")
p <- ggbarplot(tum_types, "V1", "N", fill="V2", xlab="Tumour Organ Site", palette=cols) +
ggtitle("Distribution of Tumour RNA-Seq Samples by Histological Subtype")
ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="right", font.legend = c(5, "plain", "black"), legend.title="Histological Subtype")
dev.off()

###WITHIN CANCER TYPE, REMOVE HISTO SUBTYPE IF ONLY 1-10 SAMPLES IN IT 
tum_types <- tum_types[tum_types$N>=10]

pdf("50patients_min_cancers_dist_wHISTOsubtype_TOP5_removedLOwhistoTypes.pdf")
p <- ggbarplot(tum_types, "V1", "N", fill="V2", xlab="Tumour Organ Site", palette=mypal) +
ggtitle("Distribution of Tumour RNA-Seq Samples by Histological Subtype")
ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="right", font.legend = c(5, "plain", "black"), legend.title="Histological Subtype")
dev.off()




















