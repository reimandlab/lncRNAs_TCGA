###---------------------------------------------------------------
###TCGA_cancers_survival1.R
###---------------------------------------------------------------

###October 11th, 2017
###Goal: Using list of candidate lncRNAs obtained in PCAWG analysis
#validate their survival association in this new TCGA data

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. RNA data 

#--rna (includes ALL rna)
pcg = readRDS("60483_rnas_all_tumours_9602_samples_TCGAnew.rds")
rownames(pcg) = pcg$gene
pcg$gene = NULL
pcg = t(pcg)

#--normal 
norm = readRDS("all_genes_563_matched_normal_samples_TCGA_April11.rds") #numbers are actually different for these two files
rownames(norm) = norm$gene
norm$gene = NULL
norm = t(norm)

#--metastatic
met = readRDS("all_genes_354_matched_metastatic_tumours_TCGA_april.rds") #numbers are actually different for these two files
rownames(met) = met$gene
met$gene = NULL
met = t(met)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

#let's look only at protein-coding genes (ion channels).... lncRNAs and miRNAs can also be selected for if needed later...
z <- which(ucsc$hg19.ensemblSource.source %in% c("protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#keep only genes in pcg that are labelled as PCG
z <- which(colnames(pcg) %in% ucsc$hg19.ensGene.name2)
pcg = pcg[,z]

#2. Clinical data
#clin = read.csv("all_clin_XML_tcgaSept2017.csv")
#clin = clin[,1:90]

clin = fread("mmc1_clinical_data_cellpaper2018.txt")

#3. TCGA ID cancer type conversion 
canc_conversion = readRDS("tcga_id_cancer_type_conversion_oncochannels.txt")
norm_conversion = readRDS("tcga_id_NORMAL_samples_type_conversion_oncochannels.txt")
met_conversion = readRDS("tcga_id_Metastatic_samples_type_conversion_oncochannels.txt")

#norm
z = which(colnames(norm) %in% colnames(pcg))
norm = norm[,z]

#met
z = which(colnames(met) %in% colnames(pcg))
met = met[,z]

###---------------------------------------------------------------
###Process Data 
###---------------------------------------------------------------

#Change patient ids to shorted id
change = function(rowname){
  new = canc_conversion$id[which(canc_conversion$TCGA_id %in% rowname)]
  return(new)  
}

rownames(pcg) = sapply(rownames(pcg), change)

change = function(rowname){
  new = norm_conversion$id[which(norm_conversion$TCGA_id %in% rowname)]
  return(new)  
}
rownames(norm) = sapply(rownames(norm), change)

change = function(rowname){
  new = met_conversion$id[which(met_conversion$TCGA_id %in% rowname)]
  return(new)  
}
rownames(met) = sapply(rownames(met), change)

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(pcg) %in% clin$bcr_patient_barcode)
pcg = pcg[z,] #all have clinical data - 9602 patients 
z <- which(rownames(norm) %in% clin$bcr_patient_barcode)
norm = norm[z,] #all have clinical data - 699 patients (matched normals)
z <- which(rownames(met) %in% clin$bcr_patient_barcode)
met = met[z,] #all have clinical data - 393 patients (matched normals)

#Add survival info to RNA file - normal matched 
norm = as.data.frame(norm)
norm$patient = rownames(norm)
colnames(clin)[2] = "patient"
norm = merge(norm, clin, by="patient")

saveRDS(norm, "all_genes_RNASeq_matched_normals_wclinical_data_oncochannels.rds")

#Add survival info to RNA file - metastatic matched 
met = as.data.frame(met)
met$patient = rownames(met)
colnames(clin)[2] = "patient"
met = merge(met, clin, by="patient")

saveRDS(met, "all_genes_RNASeq_matched_metastatic_wclinical_data.rds")

#Add survival info to PCG file
pcg = as.data.frame(pcg)
pcg$patient = rownames(pcg)
pcg = merge(pcg, clin, by="patient")

saveRDS(pcg, "19438_PCGS_RNASeq_tcga_all_cancers_wclinical_data.rds")





/.mounts/labs/reimandlab/private/projects/oncochannels/data



















































