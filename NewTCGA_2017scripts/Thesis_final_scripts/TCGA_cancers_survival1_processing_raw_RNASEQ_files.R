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

#--lncrna
rna = readRDS("5919_all_tumours_7501_tissues_TCGAnew.rds")
rownames(rna) = rna$gene
rna$gene = NULL
rna = t(rna)

#--pcg
pcg = readRDS("54564_PCGs_all_tumours_7501_tissues_TCGAnew.rds")
rownames(pcg) = pcg$gene
pcg$gene = NULL
pcg = t(pcg)

#--normal 
norm = readRDS("all_genes_563_matched_normal_samples_TCGA_April11.rds")
rownames(norm) = norm$gene
norm$gene = NULL
norm = t(norm)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% "protein_coding")
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

#3. Fantom data 
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

#4. List of lncRNA survival associated candidates 
#cands = fread("7tier1_35tier2_lncRNA_candidates_August28th.txt")
#cands = fread("lncRNAs_sig_FDR_0.1_Nov23.txt")

#5. TCGA ID cancer type conversion 
canc_conversion = readRDS("tcga_id_cancer_type_conversion.txt")
norm_conversion = readRDS("tcga_id_NORMAL_samples_type_conversion.txt")

#6. List of TCGA IDs used in PCAWG - to remove
ids_remove = fread("TCGA_IDs_usedinPCAWG.txt")

#norm
z = which(colnames(norm) %in% c(colnames(rna), colnames(pcg)))
norm = norm[,z]

###---------------------------------------------------------------
###Process Data 
###---------------------------------------------------------------

#Change patient ids to shorted id
change = function(rowname){
  new = canc_conversion$id[which(canc_conversion$TCGA_id %in% rowname)]
  return(new)  
}

rownames(rna) = sapply(rownames(rna), change)
rownames(pcg) = sapply(rownames(pcg), change)

change = function(rowname){
  new = norm_conversion$id[which(norm_conversion$TCGA_id %in% rowname)]
  return(new)  
}
rownames(norm) = sapply(rownames(norm), change)

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(rna) %in% clin$bcr_patient_barcode)
rna = rna[z,] #all have clinical data - 7387 patients 
#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(pcg) %in% clin$bcr_patient_barcode)
pcg = pcg[z,] #all have clinical data - 7387 patients 
z <- which(rownames(norm) %in% clin$bcr_patient_barcode)
norm = norm[z,] #all have clinical data - 563 patients (matched normals)

#Add survival info to RNA file
rna = as.data.frame(rna)
rna$patient = rownames(rna)
colnames(clin)[2] = "patient"
rna = merge(rna, clin, by="patient")

saveRDS(rna, "5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")

#Add survival info to RNA file - normal matched 
norm = as.data.frame(norm)
norm$patient = rownames(norm)
colnames(clin)[2] = "patient"
norm = merge(norm, clin, by="patient")

saveRDS(norm, "all_genes_matched_normals_563_March13_wclinical_data.rds")

#Add survival info to PCG file
pcg = as.data.frame(pcg)
pcg$patient = rownames(pcg)
pcg = merge(pcg, clin, by="patient")

saveRDS(pcg, "19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")


























































