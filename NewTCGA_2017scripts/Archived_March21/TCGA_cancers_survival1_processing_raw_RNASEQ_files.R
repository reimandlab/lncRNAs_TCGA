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
rna = readRDS("5919_lncs4cancers_TCGAnew.rds")
rownames(rna) = rna$gene
rna$gene = NULL
rna = t(rna)

#--pcg
pcg = readRDS("54564_PCGS4cancers_TCGAnew.rds")
rownames(pcg) = pcg$gene
pcg$gene = NULL
pcg = t(pcg)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% "protein_coding")
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#keep only genes in pcg that are labelled as PCG
z <- which(colnames(pcg) %in% ucsc$hg19.ensGene.name2)
pcg = pcg[,z]

#2. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

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

#6. List of TCGA IDs used in PCAWG - to remove
ids_remove = fread("TCGA_IDs_usedinPCAWG.txt")

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

#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) #600 to remove 
z <- which(rownames(rna) %in% ids_remove) #666 PCAWG samples in this TCGA RNA file
rna = rna[-z,]
z <- which(rownames(pcg) %in% ids_remove) #666 PCAWG samples in this TCGA RNA file
pcg = pcg[-z,]

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(rna) %in% clin$bcr_patient_barcode)
rna = rna[z,] #all have clinical data - 7387 patients 
#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(pcg) %in% clin$bcr_patient_barcode)
pcg = pcg[z,] #all have clinical data - 7387 patients 

#Add survival info to rna file
rna = as.data.frame(rna)
rna$canc = "" 
rna$time = ""
rna$status = ""
rna$sex = ""

for(i in 1:dim(rna)[1]){
	pat = rownames(rna)[i]
	z <- which(clin$bcr_patient_barcode == pat)
	if(length(z) >1){z = z[1]}
	status = clin$vital_status[z]
	time = clin$days_to_death[z]
	if(is.na(time)){
		time = clin$days_to_last_followup[z]
	}
	rna$status[i] = status
	rna$time[i] = time
	rna$sex[i] = clin$gender[z]
}

for(i in 1:dim(rna)[1]){
	pat = rownames(rna)[i]
	z = which(canc_conversion$id == pat)
	canc = canc_conversion$Cancer[z]
	rna$canc[i] = canc
}

saveRDS(rna, "5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

#Add survival info to rna file
pcg = as.data.frame(pcg)
pcg$canc = "" 
pcg$time = ""
pcg$status = ""
pcg$sex = ""

for(i in 1:dim(pcg)[1]){
  pat = rownames(pcg)[i]
  z <- which(clin$bcr_patient_barcode == pat)
  if(length(z) >1){z = z[1]}
  status = clin$vital_status[z]
  time = clin$days_to_death[z]
  if(is.na(time)){
    time = clin$days_to_last_followup[z]
  }
  pcg$status[i] = status
  pcg$time[i] = time
  pcg$sex[i] = clin$gender[z]
}

for(i in 1:dim(pcg)[1]){
  pat = rownames(pcg)[i]
  z = which(canc_conversion$id == pat)
  canc = canc_conversion$Cancer[z]
  pcg$canc[i] = canc
}

saveRDS(pcg, "19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

































































