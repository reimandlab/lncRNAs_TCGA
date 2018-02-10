###---------------------------------------------------------------
###TCGA_cancers_data1.R
###---------------------------------------------------------------

###October 10th, 2017
###Goal: Get distribution of patients in each cancer type with
#RNA data 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. RNA data 

#--lncrna
rna = readRDS("5919_lncLiverOvariancancers_TCGAnew.rds")
rownames(rna) = rna$gene
rna$gene = NULL
rna = t(rna)

pcg = readRDS("54564_lncLiverOvariancancers_TCGAnew.rds")
rownames(pcg) = pcg$gene
pcg$gene = NULL
pcg = t(pcg)

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

#genes = c("AC006126.4", "ADORA2A-AS1", "NEAT1", "LINC00665", "ZNF503-AS2", "GS1-251I9.4")
#z1 = which(fantom[,2] %in% genes)

#Subset to intergenic and antisense 
#z2 = which(fantom$CAT_geneClass %in% c("lncRNA_antisense", "lncRNA_intergenic"))

#fantom = fantom[c(z1,z2),]

#5. TCGA ID cancer type conversion 
canc_conversion = readRDS("Liver_Ovarian_tcga_id_cancer_type_conversion.txt")

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

saveRDS(rna, "lnc_rna_ovary_liver_plus_clinical.rds")
saveRDS(pcg, "pcg_rna_ovary_liver.rds")


