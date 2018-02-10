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

genes = c("AC006126.4", "ADORA2A-AS1", "NEAT1", "LINC00665", "ZNF503-AS2", "GS1-251I9.4")
z1 = which(fantom[,2] %in% genes)

#Subset to intergenic and antisense 
z2 = which(fantom$CAT_geneClass %in% c("lncRNA_antisense", "lncRNA_intergenic"))

fantom = fantom[c(z1,z2),]

#4. List of lncRNA survival associated candidates 
#cands = fread("7tier1_35tier2_lncRNA_candidates_August28th.txt")

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

#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) #600 to remove 
z <- which(rownames(rna) %in% ids_remove) #666 PCAWG samples in this TCGA RNA file
rna = rna[-z,]

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(rna) %in% clin$bcr_patient_barcode)
rna = rna[z,] #all have clinical data - 7387 patients 
z <- which(colnames(rna) %in% fantom$CAT_geneID)
rna = rna[,z]

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

##Only keep cancer patients corresponding to candidate lncRNAs 
rna = subset(rna, rna$canc %in% cancers)

#Sepetate into candidate-cancer datasets 
#Ovarian 
ov <- c("LINC00665", "ZNF503-AS2", "GS1-251I9.4")
ovary = subset(rna, rna$canc == "Ovarian serous cystadenocarcinoma")
genes = fantom[which(fantom[,2] %in% ov),1]
ovary = ovary[,c(7:10, which(colnames(ovary) %in% genes))]

liv = c("NEAT1", "ADORA2A-AS1")
liver = subset(rna, rna$canc == "Liver hepatocellular carcinoma")
genes = fantom[which(fantom[,2] %in% liv),1]
liver = liver[,c(7:10, which(colnames(liver) %in% genes))]

kid = c("AC006126.4")
kidney = subset(rna, rna$canc == "Kidney renal clear cell carcinoma")
genes = fantom[which(fantom[,2] %in% kid),1]
kidney = kidney[,c(7:10, which(colnames(kidney) %in% genes))]

#put together 
cands = list(ovary, liver, kidney)

#add high/low tag 
add_tags = function(data){
	genes = length(5:ncol(data))
	for(i in 1:genes){
		gene = colnames(data)[4+i]
		med = median(data[,4+i])
		for(y in 1:nrow(data)){
			check = data[y,4+i]
			if(check >=med){
				data[y,4+i] = "high"
			}
			if(check < med){
				data[y,4+i] = "low"
			}
		}
	}
	return(data)
}

cands = llply(cands, add_tags)








