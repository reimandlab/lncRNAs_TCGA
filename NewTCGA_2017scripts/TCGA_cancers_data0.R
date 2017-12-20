###---------------------------------------------------------------
###TCGA_cancers_data0.R
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

#1. List of TCGA IDs used in PCAWG - to remove
ids_remove = fread("TCGA_IDs_usedinPCAWG.txt")

#2. TCGA Tumour Codes Table
tss_codes = read.csv(" TCGA_TissueSourceSite_Codes2017 .csv"     )

#3. TCGA new clinical file
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

#4. RNA-Seq File
rna = readRDS("9246rnaSEQfiles.rds") #doesn't include GBM
rna[,1] = as.character(rna[,1])

#5. Fantom data 
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

#6. List of lncRNA survival associated candidates 
#cands = fread("7tier1_35tier2_lncRNA_candidates_August28th.txt")

#7. TCGA source codes
source = read.csv("TCGA_sample_codes.csv")

###---------------------------------------------------------------
###Process Data 
###---------------------------------------------------------------

#1. Gene IDs 
extract3 <- function(row){
	gene <- as.character(row[[2]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}

rna$gene = rna[,1]
rna[,9248] <- apply(rna[,9247:9248], 1, extract3) ; 
rna$V1 = NULL
rna = as.data.table(rna)

#2. Convert UUIDs into TCGA ids
barcode = readRDS("TCGA_barcode.RDS")
barcode$uuid = as.character(barcode$uuid)
barcode$bar = as.character(barcode$bar)

change = function(colname){
	z <- which(barcode$uuid == colname)
	return(barcode$bar[z])	
}

colnames(rna)[1:9246] = sapply(colnames(rna)[1:9246], change)

#3. Seperate into cancers 
get_cancer = function(tcga_id){
	tss = unlist(strsplit(tcga_id, "-"))[2]
	canc = tss_codes[which(tss_codes[,1] %in% tss), 3]
	if(length(canc)==0){canc = ""}
	canc = c(tcga_id, canc)
	return(canc)
}

cancers_order = unlist(llply(colnames(rna)[1:9246], get_cancer))
cancers_order = matrix(cancers_order, ncol=2, byrow=T)
colnames(cancers_order) = c("TCGA_id", "Cancer")
cancers_order = as.data.table(cancers_order)

#cancers_keep = c("Ovarian serous cystadenocarcinoma", "Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma" , "Pancreatic adenocarcinoma")
#cancers_keep = filter(cancers_order, Cancer %in% cancers_keep)

#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) 
cancers_keep = cancers_order
cancers_keep$source = ""
for(i in 1:nrow(cancers_keep)){
	source = unlist(strsplit(as.character(cancers_keep[i,1]), '-'))[4]
	source = substr(source, 1, 2)
	cancers_keep$source[i] = source
}

#Keep only tumours for now, save normal for later
cancers_keep = cancers_keep[cancers_keep$source == "11", ]

clean_tcga_id = function(row){
	id = row[[1]]
	s1 = unlist(strsplit(id, '-'))[1]
	s2 = unlist(strsplit(id, '-'))[2]
	s3 = unlist(strsplit(id, '-'))[3]
	return(paste(s1, s2, s3, sep="-"))
}

cancers_keep$id = ""
cancers_keep[,4] = apply(cancers_keep, 1, clean_tcga_id) 

#remove duplciated patient samples = all kidney clear cell for some reason
dups = cancers_keep[which(duplicated(cancers_keep[,4])),4]
z <- which(cancers_keep$id %in% dups$id)
cancers_keep = cancers_keep[-z,] #1406 unique TCGA tumour IDs 

table(cancers_keep$Cancer)

#Kidney renal clear cell carcinoma    Liver hepatocellular carcinoma 
#                              526                               371 
#Ovarian serous cystadenocarcinoma         Pancreatic adenocarcinoma 
#                              332                               177 

saveRDS(cancers_keep, file="tcga_id_cancer_type_conversion.txt")

#3. Subset RNA file 
z <- which(colnames(rna) %in% cancers_keep$TCGA_id)
rna = as.data.frame(rna)
rna = rna[,c(z,9247)]

#4. Keep only lncRNA genes 
z <- which(rna$gene %in% fantom$CAT_geneID)
lnc_rna = rna[z,]
saveRDS(lnc_rna, "5919_lncs4matched_normal_tissues_TCGAnew.rds")

#everything else (includes PCGS, other lncRNAs not in FANTOM)
pcg_rna = rna[-z,]
saveRDS(pcg_rna, "54564_PCGS4cancers_TCGAnew.rds")













