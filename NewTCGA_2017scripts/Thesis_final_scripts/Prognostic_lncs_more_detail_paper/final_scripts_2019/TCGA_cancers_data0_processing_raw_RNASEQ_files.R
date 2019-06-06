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
#ids_remove = readRDS("TCGA_IDs_usedinPCAWG.rds")
#ids_remove[,2] = sapply(ids_remove[,2], function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
#ids_remove = as.data.frame(ids_remove)
#ids_remove = unique(ids_remove[,2])

#2. TCGA Tumour Codes Table
tss_codes = read.csv(" TCGA_TissueSourceSite_Codes2017 .csv"     )

#not actually done originally 
tss_codes$TSS.Code[tss_codes$TSS.Code == "2"] = "02"
tss_codes$TSS.Code[tss_codes$TSS.Code == "6"] = "06"
tss_codes$TSS.Code[tss_codes$TSS.Code == "8"] = "08"
tss_codes$TSS.Code[tss_codes$TSS.Code == "1"] = "01"
tss_codes$TSS.Code[tss_codes$TSS.Code == "4"] = "04"
tss_codes$TSS.Code[tss_codes$TSS.Code == "3"] = "03"
tss_codes$TSS.Code[tss_codes$TSS.Code == "9"] = "09"
tss_codes$TSS.Code[tss_codes$TSS.Code == "5"] = "05"

#3. TCGA new clinical file - downloaded previously 
#clin = read.csv("all_clin_XML_tcgaSept2017.csv")
#clin = clin[,1:90]

#3B. New clinical data from CELL paper PANCAN release 
clin = fread("mmc1_clinical_data_cellpaper2018.txt")

#4. RNA-Seq File
rna = readRDS("11052rnaSEQfiles.rds") 
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
rna[,ncol(rna)] <- apply(rna[,(ncol(rna)-1):ncol(rna)], 1, extract3) ; 
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

colnames(rna)[1:(ncol(rna)-1)] = sapply(colnames(rna)[1:(ncol(rna)-1)], change)

#3. Seperate into cancers 
get_cancer = function(tcga_id){
	tss = unlist(strsplit(tcga_id, "-"))[2]
	canc = tss_codes[which(tss_codes[,1] %in% tss), 3]
	if(length(canc)==0){canc = ""}
	canc = c(tcga_id, canc)
	return(canc)
}

cancers_order = unlist(llply(colnames(rna)[1:(ncol(rna)-1)], get_cancer))
cancers_order = matrix(cancers_order, ncol=2, byrow=T)
colnames(cancers_order) = c("TCGA_id", "Cancer")
cancers_order = as.data.table(cancers_order)

#cancers_keep = c("Ovarian serous cystadenocarcinoma", "Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma" , "Pancreatic adenocarcinoma")
#cancers_keep = filter(cancers_order, Cancer %in% cancers_keep)

cancers_keep = cancers_order
cancers_keep$source = ""
for(i in 1:nrow(cancers_keep)){
	source = unlist(strsplit(as.character(cancers_keep[i,1]), '-'))[4]
	source = substr(source, 1, 2)
	cancers_keep$source[i] = source
}

#Keep only tumours for now, save normal for later
normals_keep = cancers_keep[cancers_keep$source == "11", ]
metastatic_keep = cancers_keep[cancers_keep$source == "06", ]
cancers_keep = cancers_keep[cancers_keep$source == "01", ]

clean_tcga_id = function(row){
	id = row[[1]]
	s1 = unlist(strsplit(id, '-'))[1]
	s2 = unlist(strsplit(id, '-'))[2]
	s3 = unlist(strsplit(id, '-'))[3]
	return(paste(s1, s2, s3, sep="-"))
}

cancers_keep$id = ""
cancers_keep[,4] = apply(cancers_keep, 1, clean_tcga_id) 
metastatic_keep$id = ""
metastatic_keep[,4] = apply(metastatic_keep, 1, clean_tcga_id)
normals_keep$id = ""
normals_keep[,4] = apply(normals_keep, 1, clean_tcga_id) 

#remove those patients already used in PCAWG
#ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) 
#z = which(cancers_keep$id %in% ids_remove)
#external_dataset = cancers_keep[z,]
#cancers_keep = cancers_keep[-z,]

#z = which(normals_keep$id %in% ids_remove)
#normals_keep = normals_keep[-z,]

#z = which(metastatic_keep$id %in% ids_remove)
#metastatic_keep = metastatic_keep[-z,]

#remove duplciated patient samples = all kidney clear cell for some reason
dups = cancers_keep[which(duplicated(cancers_keep[,4])),4]
z <- which(cancers_keep$id %in% dups$id)
cancers_keep = cancers_keep[-z,] #7564 unique TCGA tumour IDs 
z <- which(normals_keep$id %in% dups$id)
normals_keep = normals_keep[-z,] #563 unique TCGA normal IDs 
#z = which(external_dataset$id %in% dups$id)

table(cancers_keep$Cancer)
#table(external_dataset$Cancer)

z = which(cancers_keep$Cancer == "") #remove those patients with no cancer type
cancers_keep = cancers_keep[-z,] #7501 samples in total with both RNA_Sequencing and clinical data 

#z = which(external_dataset$Cancer == "") #remove those patients with no cancer type
#external_dataset = external_dataset[-z,] #7501 samples in total with both RNA_Sequencing and clinical data 

saveRDS(cancers_keep, file="lncRNAs_2019_manuscript/tcga_id_cancer_type_conversion.txt")
saveRDS(normals_keep, file="lncRNAs_2019_manuscript/tcga_id_NORMAL_samples_type_conversion.txt")
saveRDS(metastatic_keep, file="lncRNAs_2019_manuscript/tcga_id_Metastatic_samples_type_conversion.txt")
#saveRDS(external_dataset, file="lncRNAs_2019_manuscript/tcga_id_external_samples_type_conversion.txt")

#3. Subset RNA file 
z <- which(colnames(rna) %in% normals_keep$TCGA_id)
norm = as.data.frame(rna)
norm = norm[,c(z,ncol(norm))]

z <- which(colnames(rna) %in% metastatic_keep$TCGA_id)
rna = as.data.frame(rna)
met = rna[,c(z,ncol(rna))]

#z <- which(colnames(rna) %in% external_dataset$TCGA_id)
#external_data = rna[,c(z,ncol(rna))]

z <- which(colnames(rna) %in% cancers_keep$TCGA_id)
rna = rna[,c(z,ncol(rna))]

#4. Keep only lncRNA genes 
z <- which(rna$gene %in% fantom$CAT_geneID)
lnc_rna = rna[z,]
saveRDS(lnc_rna, "lncRNAs_2019_manuscript/5919_all_tumours_9603_tissues_TCGAnew.rds")

#everything else (includes PCGS, other lncRNAs not in FANTOM)
pcg_rna = rna[-z,]
saveRDS(pcg_rna, "lncRNAs_2019_manuscript/54564_PCGs_all_tumours_9603_tissues_TCGAnew.rds")

#4. normal patients --> lncRNA and pcg data 
saveRDS(norm, "lncRNAs_2019_manuscript/all_genes_563_matched_normal_samples_TCGA_April11.rds")
saveRDS(met, "lncRNAs_2019_manuscript/all_genes_354_matched_metastatic_tumours_TCGA_april.rds")

#saveRDS(external_data, "lncRNAs_2019_manuscript/external_data_set_pcawg_patients_all_tumours_7501_tissues_TCGAnew.rds")










