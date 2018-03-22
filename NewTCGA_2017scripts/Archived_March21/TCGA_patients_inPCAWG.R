###TCGA_patients_inPCAWG.R
###---------------------------------------------------------------

###October 10th, 2017
###Goal: Identify list of patients that are in the PCAWG RNA-Seq 
#and clinical files that also have a TCGA ID so we can exclude 
#them from further TCGA analysis since they were also studied 

###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###Load Data 
###---------------------------------------------------------------

#PCAWG Clinical file that has patient TCGA ID if that's where they 
#also from

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#PCAWG RNA-seq file to see which patients were actually studied 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

###"TUMOUR SAMPLES" - PROCESSING RNA FILE 
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z] 

##Convert RNA-Seq IDs to ICGC donor IDs
for(i in 1:ncol(rna)){
	z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna)[i])
	if(!(length(z)==0)){
		colnames(rna)[i] <- conversion$icgc_donor_id[z]
	}
}

#Change gene names to match ensembl IDs 
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1:2], 1, ensg)
z <- which(check==TRUE)
rna <- rna[z,]

##Keep only patient samples that also have clinical data
z <- which(colnames(rna) %in% clin$icgc_donor_id)
rna <- rna[,z]

#Get TCGA ids of these 1267 patients that have both clinical and RNA-Seq data in PCAWG
pats = colnames(rna)
tcga_ids = clin$tcga_donor_uuid[which(clin$icgc_donor_id %in% pats)]
tcga_ids = tcga_ids[which(!(tcga_ids==""))]
tcga_ids = unique(tcga_ids)

write.table(tcga_ids, "819_unique_TCGAids_usedbyPCAWG.txt", quote=F) #remove these from TCGA analysis 