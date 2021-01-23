library(data.table)
library(plyr)
library(dplyr)
library(stringr)

setwd("/u/kisaev/TCGA_SKCM_count_data")

#rna count matrix
rnas=readRDS("472_SKCM_counts_rnaSEQfiles.rds")

#ids conversion
ids=fread("File_metadata.txt")

#change ids in rna matrix
colnames(rnas)[2:ncol(rnas)]=sapply(colnames(rnas)[2:ncol(rnas)], function(x){filter(ids,
  file_id==x)$cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id})

#make sure all samples are tumour samples
table(sapply(colnames(rnas)[2:ncol(rnas)], function(x){unlist(strsplit(x, "-"))[4]}))

#01A 06A 06B 07A 11A
#103 365   2   1   1

#only keep th 01As
cols=colnames(rnas)[which(str_detect(colnames(rnas), "01A"))]
z=which(colnames(rnas) %in% c("V1", cols))
rnas=rnas[,z]

colnames(rnas)[2:ncol(rnas)] = sapply(colnames(rnas)[2:ncol(rnas)], function(x){
  paste(unlist(strsplit(x, "-"))[1:3], collapse="-")
})

colnames(rnas)[1]="gene"
z = which(str_detect(rnas$gene, "ENSG"))
rnas=rnas[z,]
rnas$gene=as.character(rnas$gene)
rnas$gene= sapply(rnas$gene, function(x){unlist(strsplit(x, "\\."))[1]})
rownames(rnas) = rnas$gene
rnas$gene=NULL
rnas=t(rnas)
rnas=as.data.frame(rnas)
rnas$patient = rownames(rnas)
rnas$type="SKCM"

saveRDS(rnas, file="SKCM_count_data.rds")
