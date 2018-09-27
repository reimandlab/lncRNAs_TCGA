dim(pcg_rna)

pcg_rna = t(pcg_rna)
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient = ""
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
#pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")

saveRDS(pcg_rna, "all_mrna_may8th.rds")


dim(lnc_rna)

lnc_rna = t(lnc_rna)
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient = ""
lnc_rna$patient <- rownames(lnc_rna)

#remove duplicated column names 
dups <- colnames(lnc_rna)[which(duplicated(colnames(lnc_rna)))]   
#save them in a list for future reference 
#pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")

saveRDS(lnc_rna, "all_lncRNA_may8th.rds")