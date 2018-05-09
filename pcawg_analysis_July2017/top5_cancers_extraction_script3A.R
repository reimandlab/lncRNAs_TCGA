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

saveRDS(pcg_rna, "all_rna_may8th.rds")
