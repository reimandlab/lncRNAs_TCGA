#get_TCGA_ids_fromPCAWG.R

library(data.table)
library(stringr)
library(dplyr)
library(plyr)

ids = clin$submitted_sample_id

#all clin file from TCGA 
tcga = read.csv("all_clin_XML-3.csv")

z <- which(str_detect(ids, "TCGA"))

ids_tcgaipcawg = ids[z]

change = function(id){
	sep = unlist(strsplit(id, "-"))
	tog = paste(sep[1], sep[2], sep[3], sep = "-")
	return(tog)
}

ids_tcgaipcawg = unique(unlist(llply(ids_tcgaipcawg, change)))

z <- which(tcga$bcr_patient_barcode %in% ids_tcgaipcawg)
keep = tcga[z, c(1,5)]
z = which(duplicated(keep[,1]))
keep = keep[-z,]

write.table(keep, file="TCGA_IDs_usedinPCAWG.txt", quote=F, row.names=F)




