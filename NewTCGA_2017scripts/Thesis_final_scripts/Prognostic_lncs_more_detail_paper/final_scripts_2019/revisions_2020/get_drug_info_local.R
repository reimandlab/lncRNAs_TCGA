library(TCGAbiolinks)
library(data.table)
library(plyr)

cancers = c("TCGA-COAD", "TCGA-LGG", "TCGA-BRCA", "TCGA-READ", "TCGA-LUSC", 
"TCGA-LUAD", "TCGA-ACC", "TCGA-COAD", "TCGA-OV", "TCGA-GBM",
"TCGA-COAD", "TCGA-LIHC", "TCGA-KIRC", "TCGA-KIRP",
"TCGA-MESO", "TCGA-STAD", "TCGA-SARC", "TCGA-BLCA", "TCGA-HNSC", "TCGA-UCEC",
"TCGA-ESCA")

get_drug_info = function(canc){

query <- GDCquery(project = canc, 
              data.category = "Clinical", 
              file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "drug")
clinical = as.data.table(clinical)
clinical$cancer = unlist(strsplit(canc, "-"))[2]
return(clinical)
}

drug_info = llply(cancers, get_drug_info, .progress="text")
saveRDS(drug_info, file="all_cancers_drug_info.rds")
