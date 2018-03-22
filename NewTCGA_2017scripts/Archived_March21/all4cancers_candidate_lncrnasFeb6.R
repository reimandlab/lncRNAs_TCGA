###combined results from all cancers

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)

#liver
liver = readRDS("final_candidates_10batches_of1000CV_LIHC_cancer_patientsFeb2.RDS")
liver$canc = "liver"
keep =  as.data.table(table(liver$gene))
keep = dplyr::filter(keep, N >=5)
liver = liver[which(liver$gene %in% keep$V1),]

#kidney
kidney = readRDS("final_candidates_10batches_of1000CV_KIRC_cancer_patientsJan29.RDS")
kidney$canc = "kidney"
keep =  as.data.table(table(kidney$gene))
keep = dplyr::filter(keep, N >=5)
kidney = kidney[which(kidney$gene %in% keep$V1),]

#ovary
ovary = readRDS("final_candidates_10batches_of1000CV_305ovarian_cancer_patientsJan23.RDS")
ovary$canc = "ovary"
keep =  as.data.table(table(ovary$gene))
keep = dplyr::filter(keep, N >=5)
ovary = ovary[which(ovary$gene %in% keep$V1),]

#pancreas
pancreas = readRDS("final_candidates_10batches_of1000CV_PAAD_cancer_patientsFeb6.RDS")
pancreas$canc = "pancreas"
keep =  as.data.table(table(pancreas$gene))
keep = dplyr::filter(keep, N >=5)
pancreas = pancreas[which(pancreas$gene %in% keep$V1),]

all = rbind(liver, kidney, ovary, pancreas) #all mutually exclusive 
saveRDS(all, file="36_unique_cands_4cancers_TCGA_Feb6.rds")
