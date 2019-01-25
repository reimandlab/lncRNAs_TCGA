source("source_file.R")
library(stringr)
library(VennDiagram)
library(EnvStats)

combined = readRDS("ICs_wSurvival_GTEx_info_Jan24.rds")

#1. add gene names 

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"
colnames(ucsc)[6] = "gene"

change_name = function(name){
  z = which(ucsc$gene == name)
  return(ucsc$HGNC.symbol[z])
}

combined$gene_name = unlist(llply(combined$gene, change_name))
combined$combo = NULL

#2. order by prognosis p-value and then by median difference between cancer and normal tissues 

combined = combined[order(pval)]

#3. summarize only the ones whose HR matches gtex signal 
combined = as.data.table(filter(combined, !(biological_match == ""))) #282 pairs left, 81 unique ion channels 
combined$coef = NULL
combined$ic_test_ph = NULL
combined$global_test_ph = NULL
combined$fc_mean = NULL
combined$cancer = NULL

write.csv(combined, file="ICs_sig_HR_GTEX_matches_Jan24.csv", quote=F, row.names = F)
