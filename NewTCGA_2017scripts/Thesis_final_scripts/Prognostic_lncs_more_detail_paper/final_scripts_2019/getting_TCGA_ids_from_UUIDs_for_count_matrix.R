setwd("/Users/kisaev/remote22") #<- "/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ"

library(GenomicDataCommons)

#get UUIDs used for RNA-Seq files
rna = readRDS("9246rnaSEQfiles.rds") 
uuids = colnames(rna)

library(GenomicDataCommons)
library(magrittr)

TCGAtranslateID = function(file_ids, legacy = FALSE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

res = TCGAtranslateID(uuids)
saveRDS(re, file="File_metadata.txt")