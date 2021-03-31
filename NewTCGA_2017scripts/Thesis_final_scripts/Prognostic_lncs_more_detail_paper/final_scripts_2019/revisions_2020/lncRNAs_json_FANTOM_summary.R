options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)
library("rjson")

date = Sys.Date()

setwd("/Users/kisaev/Documents/lncRNAs/Jan2021/json_files")

#figure 2 - univaraite
res = readRDS("/Users/kisaev/Documents/lncRNAs/Jan2021/final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#json files
jsons = list.files(pattern=".json")

#function that reads each file
get_json_res = function(file_name){
  print(file_name)
  read_file = fromJSON(file = file_name)
  read_file_sum = as.data.table(read_file$summary)
  gene = unlist(strsplit(file_name, "\\."))[1]
  read_file_sum$gene = gene
  read_file_sum$gene_name = res$gene_symbol[which(res$gene == gene)][1]
  read_file_sum$Pubmed = paste(read_file$Pubmed, collapse=";")
  read_file_sum = cbind(read_file_sum,as.data.table(read_file$External_ID))
  read_file_sum$FANTOM_Gene_Class = read_file$FANTOM_Gene_Class
  read_file_sum$Synonyms = paste(read_file$Synonyms, collapse=",")
  print(head(read_file_sum))
  return(read_file_sum)
}

all_genes = as.data.table(ldply(llply(jsons, get_json_res)))
all_genes$eQTL_mRNA_coexpression[all_genes$eQTL_mRNA_coexpression == "No implications in coexpression with eQTL linked mRNAs."] = NA
all_genes$sample_ontology_association[all_genes$sample_ontology_association == "Not associated with any sample ontology."] = NA
all_genes$trait_association[all_genes$trait_association == "Not associated with any trait."] = NA
all_genes$trait_vs_sample_ontology_enrichment[all_genes$trait_vs_sample_ontology_enrichment == "No implications in trait vs sample ontology pairs."] = NA

library(dplyr)
all_genes = all_genes %>%
  mutate(across(everything(), as.character))

write.table(all_genes, file="lncRNA_fantom_annotations.csv", quote=F, row.names=F, sep="}")

library(easyPubMed)
library(httr)

all_genes = unique(res$gene_symbol)

get_pubmed_info = function(gene){
  print(gene)
  myQuery = gene
  myIdList <- get_pubmed_ids(myQuery)
  if(!(myIdList$Count == "0")){
  y <- fetch_pubmed_data(pubmed_id_list = myIdList)
  yy <- table_articles_byAuth(y, included_authors = "first", getKeywords = TRUE)
  yy = yy[, c("pmid", "jabbrv", "title", "abstract", "year")]
  z = which(str_detect(yy$abstract, gene))
  yy=yy[z,]
  if(!(dim(yy)[1] == 0)){
  yy$gene = res$gene[res$gene_symbol == gene][1]
  yy$gene_symbol = gene
  return(yy)
  }
}}

all_pubmeds = as.data.table(ldply(llply(all_genes, get_pubmed_info)))
write.table(all_pubmeds, file="lncRNA_pubmed_annotations.csv", quote=F, row.names=F, sep="}")
