options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)

date = Sys.Date()

setwd("Downloads")

#data
f = fread("GSE89567_IDH_A_processed_data.txt")

#how amny are protein coding genes
library(annotables)
gene_data = as.data.table(grch37)

genes = unique(unlist(f[,1]))
rownames(f) = genes
f$V1 = NULL

check_gene_type = function(gene){
  print(gene)
  gene = unlist(strsplit(gene, "'"))[2]
  genetype = filter(gene_data, symbol == gene)$biotype[1]
  return(c(gene, genetype))
}
all_genes_types = as.data.table(ldply(llply(genes, check_gene_type, .progress="text")))

get_perc = function(gene){
  print(gene)
  gene_dat = which(rownames(f) == gene)
  gene_dat = f[gene_dat,]
  gene_dat = unlist(gene_dat)

  #how many non zeroes
  gene_dat = gene_dat[which(gene_dat >0)]
  num_cells = length(gene_dat)
  med_exp = median(gene_dat)

  return(c(gene, num_cells, med_exp))
}

all_genes = as.data.table(ldply(llply(genes, get_perc, .progress="text")))
colnames(all_genes) = c("gene", "num_cells_wexp", "median_exp")
all_genes = filter(all_genes, median_exp >0)
all_genes = all_genes[order(median_exp)]
all_genes$rank = 1:nrow(all_genes)
all_genes$perc_exp = all_genes$rank/nrow(all_genes)
all_genes$perc_cells_exp = as.numeric(all_genes$num_cells_wexp)/6341
all_genes = all_genes[order(perc_cells_exp)]
all_genes$rank_num_cells = 1:nrow(all_genes)
all_genes$perc_num_cells = all_genes$rank_num_cells/nrow(all_genes)
