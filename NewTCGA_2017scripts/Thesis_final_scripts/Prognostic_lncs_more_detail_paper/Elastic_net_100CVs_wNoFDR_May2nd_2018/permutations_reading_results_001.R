
#ran elastic net simulations for each cancer types 100 times 
#this means 100 runs of 100 elastic net cross-validations 
#should be 100 files per cancer types if all ran correctly 

#results files are here: 
#/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/elastic_net_permutations_results

library(data.table)
library(dplyr)
library(stringr)

#get cancer names 

all_files = list.files()

f = list.files()[[1]]
unlist(strsplit(f, "_"))[1]

cancers = unique(sapply(all_files, function(x){unlist(strsplit(x, "_"))[1]}))

#seperate files into cancer types




