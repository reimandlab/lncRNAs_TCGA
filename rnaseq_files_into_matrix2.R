#rnaseq_files_into_matrix2.R 

library(data.table)
library(plyr)
library(dplyr)

named.list <- readRDS(file="3228rnaSEQfilesLIST.rds")

new <- named.list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="V1"), .)

saveRDS(new, file="3228rnaSEQfiles.rds")

#columns of files.matrix = order of metafast files UUIDS


