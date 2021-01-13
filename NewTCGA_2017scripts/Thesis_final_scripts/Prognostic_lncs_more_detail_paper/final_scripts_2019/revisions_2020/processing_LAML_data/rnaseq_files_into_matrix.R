#rnaseq_files_into_matrix.R

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

setwd("/u/kisaev/TCGA_LAML_count_data")
files <- list.files()

#need to open each file folder --> unzip file --> readfile -->
temp = files
z = which(str_detect(temp, "Counts"))
#temp = temp[-z]
z = which(str_detect(temp, "rna"))
#temp = temp[-z]
z = which(str_detect(temp, "tar"))
temp = temp[-z]
z = which(str_detect(temp, ".txt"))
temp = temp[-z]

#for each f of list
readALLfiles <- function(f){
        open <- list.files(f)
        if(!(length(open)==0)){

        if(length(open) >1){
                g <- grep("gz", open)
                open <- open[g]
        }
        open <- paste(getwd(), f, open, sep="/")
        o <- read.table(open, sep="\t")
        colnames(o)[2] = f
        print(f)
        return(o)
}
}

#apply to each file in list
named.list <- llply(temp, readALLfiles, .progress = "text")

saveRDS(named.list, file="151_LAML_raw_countsLIST.rds")
print("finished saving/making list")

named.list = readRDS("151_LAML_raw_countsLIST.rds")

new <- named.list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="V1"), .)

saveRDS(new, file="151_LAML_counts_rnaSEQfiles.rds")
