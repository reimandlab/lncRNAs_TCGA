#rnaseq_files_into_matrix.R 

library(data.table)
library(plyr)
library(dplyr)

#get list of files, have 9,246 rna-seq files in total

files <- list.files()

#need to open each file folder --> unzip file --> readfile --> 

temp = files

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

saveRDS(named.list, file="11052rnaSEQfilesLIST.rds")
print("finished saving/making list")

new <- named.list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="V1"), .)

saveRDS(new, file="11052rnaSEQfiles.rds")

#columns of files.matrix = order of metafast files UUIDS







