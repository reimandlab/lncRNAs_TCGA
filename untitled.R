#rnaseq_files_into_matrix2.R 

library(data.table)
library(plyr)
library(dplyr)

named.list <- readRDS(file="167rnaSEQfilesLIST.rds")

#remove datafiles with nothing in them not sure why
#it's there in the first place 

find_0 = function(dat){
	d = length(dat)
	if(d ==0){
		return("remove")
	}
}

remove = llply(named.list, find_0)
z <- which(remove == "remove")
named.list = named.list[-z]

new <- named.list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="V1"), .)

saveRDS(new, file="167rnaSEQfiles.rds")

#columns of files.matrix = order of metafast files UUIDS








