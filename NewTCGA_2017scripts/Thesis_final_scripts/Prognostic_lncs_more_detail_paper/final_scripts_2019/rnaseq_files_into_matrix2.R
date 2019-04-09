#rnaseq_files_into_matrix2.R 

library(data.table)
library(plyr)
library(dplyr)

named.list <- readRDS(file="9314_raw_countsLIST.rds")
print("finished reading file")

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
if(!(length(z)==0)){
	named.list = named.list[-z]
}

print(paste(length(z), "removed"))

new <- named.list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="V1"), .)

print("got new list")

saveRDS(new, file="9246rnaSEQfiles.rds")

#columns of files.matrix = order of metafast files UUIDS








