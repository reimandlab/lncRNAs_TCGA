setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2018_HTSEQ_counts")

library(data.table)
library(plyr)
library(dplyr)

named.list <- readRDS(file="9320_raw_countsLIST.rds")
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

saveRDS(new, file="9320rnaSEQfiles.rds")

#columns of files.matrix = order of metafast files UUIDS
