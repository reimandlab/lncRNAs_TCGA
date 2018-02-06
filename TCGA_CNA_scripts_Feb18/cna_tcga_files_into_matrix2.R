#rnaseq_files_into_matrix2.R 

library(data.table)
library(plyr)
library(dplyr)

named.list <- readRDS(file="Panc_Liv_Ov_KIRC_list_of_CNA_data_files.rds")

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

#combine list of datafiles into one dataframe
new = do.call(rbind.data.frame, named.list)

saveRDS(new, file="OvaryLiverPancreasKIRC_CNA_TCGA_files.rds")









