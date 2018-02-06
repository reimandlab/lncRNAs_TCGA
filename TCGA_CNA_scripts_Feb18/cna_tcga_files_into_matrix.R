#rnaseq_files_into_matrix.R 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#get list of files, have 167 rna-seq files in total

files <- list.files()

#need to open each file folder --> unzip file --> readfile --> 

temp = files

#for each f of list
readALLfiles <- function(f){	
	open <- list.files(f)
	#only want the ones that have "nocnv_grch38.seg.txt in the end" 
	check = which(str_detect(open, "nocnv_grch38.seg.txt")) >= 1

	if(!(length(open)==0)){

	if(check){

	if(length(open) >1){
		open <- open[which(str_detect(open, "nocnv_grch38.seg.txt"))]
	}
	open <- paste(getwd(), f, open, sep="/")
	o <- read.table(open, sep="\t", header=T)
	colnames(o)[1] = "patient"
	o[,1] = f
	print(f)
	return(o)
}
}
}

#apply to each file in list
named.list <- llply(temp, readALLfiles, .progress = "text")

saveRDS(named.list, file="Panc_Liv_Ov_KIRC_list_of_CNA_data_files.rds")
print("finished saving/making list")

