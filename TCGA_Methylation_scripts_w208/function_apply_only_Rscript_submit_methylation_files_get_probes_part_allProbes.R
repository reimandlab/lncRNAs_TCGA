###methylation files 
###identify probes overlapping lncRNAs 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(doParallel)
library(parallel)


#############in R###################################################################################################################################


files_first = as.list(list.files(pattern="patients.txt"))
files_second = as.list(list.files(pattern="data.data.txt"))

transform_files = function(first_file, second_file){
	pats = read.table(first_file)
	pats = pats[,3:ncol(pats)]
	pats = as.character(unlist(pats))
	
	library(data.table)

	meth = fread(second_file)
	colnames(meth) = c("probe", pats)


	pats = as.list(unique(colnames(meth)[2:ncol(meth)]))

	rearrange = function(pat){
		z = which(colnames(meth)==pat)
		newdat = meth[,c(1,z), with=FALSE]
		newdat = newdat[-1,]
		colnames(newdat) = c("probe", "beta", "gene", "chr", "coord")
		newdat$patient = pat
		return(newdat)
	}

	newdats = llply(pats, rearrange, .progress="text")
	dat <- ldply(newdats, data.table)
	dat$cancer = unlist(strsplit(first_file, split="_"))[1]

	saveRDS(dat, file=paste(substr(first_file, 1, 4), "all_methylation", ".rds", sep="_"))
	print(paste(substr(first_file, 1, 4), "complete"))
	return(paste(substr(first_file, 1, 4), "complete"))
}


#mapply(transform_files, files_first, files_second)

mcmapply(transform_files, files_first, files_second, mc.cores = 4L)


print("finished yay")






























