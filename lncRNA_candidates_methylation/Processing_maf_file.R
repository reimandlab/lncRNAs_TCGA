###Processing_maf_file.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

files = list.files(pattern= "output_file") #pattern of mini files 

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")

###Data
#+++++++++++++++++++++++++++++++++

ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

for(i in 1:length(files)){

	f = fread(files[i], sep ="\t", data.table=F)
	cols_keep = c(1:11, 17, 42:43)
	f = f[,cols_keep]
	f = subset(f, f[,14] %in% ov_pats$V2)
	if(!(nrow(f) == 0)){
	name = paste("newfiles/", i, ".txt", sep="_")
	write.table(f, name, quote=F, sep = ":")
}
}

print("end")


#gc_content 
