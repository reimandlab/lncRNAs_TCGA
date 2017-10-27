###Processing_450methylationfile.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

files = list.files(pattern= "pro") #pattern of mini files 
files = files[3:5624]

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")

###Data
#+++++++++++++++++++++++++++++++++

ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

lnc_probes = read.csv("39CandidatelncRNAprobesKIoct2017.csv")

for(i in 1:length(files)){

	f = fread(files[i])
	f = subset(f, f$V1 %in% ov_pats$V2)
	if(!(nrow(f) == 0)){
	f = subset(f, f$V8 %in% lnc_probes$Name)
	if(!(nrow(f) == 0)){
	name = paste("newfiles/", i, ".txt", sep="_")
	write.table(f, name, quote=F)
}
}
}

print("end")


