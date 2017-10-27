###Processing_450methylationfile.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

files = list.files(pattern= "pro") #pattern of mini files 

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")

###Data
#+++++++++++++++++++++++++++++++++

ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

lnc_probes = read.csv("5000lncRNAbesKIoct2017.csv")

for(i in 1:length(files)){

	f = fread(files[i])
	f = subset(f, f$icgc_donor_id %in% ov_pats$V2)
	f = subset(f, f$probe_id %in% lnc_probes$Name)
	name = paste("newfiles/", i, ".txt", sep="_")
	write.table(f, name, quote=F)
}

print("end")


