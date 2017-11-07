###Processing_450methylationfile.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

files = list.files(pattern= "x") #pattern of mini files 
files = files[2:11517]

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")

###Data
#+++++++++++++++++++++++++++++++++

canc_pats = fread("485_patient_IDs_top5CancersPCAWG.txt")
canc_pats = filter(canc_pats, canc == "Liver Hepatocellular carcinoma")


lnc_probes = read.table("cand_lincs_wMethylationProbes.txt")
colnames(lnc_probes) = c("Chr_lnc", "start_lnc", "end_lnc", "ensg", "hugo", "type", 
	"transcript", "chr_probe", "start_probe", "end_probe", "probe")

for(i in 1:length(files)){

	f = fread(files[i], data.table=F)
	f = subset(f, f[,1] %in% canc_pats$patient)
	if(!(nrow(f) == 0)){
	f = subset(f, f[,8] %in% lnc_probes$probe)
	if(!(nrow(f) == 0)){
	name = paste("newfiles/", i, ".txt", sep="_")
	write.table(f, name, quote=F)
}
}
}

print("end")


