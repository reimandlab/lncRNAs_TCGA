###Processing_clinical_data_ovarian_Cancer.R

#File downloaded from firehose 
#folder downlaoded called = gdac.broadinstitute.org_OV.Merge_Clinical.Level_1.2016012800.0.0
options(stringsAsFactors=F)

library(data.table)
library(stringr)
library(plyr)
library(dplyr)

f3 = fread("OV.clin.merged.txt", data.table=F, header=F)
f3 = t(f3)
cols = f3[1,1:ncol(f3)]
colnames(f3) = cols
f3 = f3[-1,]

cols = as.list(1:ncol(f3))

#remove columns that are missing data for 90% of patients 
checkna = function(col){
	z = which(is.na(f3[,col]))
	z2 = length(z)
	if(z2>=530){
	return(col)
	}
	else{return("no")}
}

remove = llply(cols, checkna)
remove = unlist(remove)
remove = remove[-c(which(remove=="no"))]
remove = as.numeric(remove)
f3 = f3[,-c(remove)]
f3 = as.data.frame(f3)

colss = colnames(f3)
which(str_detect(colss, "stage"))
which(str_detect(colss, "grade"))

#change TCGA IDs so it's all capitals 
f3$patient.bcr_patient_barcode = toupper(f3$patient.bcr_patient_barcode)

#591 unique lncRNAs downloaded and processed on Jan 23/18 by KI
saveRDS(f3, "591_OV_pats_clinical_data_Jan23_firehose.rds")