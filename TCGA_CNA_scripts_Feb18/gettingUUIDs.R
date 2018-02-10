require(httr)
library(dplyr)
library(plyr)

data = readRDS("OvaryLiverPancreasKIRC_CNA_TCGA_files.rds")
ids = unique(data$patient)

uuid_barcodes = as.data.frame(matrix(ncol=2))
colnames(uuid_barcodes) = c("uuid", "barcode")

for(i in 1:length(ids)){

	files_endpt = "https://gdc-api.nci.nih.gov/files"

	id1 = ids[i]
	body = list("filters" = list("op" = "in", "content" = list('field'='files.file_id', 'value'=id1)), 
		'format' = 'TXT', 'fields'="file_id,file_name,cases.case_id,cases.submitter_id,cases.samples.sample_id", 'size' = 1)

	r <- POST(url = files_endpt, body = body, encode = 'json')

	content(r, 'parsed')

	barcode = content(r)$data$hits[[1]]$cases[[1]]$submitter_id
	row = c(id1, barcode)
	names(row) = colnames(uuid_barcodes)
	uuid_barcodes = rbind(uuid_barcodes, row)

}


uuid_barcodes = uuid_barcodes[-1,]

change = function(colname){
	z <- which(uuid_barcodes$uuid == colname)
	return(uuid_barcodes$barcode[z])	
}

data2 = data
#change to TCGA IDs
change_id = function(id){
	newid = uuid_barcodes$barcode[which(uuid_barcodes$uuid %in% id)]
	return(newid)
}
uuids = as.list(data2$patient)
newids = llply(uuids, change_id, .progress="text") 

data2$tcgaid = unlist(newids)

saveRDS(data2, file="OvaryLiverPancreasKIRC_CNA_TCGA_files_wTCGA_IDs_Feb6.rds")

