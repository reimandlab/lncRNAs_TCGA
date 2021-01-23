library(data.table)
library(plyr)
library(dplyr)

setwd("/u/kisaev/TCGA_SKCM_count_data")

manifest= "gdc_manifest_20210122_222008.txt" #Manifest name
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))

Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '

Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'

Part3= paste0("\"",manifest_length, "\"", "}") #just change single quote ' to double quote "

Sentence= paste(Part1,id,Part2,Part3, collapse=" ")

#write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)

#in terminal
curl --request POST --header "Content-Type: application/json" --data @Payload.txt "https://api.gdc.cancer.gov/files" > File_metadata.txt

#change UUIDs --> patient IDs
