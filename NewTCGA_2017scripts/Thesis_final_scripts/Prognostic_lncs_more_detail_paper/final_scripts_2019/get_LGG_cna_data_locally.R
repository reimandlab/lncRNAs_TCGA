library(TCGAbiolinks)
library(data.table)
query <- GDCquery(project = "TCGA-LGG",
             data.category = "Copy Number Variation",
             data.type = "Gene Level Copy Number Scores",              
             access="open")
GDCdownload(query)
data <- as.data.table(GDCprepare(query))
colnames(data)[4:ncol(data)] = sapply(colnames(data)[4:ncol(data)], function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
#data[,1] = unlist(data[,1])
data$gene_id = as.character(unlist(data[,1]))
data$gene_id = sapply(data$gene_id , function(x){unlist(strsplit(x, "\\."))[1]})

#tp53 - ENSG00000141510
z = which(data$gene_id == "ENSG00000141510")

saveRDS(data, file="LGG_copy_number_variation.rds")

maf <- GDCquery_Maf("LGG", pipelines = "muse")
maf=as.data.table(maf)
maf$Tumor_Sample_Barcode = sapply(maf$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
maf = as.data.table(subset(maf, Hugo_Symbol %in% c("IDH1", "IDH2", "TP53")))
maf = as.data.table(subset(maf, !(Variant_Classification %in% c("Intron", "Silent"))))
maf = maf[,1:25]

write.table(maf, "LGG_maf_idh_tp53.txt", quote=F, row.names=F)