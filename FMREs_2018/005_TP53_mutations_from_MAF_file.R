#f = read.table(gzfile("mc3.v0.2.8.PUBLIC.maf.gz"),header=T, sep="\t", stringsAsFactors=F)

library(data.table)

x=fread("gunzip -c mc3.v0.2.8.PUBLIC.maf.gz")

#first just get Tp53 mutations only 
z = which(x$Hugo_Symbol == "TP53")
tp53_muts = x[z,]

#save just in case
saveRDS(tp53_muts, file="TCGA_TP53_muts.rds")

#remove silent ones
z = which(tp53_muts$Variant_Classification == "Silent")
tp53_muts = tp53_muts[-z,]

#for now just start with misense mutations
z = which(tp53_muts$Variant_Classification == "Missense_Mutation")
tp53_muts = tp53_muts[z,] #2927 mutations

#change patient barcodes
change_col = function(pat){
	p1 = unlist(strsplit(pat, "-"))[1]
	p2 = unlist(strsplit(pat, "-"))[2]
	p3 = unlist(strsplit(pat, "-"))[3]
	newpat = paste(p1, p2, p3, sep="-")
	return(newpat)
}

library(plyr)
library(dplyr)
tp53_muts$Tumor_Sample_Barcode = llply(tp53_muts$Tumor_Sample_Barcode, change_col)
tp53_muts$Tumor_Sample_Barcode = unlist(tp53_muts$Tumor_Sample_Barcode)

#remove SNP ids
z = which(tp53_muts$dbSNP_RS == ".")
tp53_muts = tp53_muts[z,]
head(tp53_muts[,1:20])

saveRDS(tp53_muts, file="TCGA_TP53_muts.rds")
