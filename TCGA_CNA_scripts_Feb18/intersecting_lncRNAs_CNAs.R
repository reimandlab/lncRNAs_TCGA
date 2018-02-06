###check if any of the CNAs from TCGA actually cover any of the lncRNAs - especially the candidates

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#Data

#1. CNAs
cnas = readRDS("OvaryLiverPancreasKIRC_CNA_TCGA_files_wTCGA_IDs_Feb6.rds")

#add cancer type to each one 
#TCGA source codes
tss_codes = read.csv("TCGA_TissueSourceSite_Codes2017.csv")

#Seperate into cancers 
get_cancer = function(tcga_id){
	tss = unlist(strsplit(tcga_id, "-"))[2]
	canc = as.character(tss_codes[which(tss_codes[,1] %in% tss), 3])
	if(length(canc)==0){canc = ""}
	return(canc)
}

cancers = unlist(llply(cnas$tcgaid, get_cancer, .progress="text"))
cnas$canc = cancers
z = which(cnas$canc == "")
cnas = cnas[-z,]

#2. GENCODE V27 lncRNAs GTF file
lncrnas = fread("gencode.v27.long_noncoding_RNAs.gtf")
coords = lncrnas[,1:8]
genes = lncrnas[,9]

genelist = as.list(genes$V9)

get_gene = function(gene){
	g = unlist(strsplit(gene, ";"))[1]
	gg = unlist(strsplit(g, " "))[[2]]
	gg = str_sub(gg, 2, 16)
	return(gg)
}

genenames = llply(genelist, get_gene)

genes$V9 = unlist(genenames)
coords = cbind(coords, genes)

#Process files 
cnas_bed = cnas[,c(2:8)]
cnas_bed[,1] = paste("chr", cnas_bed$Chromosome, sep ="")

coords_bed = coords[,c(1, 4:5, 9)]
colnames(coords_bed) = c("chr", "start", "end", "gene")

#3. FANTOM lncRNAs 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]
colnames(fantom)[1] = "gene"

z =  which(coords_bed$gene %in% fantom$gene)
coords_bed = coords_bed[z,]
coords_bed$name = ""
genelist = as.list(coords_bed$gene)
get_name = function(gene){
	name = fantom$CAT_geneName[which(fantom$gene == gene)]
	return(name)
}
genenames = llply(genelist, get_name)
coords_bed$name = unlist(genenames)

#candidate lncrnas 
cands = readRDS("36_unique_cands_4cancers_TCGA_Feb6.rds")
coords_bed = coords_bed[which(coords_bed$gene %in% cands$gene)]
coords_bed = coords_bed[!duplicated(coords_bed), ]

write.table(cnas_bed, file="cnas_bed_TCGA.bed", col.names=F, row.names=F, quote=F, sep="\t")
write.table(coords_bed, file="lncrna_coords_bed.bed", col.names=F, row.names=F, quote=F, sep="\t")

#################BEDTOOLS#########################################################################

bedtools intersect -a lncrna_coords_bed.bed -b cnas_bed_TCGA.bed -f 0.50 -wa -wb > fantom_lncrnas_wTCGA_CNAs_4cancers.bed





