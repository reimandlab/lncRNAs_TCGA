###check if any of the CNAs from TCGA actually cover any of the lncRNAs - especially the candidates
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#Data
#1. CNAs
cnas = readRDS("OvaryLiverPancreasKIRC_CNA_TCGA_files_wTCGA_IDs_Feb7.rds")

#add cancer type to each one 
#TCGA source codes
tss_codes = read.csv("TCGA_TissueSourceSite_Codes2017.csv")
source_codes = source = read.csv("TCGA_sample_codes.csv")

cnas$source = ""
get_source = function(id){
	source = unlist(strsplit(as.character(id), '-'))[4]
	source = substr(source, 1, 2)
	return(source)
}
cnas$source = llply(cnas$Sample, get_source, .progress="text")
#only keep primary solid tumour samples 
cnas = dplyr::filter(cnas, source =="01")
clean_tcga_id = function(id){
	s1 = unlist(strsplit(id, '-'))[1]
	s2 = unlist(strsplit(id, '-'))[2]
	s3 = unlist(strsplit(id, '-'))[3]
	return(paste(s1, s2, s3, sep="-"))
}
cnas$patient = llply(cnas$Sample, clean_tcga_id, .progress="text")
cnas = as.data.table(cnas)
#Process files 
cnas_bed = cnas[,c(2:9)]
#keep only those with >10 probes and abs mean >0.2
cnas_bed = dplyr::filter(cnas_bed, Num_Probes>=10)
cnas_bed$patient = unlist(cnas_bed$patient)
cnas_bed$source = unlist(cnas_bed$source)

#2. GENCODE V27 lncRNAs GTF file
lncrnas = fread("gencode.v27lift37.long_noncoding_RNAs.gtf")
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
coords_bed$name = unlist(genenames) #5665 lncRNAs in file 

##keep just the start and end of the lncRNA, don't really need all transcripts 
lncs = as.list(unique(coords_bed$gene))
shorten = function(geneid){
	sub = dplyr::filter(coords_bed, gene == geneid)
	sub = sub[!duplicated(sub), ]
	if(dim(sub)[1] ==1){
		return(sub)
	}
	min = min(sub$start)
	max = max(sub$end)
	sub = sub[1,]
	sub$start = min
	sub$end = max 
	return(sub)
}
lncs_coords = llply(lncs, shorten, .progress="text")
lncs_coords <- ldply(lncs_coords, data.frame)

#candidate lncrnas 
cands = readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
lncs_coords = lncs_coords[which(lncs_coords$gene %in% cands$gene),] #21/25 are here in GENCODE 

write.table(cnas_bed, file="cnas_bed_TCGA.bed", col.names=F, row.names=F, quote=F, sep="\t")
write.table(lncs_coords, file="lncrna_coords_bed.bed", col.names=F, row.names=F, quote=F, sep="\t")

#################BEDTOOLS#########################################################################

bedtools intersect -a lncrna_coords_bed.bed -b cnas_bed_TCGA.bed -f 0.90 -wa -wb > fantom_lncrnas_wTCGA_CNAs_4cancers.bed





