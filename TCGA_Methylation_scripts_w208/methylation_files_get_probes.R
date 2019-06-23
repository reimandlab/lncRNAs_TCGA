###methylation files 
###identify probes overlapping lncRNAs 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#Data
#1. Probes 450 
#probe_coordinates = fread("hm450.hg19.manifest")
#strand is 6th column 
#probe_coordinates = probe_coordinates[,c(1:3, 5,9,4)]
#write.table(probe_coordinates, file="450meth_probe_coordinates.bed", quote=F, row.names=F, col.names=F, sep="\t")

#new probe file --> hg38 --> 450meth_probes_liftover_hglft_genome_90e_e56a20.bed

#2. GENCODE V27 lncRNAs GTF file
#lncrnas = fread("gencode.v27lift37.long_noncoding_RNAs.gtf")
#gencode.v30.long_noncoding_RNAs.gtf
lncrnas = fread("gencode.v30.long_noncoding_RNAs.gtf")

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
coords_bed = coords[,c(1, 4:5, 8,9,7)]
colnames(coords_bed) = c("chr", "start", "end", "rand", "gene", "strand")

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
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
colnames(cands)[7] = "canc"
cands = filter(cands, data == "TCGA")

lncs_coords = lncs_coords[which(lncs_coords$gene %in% cands$gene),] #154/166 are here in GENCODE 
write.table(lncs_coords, file="lncrna_coords_bed.bed", quote=F, row.names=F, col.names=F, sep="\t")

#################BEDTOOLS#########################################################################

bedtools intersect -a 450meth_probes_liftover_hglft_genome_90e_e56a20.bed -b lncrna_coords_bed.bed -wa -wb -s > fantom_lncrnas_mapped_to450_probes.bed

#3. OV Methylation 
#ov_meth = fread("OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")