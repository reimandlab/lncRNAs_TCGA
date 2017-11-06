#make_bed_files.R

##methylationLNCRNAprobes.R

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(plyr)
library(GenomicRanges)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
#z <- which(duplicated(ucsc[,6]))
#ucsc <- ucsc[-z,]

#fantom 
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

#Candidate lncRNAs 
lncs = fread("results_October12_42candsFromPCAWG.txt")
lncs = filter(lncs, pval < 0.05)

#ucsc only keep those lncrnas in fantom
lincs = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]

#maf file should contain whatever mutation data there is for 485 patients 
maf = read.table("merged_maf_file_containing_485cancerpatients.maf", sep="}")
z <- which(maf$V14 %in% "")
maf = maf[-z,]
cols = c(3, 4, 5, 15, 14)
maf = maf[,cols]
maf$V3 = paste("chr", maf$V3, sep="")

write.table(maf, file="maf_mutation_data_485patientsTOP5cancers.bed", sep="\t", quote=F, row.names=F, col.names=F)

#Methylation data (ovarian only)
meth = read.csv("HumanMethylation450probe_coordinatesIllumina.csv")
meth$Strand[meth$Strand=="F"] = "+"
meth$Strand[meth$Strand=="R"] = "-"
meth$CHR = paste("chr", meth$CHR, sep="")

meth$Strand[meth$Strand=="+"] = '+'
meth$Strand[meth$Strand=="-"] = '-'
meth$Strand[meth$Strand== ""] = '*'

#meth bed file 
cols = c(12, 13, 13, 2)
meth_bed = meth[,cols
#remove NAs 
z <- which(is.na(meth_bed$MAPINFO))
meth_bed = meth_bed[-z,]
write.table(meth_bed, file="methylation_array_probe_info.bed", sep="\t", quote=F, row.names=F, col.names=F)


#lncRNA bed file 
cols = c(2, 4, 5, 6, 8, 7, 1)
lincs = lincs[,cols]
write.table(lincs, file="8lncRNAcandidates_coordiantes.bed", sep="\t", quote=F, row.names=F, col.names=F)
