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
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

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
ucsc = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]

#Methylation data 
meth = read.csv("HumanMethylation450probe_coordinatesIllumina.csv")
meth$Strand[meth$Strand=="F"] = "+"
meth$Strand[meth$Strand=="R"] = "-"
meth$CHR = paste("chr", meth$CHR, sep="")


meth$Strand[meth$Strand=="+"] = '+'
meth$Strand[meth$Strand=="-"] = '-'
meth$Strand[meth$Strand== ""] = '*'


#---------------------------------------------------------
#Make GRanges objects	
#---------------------------------------------------------

#1. Methylation CpG coordinates 
gr_meth = GRanges(seqnames = meth$CHR[1:485512], 
	ranges = IRanges(as.numeric(meth$MAPINFO[1:485512]), end =  as.numeric(meth$MAPINFO[1:485512])), 
	strand = meth$Strand[1:485512],
	gc = meth$Name[1:485512])

#2. lncRNAs 
gr_lncs = GRanges(seqnames = ucsc$hg19.ensGene.chrom, 
	ranges = IRanges(as.numeric(ucsc$hg19.ensGene.txStart), end = ucsc$hg19.ensGene.txStart), 
	strand = ucsc$hg19.ensGene.strand,
	gc = ucsc$hg19.ensemblToGeneName.value)

#Add 500 basepairs downstream 
down = flank(gr_lncs, 2500, start = FALSE)
#Add 500 basepairs upstream 
up = flank(gr_lncs, 2500)

down = unlist(down)
up = unlist(up)
gr_lncs = append(up, down)

#how many unique lncs - 6/8
lnc = subsetByOverlaps(gr_lncs, gr_meth)

#overlaps
lnc_probes = subsetByOverlaps(gr_meth, gr_lncs)

#Extract the methylation probe ids 
start = start(lnc_probes)
chr = as.character(seqnames(lnc_probes))

get = data.frame(chr = chr, start=start)

##keep probes
final_probes = subset(meth, ((meth$CHR %in% get$chr) & (meth$MAPINFO %in% get$start)))
write.csv(final_probes, file="39CandidatelncRNAprobesKIoct2017.csv", quote=F)


