###compilingMiniMethylationFiles.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)

###Data
#+++++++++++++++++++++++++++++++++

ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

lnc_probes = read.csv("39CandidatelncRNAprobesKIoct2017.csv")

###Cat methylation results into one file 
results = fread("mergedMethylationsOvary.txt", fill=TRUE)
z <- which(results$V2 == "V2")
results = results[-z,]

#z <- which(results$V2 == "project_code")
#results = results[-z,]

###which lncRNAs are covered by probes 
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

#get coordinates of which probes overlap which lncRNA genes 
hits <- findOverlaps(gr_lncs, gr_meth)

idx <- (subjectHits(hits))
id2 <- (queryHits(hits))

values <- DataFrame(probes=gr_meth$gc[idx], lncs = gr_lncs$gc[id2])

#Merge results file with vlaues file, adding lncRNA ID to probe ID 
colnames(results)[9] = "probes"
results = merge(results, values, by = "probes")

#48 ovarian patients out of the 70 with methylation data for 6/8 lncRNAs 
write.csv(results, file="48ovaryPatients_methylation_6candidatelncRNAs.csv", quote=F)



