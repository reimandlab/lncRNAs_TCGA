###intersecting_mutations_wlncRNAscoordinates.R
#++++++++++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#++++++++++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)

###Data
#++++++++++++++++++++++++++++++++++++++++

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
lncs = filter(lncs, canc == "Ovar")

#ucsc only keep those lncrnas in fantom
ucsc = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]

#lncRNA candidate expression status 
medianexp = fread("medianScoresOvarianCancerTop3_lncRNAs.txt")
medianexp = medianexp[,-c(1, 6, 7)]
colnames(medianexp) = c("GS1", "LINC00665","ZNF503","canc","patient","status","time", "sex")

#Mutation data
muts = fread("simple_somatic_mutation.open.OV-US.tsv")

#subset to patients that are in these 70 
muts = filter(muts, icgc_donor_id %in% ov_pats$V2)

#---------------------------------------------------------
#Make GRanges objects	
#---------------------------------------------------------

#1. Methylation CpG coordinates 
muts$chromosome = paste("chr", muts$chromosome, sep="")

gr_muts = GRanges(seqnames = muts$chromosome, 
	ranges = IRanges(as.numeric(muts$chromosome_start), end =  as.numeric(muts$chromosome_end)), 
	strand = muts$chromosome_strand,
	gc = muts$icgc_mutation_id)

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
comb = append(up, down)
gr_lncs = append(gr_lncs, comb)

#how many unique lncs - 6/8
lnc = subsetByOverlaps(gr_lncs, gr_muts, ignore.strand = TRUE)

