###methylationOfCandidatesOvary2lncRNAs.R
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

#ucsc only keep those lncrnas in fantom
ucsc = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]

#lncRNA candidate expression status 
medianexp = fread("medianScoresOvarianCancerTopTwolncRNAs.txt")
medianexp = medianexp[,-c(1, 5, 6)]
colnames(medianexp) = c("LINC00665","ZNF503","canc","patient","status","time", "sex")

#Methylation data
meth = read.csv("48ovaryPatients_methylation_6candidatelncRNAs.csv")

medianexp = filter(medianexp, patient %in% meth$V2)

###Subset to methylation data of two lncs
#++++++++++++++++++++++++++++++++++++++++

meth = subset(meth, meth$lncs %in% c("LINC00665", "ZNF503-AS2"))
meth$lncs[meth$lncs=="ZNF503-AS2"] = "ZNF503"

keep = c("probes", "V2", "V10", "V11", "lncs")
meth = meth[keep]
colnames(meth) = c("probe", "patient", "value", "metric", "lncRNA")

##Add tag for high/low expression 
meth = merge(meth, medianexp, by = "patient")

#plot 

linc00665 = meth[meth$lncRNA == "LINC00665",]
p1 = ggboxplot(linc00665, x = "LINC00665", y = "value", add = "jitter", fill="LINC00665")
p1 = ggpar(p1, palette = "Dark2", xlab = "LINC00665 Expression", ylab="Methylation Beta Value (All available probes)", legend="none")
#  Add p-value
p1 = p1 + stat_compare_means()

ZNF503 = meth[meth$lncRNA == "ZNF503",]
p2 = ggboxplot(ZNF503, x = "ZNF503", y = "value", add = "jitter", fill="ZNF503")
p2 = ggpar(p2, palette = "Dark2", xlab = "ZNF503-AS2 Expression", ylab="Methylation Beta Value (All available probes)", legend="none")
p2 = facet(p2, facet.by = "probe")

#  Add p-value
p2 = p2 + stat_compare_means()
pdf("ZNF503AS2_probes_methylationDiffsOvary.pdf", width=8, height=9)
print(p2)
dev.off()