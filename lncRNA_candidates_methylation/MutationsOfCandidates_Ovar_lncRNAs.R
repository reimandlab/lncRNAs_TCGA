###MutationsOfCandidates_Ovar_lncRNAs.R
#++++++++++++++++++++++++++++++++++++++++

#Purpose: this script will find how many mutations
#each lncRNA has for each patient 
#then can compare difference between high and low expressing lncRNA group

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
muts = fread("simple_somatic_mutation.open.OV-AU.tsv")

#subset to patients that are in these 70 
muts = filter(muts, icgc_donor_id %in% ov_pats$V2)

lnc_muts = muts[(which(muts$gene_affected %in% ucsc$hg19.ensGene.name2)),]

unique(lnc_muts$gene_affected)

#quickly check if genes in here:gene_affected
colnames(lnc_muts)[2] = "patient"
lnc_muts = merge(lnc_muts, medianexp, by = "patient")


meth = read.csv("48ovaryPatients_methylation_6candidatelncRNAs.csv")
pats = unique(meth$V2)


#LINC00665 
linc00665 = lnc_muts[lnc_muts$gene_affected == "ENSG00000232677",]
patswmutations = unique(linc00665$patient)
#three possible outcomes, no info, have mutation, don't have mutation 
medianexp$linc00665mutation = ""
medianexp$linc00665mutation[medianexp$patient %in% patswmutations] = 1
medianexp$linc00665mutation[!(medianexp$patient %in% patswmutations)] = 0
medianexp$linc00665mutation[!(medianexp$patient %in% pats)] = "na"

#ZNF503
ZNF503 = lnc_muts[lnc_muts$gene_affected == "ENSG00000237149",]
patswmutations = unique(ZNF503$patient)
#three possible outcomes, no info, have mutation, don't have mutation 
medianexp$znf503mutation = ""
medianexp$znf503mutation[medianexp$patient %in% patswmutations] = 1
medianexp$znf503mutation[!(medianexp$patient %in% patswmutations)] = 0
medianexp$znf503mutation[!(medianexp$patient %in% pats)] = "na"


#GS1
GS1 = lnc_muts[lnc_muts$gene_affected == "ENSG00000253738",]
patswmutations = unique(GS1$patient)
#three possible outcomes, no info, have mutation, don't have mutation 
medianexp$GS1mutation = ""
medianexp$GS1mutation[medianexp$patient %in% patswmutations] = 1
medianexp$GS1mutation[!(medianexp$patient %in% patswmutations)] = 0
medianexp$GS1mutation[!(medianexp$patient %in% pats)] = "na"


##fisher's test 



##heatmap 
heat = medianexp[,c(1,2,3,5,9,10,11)]
heat = as.data.frame(heat)
rownames(heat) = heat$patient
heat$patient = NULL

heat$GS1[heat$GS1 == "high"] = 1
heat$GS1[heat$GS1 == "low"] = 0

heat$LINC00665[heat$LINC00665 == "high"] = 1
heat$LINC00665[heat$LINC00665 == "low"] = 0

heat$ZNF503[heat$ZNF503 == "high"] = 1
heat$ZNF503[heat$ZNF503 == "low"] = 0

heat = heat[rownames(heat) %in% pats,]
distance = dist(heat, method = "euclidian")
cluster = hclust(distance, method = "ward.D2")

heat[,1:6] = as.numeric(heat[,1:6])

heatmap.2(heat,
  cellnote = heat,  # same data set for cell labels
  main = "Ovarian Cancer Data", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  Rowv = as.dendrogram(cluster), # apply default clustering method
  Colv = as.dendrogram(cluster)) # apply default clustering method














