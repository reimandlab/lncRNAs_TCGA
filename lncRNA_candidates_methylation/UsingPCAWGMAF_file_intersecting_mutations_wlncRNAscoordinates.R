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

#mutations from MAF file 
maf = fread("merged_mutations_ovarian_cancer_patients.txt", fill=TRUE)
z = which(maf$V43 == "V43")
maf = maf[-z,]

colnames(maf)[15] = "patient"
maf = filter(maf, patient %in% ov_pats$V2)

#which patients US which patients AUS
us = fread("sample.OV-AU.tsv")
us = unique(us$icgc_donor_id)
aus = fread("sample.OV-US.tsv")
aus = unique(aus$icgc_donor_id)

maf$country = ""
maf$country[which(maf$patient %in% us)] ="US"
maf$country[which(maf$patient %in% aus)] = "AUS"
maf = maf[,-1]

colnames(maf) = c("Hugo", "Chr", "Start", "End", "Strand", "Variant_Classification", "Variant_Type", "Reference_allele", "Reference_allele2", "alt_allele", "dbSNP_RS", "gc_content",
	"Project_Code", "Donor_ID", "country")

#---------------------------------------------------------
#Make Bed files 
#---------------------------------------------------------
maf$Chr = paste("chr", maf$Chr, sep="")
maf = maf[,c(2:4, 1, 5:15)]
#maf = maf[,1:3]
write.table(maf, file= "maf_ovarian.bed" , sep = "\t", quote=F, row.names=F, col.names=F)

ucsc = ucsc[,c(2,4,5,1,3,6:8)]
#ucsc = ucsc[,1:3]
write.table(ucsc, file="3ovarianLNCs.bed",sep = "\t", quote=F, row.names=F, col.names=F)
