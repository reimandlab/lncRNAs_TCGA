###cluster_OV_pats_by_CNAs.R

###Load Libraries
#++++++++++++++++++++++++++++++++++++++++
source("source_file.R")
library(GenomicRanges)
library(stringr)
#Load necessary packages
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

###Data
#++++++++++++++++++++++++++++++++++++++++
ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

##which lncRNAs are covered by probes 
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
lncs_ucsc = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]

#lncRNA candidate expression status 
medianexp = read.table("medianScoresOvarianCancerTop3_lncRNAs.txt", sep=";", header=T)

#copy number data for ovarian patients from JR's file obtained on October 30th 
cna = readRDS("ovarianPatients_CNA_data.rds") 

###CNA figure out how segments match to gene 
segments = rownames(cna) 
segments = str_split_fixed(segments, "::", n=2)

#make sure it is still atttached to patients so that later doesn't get lost, order of cnas...
cna = cbind(segments, cna)

#some rows are segments and some are gene IDs, let's divide those into seperate dataframes
segments1 = cna[which(cna[,1] == "ENCODEmerge"),]
segments2 =  cna[which(!(cna[,1] == "ENCODEmerge")),]

#segments 1
chr = str_split_fixed(segments1[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments1 = cbind(segments1, chr, coordinates)
segments1 = as.data.frame(segments1)
segments1$strand = "*"

#segments 2
chr = str_split_fixed(segments2[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments2 = cbind(segments2, chr, coordinates)
segments2 = as.data.frame(segments2)
segments2$strand = "*"

#seperate gene ID from the rest 
gene = str_split_fixed(segments2[,2], "::", n=3)
gene = cbind(segments2, gene)
gene = gene[which(gene[,73] %in% c("gencode")),]

g = function(row){
	gg = row[[75]]
	gg = gsub("\\..*","",gg)
	return(gg)
}
gene$ensg = apply(gene, 1, g)

#-----------------------------------------------------------------------------------------------
###Data
#++++++++++++++++++++++++++++++++++++++++

#Using copy number matrix, cluster patients based on H/L tag of lncRNA expression 

#remove duplicates
z <- which(duplicated(gene$ensg))
gene = gene[-z,]

rownames(gene) = gene$ensg
gene = gene[,-(c(1:2, 69:76))]

###add patient info - expression status-----------------------------------------------------
exp = c(1:66)
names(exp) = colnames(gene)[1:66]

for(i in 1:66){
	pat = colnames(gene)[i]
	z <- which(medianexp$patient %in% pat)
	if(!(length(z)==0)){
		exp[i] = medianexp$GS1[z]
	}
}

gene = rbind(gene, exp)
rownames(gene)[22873] = "GS1Expression"

#Remove all genes that are just NAs
rm = c()
for(i in 1:nrow(gene)){
	z <- which(is.na(gene[i,]))
	if(length(z) > 0){
		rm = c(rm, i)
	}
}

gene = gene[-rm, ]

#Remove all genes that have no copy number aberations
rm = c()
for(i in 1:nrow(gene)){
	z <- length(which(gene[i,] == 0))
	if(z > 40){
		rm = c(rm, i)
	}
}

length(rm)
gene = gene[-rm, ]

#make colour labels for the lncRNA expression status 
clab = as.data.frame(matrix(nrow=ncol(gene)))
colnames(clab)[1] = "lncRNA_Expression"
clab$lncRNA_Expression = as.character(gene[nrow(gene),])
clab$lncRNA_Expression[clab$lncRNA_Expression=="high"] = "darkred" 
clab$lncRNA_Expression[clab$lncRNA_Expression=="low"] = "darkorchid"
clab = as.matrix(clab)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create heatmap using custom heatmap.3 source code loaded above
par(cex.main=1)
cna_matrix = as.matrix(gene[1:21495,])
cna_matrix = apply(cna_matrix, 2, as.numeric)
rownames(cna_matrix) = rownames(gene)[1:21495]

#calculate variance for each gene
#keep top variable genes for now
vars = apply(cna_matrix, 1, var)
variance = as.data.frame(vars)
variance$gene = rownames(cna_matrix)
variance = as.data.table(variance)
variance = variance[order(-vars)]

top500 = variance[1:50,]
cna_matrix = cna_matrix[which(rownames(cna_matrix) %in% top500$gene),]

heatmap.3(cna_matrix, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
          Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", labCol=FALSE, labRow=FALSE, cexRow=1, col=rev(heat.colors(75)),
          ColSideColorsSize=7)



heatmap.2(cna_matrix, col=redgreen(200), ColSideColors=clab[,1],scale="row",
           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


#the above don't really work

#try counting the number of CNAs each gene has 
#postivie versus negative and compare this way 

cna_results = as.data.frame(matrix(ncol=2))
colnames(cna_results) = c("gene", "pval")

for(i in 1:(nrow(gene)-1)){
	pcg = gene[c(i, nrow(gene)), ]
	pcg = t(pcg)
	colnames(pcg) = c("CNA_gene", "Expression_GS1")
	pcg = as.data.frame(pcg)
	pcg$patient = rownames(pcg)
	pcg$CNA_gene = as.numeric(pcg$CNA_gene)
	w = wilcox.test(pcg$CNA[pcg$Expression_GS1=="high"],pcg$CNA[pcg$Expression_GS1=="low"], distribution = "exact")
	pval = w$p.value
	
	add = c(rownames(gene)[i], pval)
	names(add) = colnames(cna_results)
	cna_results = rbind(cna_results, add)
	
}


#plot the ones that are signficiant 




































