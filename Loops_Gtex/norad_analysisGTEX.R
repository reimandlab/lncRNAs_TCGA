#norad_analysisGTEX.R

###GTEX_42_candidatelncs_limma_diffCoexpression.R

#Author: Karina_Isaev
#Date_started: September 6ths 2017
#dir: Thesis/GTEx_data/data/raw

#Description:
#Processing RNA-Seq file from GTEx
#Reformating loops file
#Note: these values are RPKM 
#There are 4,859 unique tissue samples with RNA-Seq data

###Preamble###############################################
options(stringsAsFactors=F)

###Libraries##############################################
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
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(limma)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Multiplot function
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

###Data#####################################################

#RNA-Seq file 
rna <- readRDS("gtex_expression_file.rds")
rna <- as.data.frame(rna)

#Clinical file 
clin <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt") ; clin <- as.data.frame(clin)

#------------------------------------------------------------------
#Processing clinical file - inlcude samples that have expression
#------------------------------------------------------------------

tissues <- as.data.frame(clin[,c(1,6)])
c <- as.data.frame(table(tissues$SMTS))
c[,1] <- as.character(c[,1])
c[,2] <- as.numeric(c[,2])

#PCAWG Cancers top 5 
cancers <- c("Breast" , "Kidney" ,  "Liver" , "Ovary" , "Pancreas")	

z <- which(c[,1] %in% cancers)
cancers_keep <- c[z,1]

#Filter clinical file
z <- which(clin[,6] %in% cancers_keep)
clin <- clin[z,]

#keep only patients in gene expression file that are part of the tissues
#we want to study 
z <- which(colnames(rna) %in% clin[,1])
rna <- rna[,c(1,2,z)]

#keep only patients in clinical file that have gene-Expression 
z <- which(clin[,1] %in% colnames(rna))
clin <- clin[z,]

####-------------------
### 633 samples TOTAL
####-------------------

##PCAWG top 5 high cancers 
high_lncs <- fread("high_lncsmed4top5cancersPCAWG.txt", sep=";")

##PCAWG sig lncRNAs < 0.05 pvalue from high lncs
sig <- fread("42sig_lncRNACancerAssociations.txt", sep=";")
sig$canc <- lapply(sig$canc, function(x) unlist(strsplit(x, " "))[1])
fdr_sig <- sig[fdr <0.1]

allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")

#list of functional lncRNAs from FANTOM5 paper (will assume that all other genes othan than these are protein-coding for now)
#can always subset the protein coding gene expression file later 
lncs <- fread("lncs_ensg_fromFantom5paper_downjune7th")

#ucsc genes
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

###Processing#################################################head()

#1. Want to only look at ENSG genes in rna file
#split feature column and extract third component, write function and apply

#seperate first by "."
extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "\\..*"))
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 
lncs[,1] <- apply(lncs[,1], 1, extract) ; 

#2. remove duplicates 
rna <- rna[! rna$Description %in% unique(rna[duplicated(rna$Description), "Description"]), ]

#ALL GENES USED IN PCAWG
genes <- fread("all_genes_used_inRankingAnalysis.txt", sep=";")

z <- which(rna[,2] %in% genes$x)
rna <- rna[z,] #25125 both pcgs and lncrnas 

rownames(rna) <- rna[,2]
rna <- rna[,-c(1,2)]
rna <- t(rna)
rna <- as.data.frame(rna)
rna$tissue <- ""
rna$patient <- ""
rna$patient <- rownames(rna)

for(i in 1:nrow(rna)){
	z <- which(clin$SAMPID %in% rna$patient[i])
	rna$tissue[i] <- clin[z,6]
}

z <- which(rna$tissue %in% c(""))

#need to seperate into lncrnas and pcgs 
lncsPcawgList <- fread("lncsPcawgList.txt", sep=";")
pcgsPcawgList <- fread("pcgsPcawgList.txt", sep=";")

lnc_rna <- rna[,c(which(colnames(rna) %in% lncsPcawgList$x), 25126)]
pcg_rna <- rna[,c(which(colnames(rna) %in% pcgsPcawgList$x), 25126)]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------

#List of canddidates and cox results
allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
allCands$canc <- unlist(lapply(allCands$canc, function(x) unlist(strsplit(x, ' '))[1]))

#[1] - for each lncRNA, divide expression data by patient high or low lncRNA expression
#for each row of allCands

#FUNCTION1 - divide expression into high and low for each lnc within a cancer type 
#Apply this function to every row of allCands
getExpression <- function(row){
	lnc <- row[[1]]
	canc <- row[[5]]
	#First divide patients into high and low based on median expression of lnc
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5488:5489)]
	lncData <- lncData[lncData$tissue==canc,]
	med <- median(as.numeric(lncData[,1]))
	lncData$exp[lncData[,1] >= med] <- 1
	lncData$exp[lncData[,1] < med] <- 0 
	return(lncData)
}

divided <- apply(allCands, 1, getExpression)
print("pass")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  
getPCGS <- function(df){
	all <- merge(df, pcg_rna, by="patient")
	all <- all[,-(ncol(all))]
	colnames(all)[3] <- "canc"
	#col 2 is the lncRNA expression data 
	#cols 1, 3, 4 are patient data 

	#Remove PCGs with median E < 5 FPKM 	
	#get medians of all PCGs
	meds <- apply(all[,5:ncol(all)], 2, median)

	#names of pcgs with median <5 
	low <- names(meds[which(meds <5)]) 
	all <- all[,-(which(colnames(all) %in% low))] 
	#log all gene expression values
	all[,c(2,5:(ncol(all)))] <- log1p(all[,c(2,5:(ncol(all)))])

	return(all)
}

dividedWpcgs <- llply(divided, getPCGS, .progress = "text")
print("pass2")

#FUNCTION3 - limma differential expression between high and low lncRNA groups
diffE <- function(d){
	design <- model.matrix(~ 0 + factor(d$exp))
	colnames(design) <- c("Low", "High")
	rownames(d) <- d$patient

	expression <- t(d[,5:ncol(d)])
	fit <- lmFit(expression, design)
	cont.matrix <- makeContrasts(LowvsHigh=Low-High, levels=design)

	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)

	ps <- fit2$p.value
	ps <- p.adjust(ps, method="fdr")
	numGenes <- length(which(ps <= 0.05))

	par(mfrow=c(2,1))
	
	pdf(paste(colnames(d)[2], d[1,3], "volcano.pdf", sep="_"), pointsize=8, height=9, width=10)

	for (i in 1:ncol(fit2$p.value)) {
	print(hist(p.adjust(fit2$p.value[,i], method="fdr"), main=colnames(fit2$p.value)[i]))
		}

	genes=rownames(expression)
    t <- topTable(fit2,coef=1,adjust.method="fdr",n=numGenes,p.value=0.05,genelist=genes)

    if(dim(t)[1] > 1){

    #save top 100 gene names 
    top <- c(paste(colnames(d)[2], d[1,3]), t$ID)

    #generate volcano plot
    print(ggplot(aes(x=logFC, y= -log10(P.Value)), data=t) + geom_point(aes(colour = -log10(adj.P.Val))) + scale_colour_gradient(low = "blue", high="red") +
    	labs(colour = "-log10(fdr)", x = "Limma logFC", y= "-log10(p-value)") + ggtitle(paste(colnames(d)[2], d[1,3]), "Differential Expression"))
    
    #generate heatmap 
    heat <- expression[which(rownames(expression) %in% top),]
    tags <- d$exp
	color.map <- function(tags) { if (tags==1) "#FF0000" else "#0000FF" }
    patientcolors <- unlist(lapply(tags, color.map))
	#heatmap.2(heat, col=topo.colors(200), scale="row", ColSideColors=patientcolors,
          #key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5, keysize=1.05)
	
	# cluster on correlation
	hc <- hclust(as.dist(1 - cor(t(heat))), method="ward.D2")
	# draw a heatmap
	heatmap(as.matrix(heat),
	Rowv=as.dendrogram(hc),
	col=greenred(100),ColSideColors= patientcolors, cexRow=0.5, cexCol=0.6)

	#pathway enrichment

	#split into positive vs negative 
	t <- as.data.table(t)
	pos <- t[logFC >0]
	neg <- t[logFC <0]

	posPaths <- gprofiler(pos$ID, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr", include_graph=T)
	negPaths <- gprofiler(neg$ID, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr", include_graph=T)

	#GMT file 
	gmt <- fread("pathways_names.gmt", data.table=F, header=F)
	colnames(gmt) <- c("Pathway", "Name")

	#Output needs to look like
	#GO.ID      {tab} Description                     {tab} p.Val {tab} FDR  {tab} Phenotype
	#GO:0000346 {tab} transcription export complex    {tab} 0.01  {tab} 0.02 {tab} +1
	#GO:0030904 {tab} retromer complex                {tab} 0.05  {tab} 0.10 {tab} +1
	#GO:0008623 {tab} chromatin accessibility complex {tab} 0.05  {tab} 0.12 {tab} -1
	#GO:0046540 {tab} tri-snRNP complex               {tab} 0.01  {tab} 0.03 {tab} -1

	pathsPosnew <- posPaths[, c(9, 12, 3)]
	colnames(pathsPosnew) <- c("GO.ID", "Description", "p.Val")
	#only keep GO or REACTOME
	reac <- grep("REAC", pathsPosnew$GO.ID)
	go <- grep("GO", pathsPosnew$GO.ID)
	pathsPosnew <- pathsPosnew[c(reac, go), ]

	pathsNegnew <- negPaths[, c(9, 12, 3)]
	colnames(pathsNegnew) <- c("GO.ID", "Description", "p.Val")
	#only keep GO or REACTOME
	reac <- grep("REAC", pathsNegnew$GO.ID)
	go <- grep("GO", pathsNegnew$GO.ID)
	pathsNegnew <- pathsNegnew[c(reac, go), ]

	write.table(pathsPosnew, sep= "\t", file=paste(colnames(d)[2], d[1,3], "positivegProfiler.txt", sep="_"), quote=F, row.names=F)
	write.table(pathsNegnew, sep= "\t", file=paste(colnames(d)[2], d[1,3],"negativegProfiler.txt", sep="_"),  quote=F, row.names=F)

	dev.off()

	}

    if(dim(t)[1] <= 1){
    	top <- c(paste(colnames(d)[2], d[1,3]), "none")
   	}

    return(top)   	
}

diffExpressed <- llply(dividedWpcgs, diffE, .progress = "text")

saveRDS(diffExpressed, file="42candidatesWithDEpcgsSept5.rds")










