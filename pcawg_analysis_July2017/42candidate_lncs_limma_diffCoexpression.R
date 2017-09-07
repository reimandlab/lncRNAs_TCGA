#42candidate_lncs_limma_diffCoexpression .R

#Karina Isaev
#September 5th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run co-expression analysis with cancer type and sex? as confounders
#to get list of co-expressed PCGs with NEAT1 to run through m:Explorer  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(parallel)
library(limma)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data#-------------------------------------------------

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

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin_pats <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin_pats$icgc_donor_id %in% rownames(lnc_rna))
clin_pats <- clin_pats[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin_pats$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin_pats$icgc_donor_id),] #485 patients remain 

#make sure there aren't any lncRNAs in protein coding file 
table(ucsc[,7][ucsc[,8] %in% colnames(pcg_rna)])
#Remove 
z <- which(colnames(pcg_rna) %in% fantom[,2])
pcg_rna <- pcg_rna[,-z]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------

#List of canddidates and cox results
allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")

#[1] - for each lncRNA, divide expression data by patient high or low lncRNA expression
#for each row of allCands

#FUNCTION1 - divide expression into high and low for each lnc within a cancer type 
#Apply this function to every row of allCands
getExpression <- function(row){
	lnc <- row[[1]]
	canc <- row[[5]]
	#First divide patients into high and low based on median expression of lnc
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5608:5609)]
	lncData <- lncData[lncData$canc==canc,]
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

#make scatter plot showing patient vs lncRNA candidate expression
#what does the curve look like? how close are the points to each other?

simplePlot <- function(d){
	toplot <- as.data.table(d[,1:4])
	colnames(toplot)[2] <- "lnc"
	toplot <- toplot[order(lnc)]
	toplot$id <- rownames(toplot)
	toplot$exp[toplot$exp==0] <- "Low lncRNA"
	toplot$exp[toplot$exp==1] <- "High lncRNA"

	toplot$exp <- as.factor(toplot$exp)
	g <- ggscatter(toplot, x = "id", y = "lnc", color = "exp",
  				palette = mypal)
	g <- ggpar(g, xlab=FALSE,  legend = "right", legend.title = "lncRNA Group", ylab="lncRNA log1p(FPKM) Expression",
		main = paste(colnames(d)[2], "Expression Distribution"))
	g <- g + rremove("x.ticks")
	print(g + rremove("x.text"))
}

pdf("42lncRNAcandidates_distributionofExpSept6.pdf", pointsize=6)
llply(dividedWpcgs, simplePlot, .progress = "text")
dev.off()


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

    #rank list of genes before making heatmap
    t <- as.data.table(t)
    #first by adj p values then by decreasing FC
    t <- t[order(adj.P.Val)]
    t <- t[order(-abs(as.numeric(logFC)))]

    #save top gene names 
    top <- c(paste(colnames(d)[2], d[1,3]), t$ID)

    if(dim(t)[1] >= 200){
	#save 200 most DE genes 
    t <- t[1:200,]
	}

    #generate volcano plot
    print(ggplot(aes(x=logFC, y= -log10(P.Value)), data=t) + geom_point(aes(colour = -log10(adj.P.Val))) + scale_colour_gradient(low = "blue", high="red") +
    	labs(colour = "-log10(fdr)", x = "Limma logFC", y= "-log10(p-value)") + ggtitle(paste(colnames(d)[2], d[1,3]), "Differential Expression"))
    
    #generate heatmap 
    heat <- expression[which(rownames(expression) %in% t$ID),]
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

	#split into positive and nagative co-experessed genes 
	neg <- t[logFC <0]
	#pathway enrichment looking at top 100 coexpressed genes based on the above ranking
	negGenes <- neg[,1]
	negNoradpaths <- gprofiler(negGenes, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr")
	if(!(dim(negNoradpaths)[1]==0)){
	negNoradpaths <- negNoradpaths[, c(9, 12, 3)]
	colnames(negNoradpaths) <- c("GO.ID", "Description", "p.Val")
	negNoradpaths$FDR <- negNoradpaths$p.Val
	negNoradpaths$Phenotype <- "-1"
	}
	pos <- t[logFC >0]
	#pathway enrichment looking at top 100 coexpressed genes based on the above ranking
	posGenes <- pos[,1]
	posNoradpaths <- gprofiler(posGenes, organism = "hsapiens", ordered_query= TRUE, min_set_size=20, max_set_size = 300, min_isect_size=5, correction_method="fdr")
	if(!(dim(posNoradpaths)[1]==0)){
	posNoradpaths <- posNoradpaths[, c(9, 12, 3)]
	colnames(posNoradpaths) <- c("GO.ID", "Description", "p.Val")
	posNoradpaths$FDR <- posNoradpaths$p.Val
	posNoradpaths$Phenotype <- "+1"
	}

	#combine 
	allPaths <- rbind(negNoradpaths, posNoradpaths)

	if(!(dim(allPaths)[1]==0)){
	#only keep GO or REACTOME
	reac <- grep("REAC", allPaths$GO.ID)
	go <- grep("GO", allPaths$GO.ID)
	allPaths <- allPaths[c(reac, go), ]

	write.table(allPaths, sep= "\t", file=paste(colnames(d)[2], "PathwaysUsingtTop200DEgenesSept7.txt", sep="_"), quote=F, row.names=F)
	}
	
	dev.off()

	}

    if(dim(t)[1] <= 1){
    	top <- c(paste(colnames(d)[2], d[1,3]), "none")
   	}

    return(top)   	
}

diffExpressed <- llply(dividedWpcgs, diffE, .progress = "text")

saveRDS(diffExpressed, file="Sept7updated_42candidatesWithDEpcgsSept7TOP200genesUsed.rds")





































