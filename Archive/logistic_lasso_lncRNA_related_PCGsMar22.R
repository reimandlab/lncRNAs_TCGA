library(glmnet)
library(Hmisc)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

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
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
#pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#make sure there aren't any lncRNAs in protein coding file 
table(ucsc[,7][ucsc[,8] %in% colnames(pcg_rna)])
#Remove 
z <- which(colnames(pcg_rna) %in% fantom[,2])
#pcg_rna <- pcg_rna[,-z]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------

#List of canddidates and cox results
allCands <- readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
allCands$Cancer[allCands$Cancer =="lihc"] = "Liver hepatocellular carcinoma"
allCands$Cancer[allCands$Cancer =="ov"] = "Ovarian serous cystadenocarcinoma"
allCands$Cancer[allCands$Cancer =="kirc"] = "Kidney renal clear cell carcinoma"
allCands$Cancer[allCands$Cancer =="kirc"] = "Kidney renal clear cell carcinoma"
allCands$Cancer[allCands$Cancer =="paad"] = "Pancreatic adenocarcinoma"

z = which(duplicated(allCands$gene))

lnc_rna[,1:5919] = log1p(lnc_rna[,1:5919])

getExpression <- function(row){
	lnc <- row[[1]]
	canc <- row[[3]]
	#First divide patients into high and low based on median expression of lnc
	lncData <- lnc_rna[, c(which(colnames(lnc_rna) %in% lnc), 5920:5924)]
	lncData <- lncData[lncData$canc==canc,]
	lncData$exp = ""
	med = median(lncData[,1])
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(lncData[,1] > 0)
    l2 = which(lncData[,1] == 0)
    lncData$exp[l1] <- 1
	lncData$exp[l2] <- 0 
    }

    if(!(med ==0)){
    l1 = which(lncData[,1] >= med)
    l2 = which(lncData[,1] < med)
    lncData$exp[l1] <- 1
	lncData$exp[l2] <- 0 
    } 
	return(lncData)
}

divided <- apply(allCands, 1, getExpression)
print("pass")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  
getPCGS <- function(df){
	all <- merge(df, pcg_rna, by=c("patient",
		"canc", "time", "status", "sex"))
	#col 2 is the lncRNA expression data 
	#cols 1, 3, 4 are patient data 

	#Remove PCGs with median E < 5 FPKM 	
	#get medians of all PCGs
	meds <- apply(all[,8:ncol(all)], 2, median)

	#names of pcgs with median <5 
	low <- names(meds[which(meds <1500)]) #1st qu of median expression 
	all <- all[,-(which(colnames(all) %in% low))] 
	#log all gene expression values
	all[,c(6,8:(ncol(all)))] <- log1p(all[,c(6,8:(ncol(all)))])

	return(all)
}

dividedWpcgs <- llply(divided, getPCGS, .progress = "text")
print("pass2")

#make scatter plot showing patient vs lncRNA candidate expression
#what does the curve look like? how close are the points to each other?

	simplePlot <- function(d){
		toplot <- as.data.table(d[,1:7])
		colnames(toplot)[6] <- "lnc"
		toplot <- toplot[order(lnc)]
		toplot$id <- rownames(toplot)
		toplot$exp[toplot$exp==0] <- "Low lncRNA"
		toplot$exp[toplot$exp==1] <- "High lncRNA"

		toplot$exp <- as.factor(toplot$exp)
		g <- ggscatter(toplot, x = "id", y = "lnc", color = "exp",
	  				palette = mypal)
		g <- ggpar(g, xlab=FALSE,  legend = "right", legend.title = "lncRNA Group", ylab="lncRNA log1p(FPKM) Expression",
			main = paste(colnames(d)[6], d$canc[1], "Expression Distribution"))
		g <- g + rremove("x.ticks")
		print(g + rremove("x.text"))
	}

	pdf("25lncRNAcandidates_distributionofExp_Mar22.pdf", pointsize=6, width=10, height=7)
	llply(dividedWpcgs, simplePlot, .progress = "text")
	dev.off()


find_correlated_PCGS = function(d){
#1. Subset to the ones that are significantly correlated with lncRNA 
z = which(is.na(d$time))
if(!(length(z)==0)){
	d = d[-z,]
}
res2 = rcorr(as.matrix(d[,c(6, 8:ncol(d))]))
res2 = flattenCorrMatrix(res2$r, res2$P)
res2 = subset(res2, row == colnames(d)[6])
res2$fdr = p.adjust(res2$p, method="fdr")
res2 = subset(res2, fdr <=0.001)
res2 = subset(res2, abs(cor) >=0.3)

sort.field = "cor"
sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
res2 = as.data.table(res2)
res2 = sortme(res2, sort.field)
if(dim(res2)[1]>10){
genes = res2$column
d = d[,c(1:7, which(colnames(d) %in% genes))]
gene_data = d[,8:ncol(d)]
x <- model.matrix( ~., gene_data)

#y <- as.factor(d$exp)
#fit = glmnet(x, y, family = "binomial")
#cvfit = cv.glmnet(x, y, family = "binomial", alpha =1) #uses cross validation to select
y = as.numeric(d[,6])
cvfit = cv.glmnet(x, y, alpha =1)
#the best lamda and then use lambda to see which features remain in model 
cvfit$lambda.min #left vertical line
cvfit$lambda.1se #right vertical line 

#active covariates and their coefficients 
coef.min = coef(cvfit, s = "lambda.min") 
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

genes_keep = rownames(coef.min)[active.min]
print(genes_keep)
if(!(length(genes_keep)==0)){
d = d[,c(1:7, which(colnames(d) %in% genes_keep))]
res = cor(as.matrix(d[,c(6, 8:ncol(d))]))
library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
res2 = rcorr(as.matrix(d[,c(6, 8:ncol(d))]))
res2 = flattenCorrMatrix(res2$r, res2$P)
res2 = subset(res2, row == colnames(d)[6])

#give gprofiler two lists
list <- list()
l1 <- filter(res2, cor > 0)
list[[1]] <- l1$column #negative fold change, more expressed in the high lncRNA group 
l2 <- filter(res2, cor < 0)
list[[2]] <- l2$column #postivie fold change, more expressed in the low lncRNA group 
combined_paths <- gprofiler(list, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 300, min_isect_size=2, correction_method="fdr")
print(dim(combined_paths)[1])

if(!(dim(combined_paths)[1]==0)){
#only keep GO or REACTOME
reac <- grep("REAC", combined_paths$term.id)
go <- grep("GO", combined_paths$term.id)
combined_paths <- combined_paths[c(reac, go), ]
combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
combined_paths$Phenotype[combined_paths$Phenotype==1] = "1"
combined_paths$Phenotype[combined_paths$Phenotype==2] = "-1"
write.table(combined_paths, sep= "\t", file=paste(colnames(d)[6], d$canc[1], "PathwaysUsingtALL_DEgenesFeb16_linearLASSO.txt", sep="_"), quote=F, row.names=F)
}

#heatmap 
heat <- d[,8:ncol(d)]
rownames(heat) = d$patient
tags <- d$exp
color.map <- function(tags) { if (tags==1) "#FF0000" else "#0000FF" }
patientcolors <- unlist(lapply(tags, color.map))
	#heatmap.2(heat, col=topo.colors(200), scale="row", ColSideColors=patientcolors,
          #key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5, keysize=1.05)
	
	# cluster on correlation
heat = t(heat)
hc <- hclust(as.dist(1 - cor(t(heat))), method="ward.D2")
# draw a heatmap
my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 100)
heatmap.2(as.matrix(heat), main = paste(colnames(d)[6], d$canc[1]),col=my_palette, ColSideColors= patientcolors, cexRow=0.5, cexCol=0.6, Rowv=as.dendrogram(hc), trace="none", scale="row")
#heatmap.2(as.matrix(heat), main = paste(colnames(d)[6], d$canc[1]),col=my_palette, ColSideColors= patientcolors, cexRow=0.5, cexCol=0.6, Rowv=as.dendrogram(hc), trace="none")
}
}
}

for(i in 1:length(dividedWpcgs)){
	d = dividedWpcgs[[i]]
	pdf(paste(colnames(d)[6], d$canc[1], "linear_lasso_PCG_analysis_Feb16.pdf", sep="_"))
	find_correlated_PCGS(d)
	dev.off()
	print(i)
}





















