#---------------------------------------------------------------------------
#54_lncRNAs_direct_comparison_PCAWG_PanCanAtlas.R
#---------------------------------------------------------------------------

#Data: August 9th

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(colorout)
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
library(factoextra)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------------------------
#Data
#---------------------------------------------------------------------------

#PanCanAtlas Expression Data-Set 
pan <- readRDS("panAtlas_65lncs_RNA-Seq_data.rds")

#PCAWG Expression Data-Set 
pcawg <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")

#---------------------------------------------------------------------------
#Processing - subset to lncRNAs in common, n = 54
#---------------------------------------------------------------------------

z <- which(colnames(pcawg) %in% colnames(pan))
pcawg <- pcawg[,z]

z <- which(colnames(pan) %in% colnames(pcawg))
pan <- pan[,z]
#change cancer names so they match cancer names in pcawg
pcawg$canc[pcawg$can== "Kidney Adenocarcinoma, papillary type"] <- "Kidney renal papillary cell carcinoma"
pcawg$canc[pcawg$can== "Kidney Adenocarcinoma, clear cell type" ] <- "Kidney renal clear cell carcinoma" 
pcawg$canc[pcawg$can== "Breast Infiltrating duct carcinoma"] <-  "Breast invasive carcinoma"
pcawg$canc[pcawg$can== "Kidney Adenocarcinoma, chromophobe type" ] <- "Kidney Chromophobe"
pcawg$canc[pcawg$can== "Ovary Serous cystadenocarcinoma"] <- "Ovarian serous cystadenocarcinoma"  
pcawg$canc[pcawg$can== "Pancreas Pancreatic ductal carcinoma" ] <- "Pancreatic adenocarcinoma"    
pcawg$canc[pcawg$can== "Liver Hepatocellular carcinoma" ] <-  "Liver hepatocellular carcinoma"  

#---------------------------------------------------------------------------
#Plot each lncRNA's logged expression distribution side by side
#comparing pcawg and pancanatlas 
#---------------------------------------------------------------------------

#54 genes 
pdf("Expression_of_54_lncRNAs_sequencedinBOTHdatesets.pdf", pointsize=8, width=14, height=13)

for(i in 1:54){
	
	#[1] - Add data for the gene 
	gene <- colnames(pan)[i]
	#rows = #patients in pan (2546) + #patients in pcawg (497) = 3043
	all_gene_data <- as.data.frame(matrix(ncol=5, nrow=3043))
	colnames(all_gene_data) <- c("Gene", "Dataset", "Patient", "Cancer", "Expression")
	all_gene_data[,1] <- gene
	#add pancan data first
	z <- which(is.na(all_gene_data[,2]))[1]
	all_gene_data[z:(z+2545),2] <- "PanCan"
	dat <- pan[,c(i,55)]
	all_gene_data[z:(z+2545),3] <- rownames(dat)
	all_gene_data[z:(z+2545),4] <- dat$canc
	all_gene_data[z:(z+2545),5] <- dat[,1]
	#add pcawg data second
	z <- which(is.na(all_gene_data[,2]))[1]
	all_gene_data[z:(z+496),2] <- "PCAWG"
	dat <- pcawg[,c(i,55)]
	all_gene_data[z:(z+496),3] <- rownames(dat)
	all_gene_data[z:(z+496),4] <- dat$canc
	all_gene_data[z:(z+496),5] <- dat[,1]
	#log gene expression
	all_gene_data$Expression <- log1p(all_gene_data$Expression)

	#[2] - Plot boxplot of gene's Expression for each dataset - side by side
	g <- ggviolin(all_gene_data, x="Cancer", y="Expression", fill="Cancer", palette=mypal, facet.by="Dataset", draw_quantiles = 0.5,
		order = unique(all_gene_data$Cancer)[c(1,3,7,2,4,5,6)])
	g <- ggpar(g, font.legend = c(6, "plain", "black")) 
	g <- g + labs(title = paste(gene, "Expression Distribution"), y="log1p(Expression)") + 
     	 theme(plot.title = element_text(hjust = 0.5))
	g <- g + rremove("x.text")
	g <- g + rremove("x.ticks")  
	g <- g  + stat_summary(fun.y=median, geom="smooth", aes(group=1))

	print(g)
	
}#end loop

dev.off()

