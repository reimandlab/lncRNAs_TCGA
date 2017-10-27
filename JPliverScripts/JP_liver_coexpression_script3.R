#---------------------------------------------------------
#JP_liver_coexpression_script3.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: July 4th 2017
#dir: Thesis/pcawg_liver_JP

#------------
#Description:
#------------

#Have two RNA-Seq files, one for lncRNAs (all from UCSC including low confidence) 
#and one for PCGs. Conduct here LM and NGB regression to identify 
#signficantly co-expressed lncRNA-PCGs in liver cancer patients 
#using an array of confounders 
#------------
#confounders:
#------------
#[1]. PCG CNA
#[2]. lncRNA CNA
#[3]. Clinical features 
#[4]. is methylation available?

#script3 - identify mutation status of lncRNAs/pcgs 


#---------------------------------------------------------
#Preamble
#---------------------------------------------------------

options(stringsAsFactors=F)

#---------------------------------------------------------
#Libraries
#---------------------------------------------------------
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

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#lncRNAs from Script 2
lnc <- readRDS("lnc_liver_expression_file.rds")

#PCGs from Script 2 
all <- readRDS("pcg_liver_expression_file.rds")

#clinical data
clin <- fread("300wgs_patient_data.txt", data.table=F)
#subset to include 68 patients that have gene expression data 
pats <- colnames(lnc[,1:68])
z <- which(clin$ID %in% pats)
clin <- clin[z,]
clin[1,23] <- 1 #fix reading in error 

#---------------------------------------------------------
#Co-Expression Analysis - Spearman
#---------------------------------------------------------

#set up matrix with columns (possible predictors) as lncRNAs with another column for a single PCG 
#also add clinical features histology and virus present 

#tier1 lincRNAs first 
lncs <- subset(lnc, lnc$tier == "TIER1" & lnc$lnc_type=="lincRNA")

#set up data frame that will be used for analysis 
genes <- t(lncs[,1:68])
genes <- as.data.frame(genes)

#start with tier1 PCGs
pcgs_tier1 <- all[all$tier=="TIER1",]

##PCA plot visualize how lncRNA expression variables 
##predict seperation of various clinical features 
library("FactoMineR")
library("factoextra")

#for pca, columns 1:723
pca <- PCA(log1p(genes[,1:293]), graph = FALSE)
fviz_pca_ind(pca, label="none")
fviz_pca_ind(pca,  label="none", habillage=as.factor(genes$histo))
fviz_pca_ind(pca, label="none", habillage=as.factor(genes$histo),
             addEllipses=TRUE, ellipse.level=0.55)

pcgs_tier1$pcg <- rownames(pcgs_tier1)

spear_pcg <- function(row){

	##Data-SetUp
	genes$PCG <- as.numeric(row[1:68])
	pcg_name <- row[[70]]
	genes$virus <- ""
	genes$histo <- ""

	for(i in 1:nrow(genes)){
		pat <- rownames(genes)[i]
		z <- which(clin$ID %in% pat)
		genes$virus[i] <- clin[z,6]
		genes$histo[i] <- clin[z,12]
		}

	genes$tag <- ""
	med <- median(genes$PCG)
	for(i in 1:nrow(genes)){
		h <- genes$PCG[i]
		z <- h >= med
		if(z){
			genes$tag[i] <- "HIGH"
		}
		if(!(z)){
			genes$tag[i] <- "LOW"
		}
		}

	##Spearman-Correlations
	##One PCG * ALL lncRNAs
	m <- rcorr(as.matrix(genes[,1:294]), type="spearman")
	pcg <- m$r[,294]
	pcg_p <- m$P[,294]
	pcg_p <- p.adjust(pcg_p, method="fdr")
	pcg <- cbind(pcg, pcg_p)
	pcg <- as.data.frame(pcg)
	pcg$sig <- ""
	z <- which(pcg$pcg_p <= 0.05)
	pcg$sig[z] <- "Significant"
	pcg$sig[-z] <- "Not"
	pcg <- pcg[-294,] 
	pcg$id <- 1:293
	pcg$sig <- as.factor(pcg$sig)

	##Plot correlations 
	mypal = pal_npg("nrc", alpha = 0.7)(10)
	name <- paste("sig_lncs_associated_wPCGs_plots/", pcg_name, "lncs.pdf", sep="_")
	pdf(name, pointsize=6)
	print(ggdotchart(pcg, x="id", y="pcg", ggtheme=theme_bw(), group="sig", color="sig", rotate=TRUE, ylab=FALSE,
	palette = mypal))
	dev.off()

	pos <- pcg[pcg$pcg>0,]
	neg <- pcg[pcg$pcg<0,]
	
		print(ggdotchart(pos, x="id", y="pcg", ggtheme=theme_bw(), group="sig", color="sig", rotate=TRUE, ylab=FALSE,
	palette = mypal))
			print(ggdotchart(neg, x="id", y="pcg", ggtheme=theme_bw(), group="sig", color="sig", rotate=TRUE, ylab=FALSE,
	palette = mypal))


	#save list of lncRNAs and results for the PCG
	sig_genes <- rownames(pcg)[pcg$sig=="Significant"]
	data_add <- c(pcg_name, length(sig_genes), paste(pcg_name, "lncs", sep="_"))
	data_add <- as.data.frame(matrix(data=data_add, ncol=3))

	name <- paste(pcg_name, "spear_results.bed", sep="_")
	write.table(data_add, name, quote=F)
	#save list of sig_genes
	name <- paste("sig_lncs_associated_wPCGs_lists/", pcg_name, "lncs.txt", sep="_")
	write.table(sig_genes, name, quote=F, row.names=F)
	return(pcg_name)

}#function end 

apply(pcgs_tier1, 1, function(x) spear_pcg(x))

#coerce all bed files together into one data.frame
#result indicating how many sig spearman correlations 
#lncRNA to PCG
p <- fread("all_pcg_lncs_results.txt"  , data.table=F, header=F)
z <- which(p[,1] == "V1")
p <- p[-z,]
p <- p[,-1]
colnames(p) <- c("PCG", "num_lncs_assoc", "id")
#remove duplicates
z <- which(duplicated(p$PCG))
p <- p[-z,]

##List of Transcription factors from consensus paper nature 2009
tfs <- fread("TF_censushumantranscription.txt", data.table=F)
#keep only a and b
tfs <- subset(tfs, tfs[,1] %in% c("a", "b"))

z <- which(p$PCG %in% tfs[,2])
p$type <- ""
p$type[z] <- "TF"
p$num_lncs_assoc <- as.numeric(p$num_lncs_assoc)	

mypal = pal_npg("nrc", alpha = 0.7)(10)

liver_tf <- fread("TFs_liver.txt", data.table=F)
z <- which(p$PCG %in% liver_tf[,2])
p$type[z] <- "Liver_TF"

pdf("spearman_correlations_result_distribut.pdf", pointsize=6, width=8, height=9)
 gghistogram(p, x="num_lncs_assoc", add="median", rug=FALSE, color="type", palette=mypal, bins=50, xlab= "Number of Signficiant Co-Expression lncRNAs", 
 	ylab="Count")
dev.off()


