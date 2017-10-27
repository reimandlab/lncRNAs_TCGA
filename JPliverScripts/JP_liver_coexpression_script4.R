#---------------------------------------------------------
#JP_liver_coexpression_script2.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: June 28th 2017
#dir: /pcawg_liver_JP/sig_lncs_associated_wPCGs_lists

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

#script4 - make node and edge table using lists of lncrnas 
#obtained from spearman correlation analysis of each PCG 

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
#Co-Expression Analysis - lncRNAs results 
#---------------------------------------------------------

#1. Make node file - aka - obtain list of unique lncRNAs (in lists)
#and PCGS (names of files)

all <- list.files()

edges <- as.data.frame(matrix(ncol=3))
colnames(edges) <- c("source", "target", "interaction")

for(i in 1:length(all)){
	file <- all[i]
	pcg <- substr(file, 2, 16)
	f <- fread(file, data.table=F)
	if(!(dim(f)[1]==0)){
	f$target <- pcg
	colnames(f)[1] <- "source" 
	f$interaction <- "co-expression"
	edges <- rbind(f, edges)
}}

edges <- edges[-882611,] #remove NA line 

lncs <- unique(edges$source)
lncs <- as.data.frame(lncs)
lncs$type <- "lncRNA"
pcgs <- unique(edges$target)
pcgs <- as.data.frame(pcgs)
pcgs$type <- "PCG"
colnames(pcgs)[1] <- "gene"
colnames(lncs)[1] <- "gene"
nodes <- rbind(pcgs, lncs)
colnames(nodes) <- c("id", "group")

##smaller subset 
edges <- edges[1:5000,]

#---------------------------------------------------------
#Get number of PCGs per lncRNA then check how many lncRNAs
#have PCGs in common
#---------------------------------------------------------

lnc_counts <- as.data.frame(table(edges$source))
colnames(lnc_counts) <- c("lnc", "Num_PCGs")

mypal = pal_npg("nrc", alpha = 0.7)(10)

pdf("spearman_correlations_result_lncs_distribut.pdf", pointsize=6, width=8, height=9)
gghistogram(lnc_counts, x="Num_PCGs", add="median", rug=FALSE, bins=50, xlab= "Number of Signficiant Co-Expressed PCGs", ylab="Number of lncRNAs", palette=mypal, color=mypal[3], fill=mypal[2])
dev.off()





































#+++++++++
#TEST - LM------------------------------------------------
#+++++++++

#++++
#STEPS----------------------------------------------------
#++++

#1.
#Data standardzation: center response y and center ands scale each explanatory variable to have meu=0 and sd=1
#response = PCG E, explanatory = lncRNA E (one at a time for now)

#For every lncRNA, all PCGs 
lncE <- lnc[1,1:68]
lncEE <- scale(as.numeric(lncE))
pcg <- all[1,1:68]
pcg <- scale(as.numeric(pcg), center=T, scale=F)

#2. 
#Linear Regression 
test <- cbind(lncEE, pcg)
test <- as.data.frame(test)
rownames(test) <- names(lncE)
colnames(test) <- c("lncE", "pcgE")

#what other predictors do we want to include? 
#virus present and histological subtype 
test$virus <- ""
test$histo <- ""

for(i in 1:nrow(test)){
	pat <- rownames(test)[i]
	z <- which(clin$ID %in% pat)
	test$virus[i] <- clin[z,6]
	test$histo[i] <- clin[z,12]
}

t <- lm(test[,2] ~ test[,1] , data= test)

#Plotting assumptions 
par(mar = c(4, 4, 2, 2), mfrow = c(1, 2)) 
plot(t, which = c(1, 2)) 

#+++++++++
#TEST - NB------------------------------------------------
#+++++++++

#++++
#STEPS----------------------------------------------------
#++++

#For every lncRNA, all PCGs 
lncE <- lnc[1,1:68]
lncEE <- as.numeric(lncE)
pcg <- as.numeric(all[1,1:68])

#2. 
#glm.nb Regression 
test <- cbind(lncEE, pcg)
test <- as.data.frame(test)
rownames(test) <- names(lncE)
colnames(test) <- c("lncE", "pcgE")

#what other predictors do we want to include? 
#virus present and histological subtype 

glm.nb(floor(test$pcgE) ~ floor(test$lncE), data=test)