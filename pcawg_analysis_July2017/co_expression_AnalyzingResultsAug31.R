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

mypal = pal_npg("nrc", alpha = 0.7)(10)


#Data#----------------------------------------------------

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

all <- readRDS(file="all_coexpression_for42candidateLNCS.rds")
high <- readRDS(file="HIGH_coexpression_for42candidateLNCS.rds")
low <- readRDS(file="LOW_coexpression_for42candidateLNCS.rds")

#list of transcription factors 
tfs <- fread("TF_censushumantranscription.txt", fill=TRUE, header=F)
tfs <- tfs$V6


#What do I want to know: Differences between high and low lncRNA expressing patients 
#So within each cancer type I want to look at list of genes
#Obtained to be co-expressed with lncRNA when all patients grouped together, high lncRNA expressing group only or low lncRNA expressing group only 
#Can also get list of pathways 

#[1] COMBINE ALL RESULTS
results <- rbind(all, high, low)

#add coef status, positive or negative
results$lm_coefficient <- as.numeric(results$lm_coefficient)
results$cor[results$lm_coefficient > 0] <- "Positive"
results$cor[results$lm_coefficient < 0] <- "Negative"

#Get some more info on lncRNA in this study
cands <- fantom[which(fantom$CAT_geneName %in% results$lnc), ]

#[2] QUESTIONS TO ANSWER
# - which cancer had the most significantly co-expressed genes 
# - how many was each candidate lncRNA significantly associated with?
	#-when all patients combined
	#-when only high
	#-when only low 
# - within a cancer type
	#- was there overlap between candidate lncRNAs co-expressed genes?


#Function 1 
#Divide resutls by cancer type

tissues <- unique(results$canc)

cancerResult <- function(tissue){
	df <- results[results$canc == tissue, ]
	return(df)
}

tissueData <- llply(tissues, cancerResult, .progress = "text")

#Function 2
#Make three rows of plots for each Patients group type, showing coefficient versus -log10pvalue 

plotDist <- function(df){
	df$lm_coefficient <- as.numeric(df$lm_coefficient)
	df$logged <- -log10(as.numeric(df$lm_anov_pval))
	plotsR <- df %>% group_by(lnc) %>% do(plots=ggplot(data=.) + 
         aes(x=lm_coefficient, y=logged) + geom_point(aes(colour = logged)) + scale_colour_gradient(low = "blue", high="red") +
         labs(colour = "-log10(p-value)", x = "lncRNA Linear Model Coefficient", y= "-log10(p-value)") + ggtitle(paste(unique(.$lnc), unique(.$canc))) + facet_grid(. ~ Patients))

	invisible(lapply(plotsR$plots, print))

}

#pdf("all_cancers_lncRNAcoExpressionResultsAugust31.pdf", pointsize=8, width=15, height=8)
l <- llply(tissueData, plotDist, .progress = "text")
dev.off()

#Function 3
#Barplot for each lncRNA-cancer facetted by Patients 
#y = Number of PCGs per lncRNA 

plotBar <- function(df){

	plot <- as.data.table(table(df$lnc, df$Patients))  
	plot <- plot[order(N)]
	p <- ggbarplot(plot, x= "V1", y ="N", color = "V2", fill = "V2", palette=mypal, position = position_dodge(0.9))	
	p <- ggpar(p, xlab="Candidate lncRNAs", ylab="Number of Co-Expressed PCGs", x.text.angle=65, font.tickslab=c(10, "plain", "black"), 
		legend="right", legend.title = "Cohort", title=paste(df$canc[1], "Co-Expression Results"),  yticks.by=200)
	print(p)
}

#pdf("all_cancers_Barplots_lncRNAcoExpressionResultsAugust31.pdf", pointsize=8, width=14, height=12)
l2 <- llply(tissueData, plotBar, .progress = "text")
dev.off()

#Function 4
#Barplot for each lncRNA-cancer facetted by Patients 
#y = Number of PCGs per lncRNA , group by positive vs negative relationships

plotBarwCors <- function(df){
	plot <- as.data.table(table(df$lnc, df$cor, df$Patients, df$cor))  
	#plot <- plot[order(N)]
	plot <- plot[order(N)]
	p <- ggbarplot(plot, x= "V3", y ="N", color = "V2", fill = "V2", palette=mypal[c(2,1)])
	p <- facet(p, facet.by = "V1", ncol=4)	
	p <- ggpar(p, xlab="Candidate lncRNAs", ylab="Number of Co-Expressed PCGs", x.text.angle=65, font.tickslab=c(10, "plain", "black"), 
		legend="right", legend.title = "Correlations", title=paste(df$canc[1], "Co-Expression Results"),  yticks.by=500)
	print(p)

}

pdf("SplitbyCorrelationall_cancers_Barplots_lncRNAcoExpressionResultsAugust31.pdf", pointsize=8, width=23, height=18)
l3 <- llply(tissueData, plotBarwCors, .progress = "text")
dev.off()


#Function 5
#Barplot for each lncRNA-cancer facetted by Patients 
#Find how many PCGs overlap between groups 

OverlapplotBar <- function(df){

	plot <- as.data.table(table(df$lnc, df$PCG))  
	plot <- plot[order(N)]
	plot <- filter(plot, N >=1)
	p <- ggbarplot(plot, x= "V1", y ="N", color = "V2", fill = "V2", palette=mypal, position = position_dodge(0.9))	
	p <- ggpar(p, xlab="Candidate lncRNAs", ylab="Number of Co-Expressed PCGs", x.text.angle=65, font.tickslab=c(10, "plain", "black"), 
		legend="right", legend.title = "Cohort", title=paste(df$canc[1], "Co-Expression Results"),  yticks.by=200)
	print(p)
}

#pdf("all_cancers_Barplots_lncRNAcoExpressionResultsAugust31.pdf", pointsize=8, width=14, height=12)
l4 <- llply(tissueData, OverlapplotBar, .progress = "text")
dev.off()

















