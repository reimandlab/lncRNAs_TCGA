###coexpression_script0_extraction.R

#Author: Karina_Isaev
#Date_started: July 11th 2017
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

#------------------------------------------------------------------
#Within each tissue type, score genes 
#------------------------------------------------------------------

#1. log1p
rna[,1:(ncol(rna)-2)] <- log1p(rna[,1:(ncol(rna)-2)])

#2. Get lncRNA - median within each tissue type
tissues <- unique(rna$tissue)

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- rna[rna$tissue==tissue,]
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 
getScores <- function(row){
	score=""
	expression <- data.frame(exp=as.numeric(row[1:(length(row)-2)]), gene=names(row)[1:(length(row)-2)])
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lnc candidates - we just want their score 
	z <- which(expression$gene %in% c(allCands$gene, "TP53"))
	expression <- expression[z, ]
	return(expression)
}

addScores <- function(dataframe){
	patients <- apply(dataframe, 1, getScores) #list of dataframes, need to coerce together
	names <- rownames(dataframe)
	patients <- rbindlist(patients)
	patients$patient <- rep(names, each=39) #38 lncRNA candidates + TP53
	patients <- as.data.frame(patients)
	patients$canc <- dataframe$tissue[1]
	patients$data <- "GTEX"
	return(patients)
}	

scored <- lapply(tissues_data, addScores) #list of dataframes
all_tissues_scored <-  rbindlist(scored)

#bind all tissues scored from GTEX with all cancers scored from pcawg
all_cancers_scored <- readRDS("all_cancers_scoredPLUS_TP53.rds")

allScoredLncsgtexANDpcawg <- rbind(all_cancers_scored, all_tissues_scored)
allScoredLncsgtexANDpcawg <- as.data.frame(allScoredLncsgtexANDpcawg)
allScoredLncsgtexANDpcawg$canc <- as.character(allScoredLncsgtexANDpcawg$canc)

#Make subsetted version that only has the lnc candidate in respective cancer type 
#plot the gtex equivalent, gene and Canc need to match 
#each gene has just one HR and pvalue 

#Significant 
row <- as.data.frame(matrix(nrow=1, ncol=6))
row[1,] <- c("TP53", 0.4367, 0.6462,0.0983, "Ovary", 0.15)
colnames(row) <- colnames(sig)
sig <- rbind(row, sig)

cands_for_plottingV3 <- merge(allScoredLncsgtexANDpcawg, sig, by=c("gene", "canc"))

#------------------------------------------------------------------
#PLOTS 
#------------------------------------------------------------------

#PLOT 1 - x-axis lncRNA candidate genes , y-axis score, stratify by cancer type 
#all genes/cancer type 
#pdf("plot1_pcawgVSgtex_38candidates.pdf", pointsize=8, width=25, height=25)
#f <- ggboxplot(allScoredLncsgtexANDpcawg, x="gene", y="score", color="data", palette=mypal)
#f <- facet(f, facet.by="canc", nrow=7, ncol=1)
#f <- ggpar(f, xlab="Candidate lncRNAs, n=38", ylab="Score", x.text.angle=65, font.tickslab=c(10, "plain", "black"), legend="right", ylim=c(0,1))
#f
#dev.off()

#plot 1 version 2
#plot only candidate lncRNAs within that cancer type
allCands$canc <- lapply(allCands$canc, function(x) unlist(strsplit(x, " "))[1])
row <- as.data.frame(matrix(nrow=1, ncol=6))
row[1,] <- c("TP53", 0.4367, 0.6462,0.0983, "Ovary", 0.15)
colnames(row) <- colnames(allCands)
allCands <- rbind(row, allCands)

tissues <- unique(allScoredLncsgtexANDpcawg$canc) #keep only candidates pancreas, kidney, liver and ovary 
tissues <- tissues[c(3,1,4,5)]
plots <- list()
wilcoxon_results <- list()

for(i in 1:length(tissues)){
	allScored <- allScoredLncsgtexANDpcawg[allScoredLncsgtexANDpcawg$canc == tissues[i],]
	cands <- filter(allCands, canc==tissues[i])
	allScored <- allScored[which(allScored$gene %in% cands$gene),]
	#get median score for each gene 
	pscore <- filter(allScored, data=="PCAWG")
	meds <- as.data.table(pscore %>%
  		group_by(gene) %>%
  		summarise_each(funs(median), score))

  	meds <- meds[order(score)]

  	#wilcoxon test, gene score between gtex and pcawg 
  	wil <- as.data.frame(allScored %>% group_by(gene) %>% do(tidy(wilcox.test(score ~ data, data=.))))
  	wil$fdr <- p.adjust(wil$p.value, method="fdr")
  	wil$sig <- "" 
  	wil$sig[wil$fdr <=0.05] <- "Yes"
  	wil$sig[wil$fdr > 0.05] <- "No"
  	wil$canc <- allScored$canc[1]

	f <- ggboxplot(allScored, x="gene", y="score", color="data", palette=mypal, order=as.character(meds$gene))
	f <- ggpar(f, xlab="Candidate lncRNAs", main= paste(length(unique(allScored$gene)), "Candidate Genes in", allScored$canc[1], "Cancer") ,ylab="Score", x.text.angle=65, font.tickslab=c(10, "plain", "black"), legend="right", ylim=c(0,1))
	plots[[i]] <- f 
	wilcoxon_results[[i]] <- wil	
}

#pdf("plot1_Version3_pcawgVSgtex_38candidatesWithinEachCancer.pdf", pointsize=8, width=30, height=20)
#multiplot(plotlist = plots, cols = 2)
#dev.off()

wilcoxon_results <- rbindlist(wilcoxon_results)
wilcoxon_results <- as.data.table(wilcoxon_results)
wilcoxon_results$fdr <- as.numeric(wilcoxon_results$fdr)
wilcoxon_results <- wilcoxon_results[order(fdr)]

pdf("WithP53_42_sig_lncRNA_predictors_usingMedian5WILCOXON_gtexVSpcawg.pdf", pointsize=8, width=12, height=14)
p<-tableGrob(wilcoxon_results)
grid.arrange(p)
dev.off()

#plot 1 version 3
#only plot candidate lncs within that cancer type and facet with xspace = free 
#the next plot underneath would show HRs 

#Order cands_for_plottingV3 by median score of each lncRNA with PCAWG cancer 
tissues <- unique(cands_for_plottingV3$canc) #keep only candidates pancreas, kidney, liver and ovary 
new_matrix <- as.data.frame(matrix(ncol=15))
colnames(new_matrix) <- c(colnames(cands_for_plottingV3), "Hazard", "medianScore", "high", "Prediction", "Tier")

#change names of duplicated candidate gens
dup_data <- sig[which(duplicated(sig$gene)),]
nonSigWilcox <- filter(wilcoxon_results, fdr > 0.05)

for(i in 1:length(tissues)){
	allScored <- cands_for_plottingV3[cands_for_plottingV3$canc == tissues[i],]
	#get median score for each gene 
	pscore <- filter(allScored, data=="PCAWG")
	pscore$score <- as.numeric(pscore$score)
	meds <- as.data.table(pscore %>%
         group_by(gene) %>%
         summarise_each(funs(median), score))
  	meds <- meds[order(score)]
  	#select all rows from original dataset that correpond to these genes and cancer type
  	z <- which(cands_for_plottingV3$canc == tissues[i])
	df <- cands_for_plottingV3[z,]
	df$Hazard <- ""
	df$Hazard[df$HR <1] <- "TS"
	df$Hazard[df$HR >1] <- "OG"

	#compare median of GTEX to median of PCAWG
	#if median in GTEX > median PCAWG and HR <1, prediction = TS
	predictions <- as.data.frame(df %>% group_by(gene, data) %>%  summarise_each(funs(median), score))
	colnames(predictions)[3] <- "medianScore"
	predictions$high <- ""
	predictions$HR <- ""
	for(j in 1:length(unique(predictions$gene))){
		check <- predictions[predictions$gene == unique(predictions$gene)[j],]
		z <- check$data[which(check$medianScore == max(check$medianScore))]
		predictions$high[predictions$gene == unique(predictions$gene)[j]] <- z
		z2 <- df$HR[which(df$gene == unique(predictions$gene)[j])][1]
		predictions$HR[predictions$gene == unique(predictions$gene)[j]] <- z2
	}
	predictions$Prediction <- ""
	predictions$Prediction[(predictions$high == "GTEX") & (predictions$HR < 1)] <- "Predicted TS"
	predictions$Prediction[(predictions$high == "PCAWG") & (predictions$HR > 1)] <- "Predicted OG"
	predictions$Prediction[(predictions$high == "PCAWG") & (predictions$HR < 1)] <- "No Prediction"
	predictions$Prediction[(predictions$high == "GTEX") & (predictions$HR > 1)] <- "No Prediction"
	predictions <- predictions[,-2]

	#add to df
	df <- merge(df, predictions, by=c("gene", "HR"))

	#if fdr < 0.1, tier == 1
	df$Tier[df$fdr <=0.1] <- 1
	df$Tier[df$fdr >0.1] <-  2

	df <- df[ order(match(df$gene, meds$gene)), ]

	#if gene name is duplicated change it to something else 
	z <- which((df$gene %in% dup_data$gene )& (df$canc %in% dup_data$canc))
	replace <- paste(df$gene[z], "*")
	df$gene[z] <- replace

	new_matrix <- rbind(new_matrix, df)

  }

new_matrix <- new_matrix[-1,]
new_matrix$pval <- as.numeric(new_matrix$pval)
new_matrix$pval <- -log10(new_matrix$pval)

#BOXPLOTS

#f <- ggboxplot(new_matrix, x="gene", y="score", color="data", fill="data", palette=mypal[c(3,4)], ggtheme=theme_bw())
#f <- f + facet_grid (.~ canc, scales = "free_x", space = "free_x") + 
 	#theme(strip.background =element_rect(fill=mypal[9]))+
  	#theme(strip.text = element_text(colour = 'white'))

#f <- f + stat_compare_means(aes(group=data, label = ..p.signif..))
#f <- facet(f, facet.by = "canc", scales="free_x")
#f <- ggpar(f, xlab="Candidate lncRNAs", ylab="Rank", x.text.angle=35, font.tickslab=c(10, "plain", "black"), legend="right", ylim=c(0,1))
#f <- f + rremove("y.grid") 
#+ rremove("xlab") + rremove("x.text")

#old boxplots
f <- ggboxplot(new_matrix, x="gene", y="score", color="data", fill="data", palette=mypal[c(3,4)], ggtheme=theme_bw())
f <- f + facet_grid (.~ canc, scales = "free_x", space = "free_x") + 
 	theme(strip.background =element_rect(fill=mypal[9]))+
  	theme(strip.text = element_text(colour = 'white'))
f <- ggpar(f, xlab="Candidate lncRNAs", ylab="Rank of Expression Relative to All Genes", x.text.angle=65, font.tickslab=c(10, "plain", "black"), legend="top", ylim=c(0,1))
f <- f + rremove("y.grid")
pdf("justBoxplots_comparingGTEX_PCAWG_inclTP53.pdf", pointsize=8, width=18, height=10)
f
dev.off()

#Instead of boxplots, label each point as gtex or pcawg, showing in which dataset it's more highlyE

#HAZARD RATIOS
new_matrix$HR <- as.numeric(new_matrix$HR)
colnames(new_matrix)[13] <- "HigherExpression"
p <- ggscatter(new_matrix, x="gene", y="HR", palette= mypal[c(2,1)], ggtheme=theme_bw(), size="pval", color="HigherExpression")
p <- p + facet_grid (.~ canc, scales = "free_x", space = "free_x") + 
 	theme(strip.background =element_rect(fill=mypal[9]))+
  	theme(strip.text = element_text(colour = 'white'))
p <- ggpar(p, xlab="Candidate lncRNAs", ylab="Hazard Ratio", x.text.angle=65, font.tickslab=c(10, "plain", "black"), legend="right", ylim=c(0,4))
p <- p + rremove("y.grid") + rremove("xlab") + rremove("x.text")
p <- p + geom_hline(yintercept = 1)


#TIER AND PREDICTION STATUS 
new_matrix$Tier <- as.factor(new_matrix$Tier)

g <- ggscatter(new_matrix, x="gene", y="Tier", palette= c("grey", mypal[c(1,2)]), ggtheme=theme_bw(), color="Prediction")
g <- g + facet_grid (.~ canc, scales = "free_x", space = "free_x") + 
 	theme(strip.background =element_rect(fill=mypal[9]))+
  	theme(strip.text = element_text(colour = 'white'))
g <- ggpar(g, xlab="Candidate lncRNAs", ylab="Tier", x.text.angle=55, font.tickslab=c(11, "plain", "black"), legend="right")
g <- g + rremove("y.grid")


pdf("WITHp53_newplots_plot1_pcawgVSgtex_38candidatesV4ordered_justcandsPLUShazardratios.pdf", pointsize=8, width=18, height=8)
plot_grid(p, g,  labels = c("A", "B"), align = "v", nrow = 2, rel_heights = c(1.6,1))
dev.off()























