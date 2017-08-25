#Karina Isaev
#August 25th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#wd - Thesis/top5_cancers_pcawg_thesis_jult2017/raw/binned_lncRNAs_wGTEX/lncs_binned_ranked

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library("colorout")
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
library(ggthemes)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------
pcawg <- list.files(pattern="lncRNAsPCAWG.txt")
gtex <- list.files(pattern="lncRNAs_GTEX.txt")

pdf("binned_lncRNAs_compared_GTEX_PCAWG_1percentile_at_aTime.pdf", pointsize=8, width=16, height=13)

for(i in 1:length(pcawg)){
	dat <- fread(pcawg[i], sep=";")
	dat <- as.data.frame(dat)
	dat$type <- "PCAWG"
	tis <- unlist(strsplit(dat$Tissue, " "))[1]
	dat_g <- fread(gtex[grep(tis, gtex)], sep=";")
	dat_g <- as.data.frame(dat_g)
	dat_g$type <- "GTEX"
	
	dat_g <- subset(dat_g, dat_g$Gene %in% dat$Gene)
	dat <- subset(dat, dat$Gene %in% dat_g$Gene)
	
	#divide by groups of percentiles
	for(y in 1:100){
		#percs <- (y*10-9):(y*10)
		percs <- y
		dat_pcawg_plot <- dat[dat$Percent %in% percs,]
		dat_g_plot <- subset(dat_g, dat_g$Gene %in% dat_pcawg_plot$Gene)
		to_plot <- rbind(dat_pcawg_plot, dat_g_plot)
		order <- dat_pcawg_plot$Gene
		g <- ggline(to_plot, x="Gene", y="Percent", linetype="type", palette=mypal, color="type", order=order, title=paste("Percentile Group", y, dat_pcawg_plot$Tissue[1]), ylab="Percentile")
		g <- ggpar(g, font.tickslab = c(8,"plain", "black"), xtickslab.rt = 75)
		print(g)
	}
}

dev.off()