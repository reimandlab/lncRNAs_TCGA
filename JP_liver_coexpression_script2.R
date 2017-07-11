#---------------------------------------------------------
#JP_liver_coexpression_script2.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: June 28th 2017
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

#script2 - more processing and visualization via HEATMAPS

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

#These are lists of genes 
tier1_lncs <- fread("tier1_lncs.txt", data.table=F)
tier2_lncs <- fread("tier2_lncs.txt", data.table=F)
tier1_pcgs <- fread("tier1_pcgs.txt", data.table=F)
tier2_pcgs <- fread("tier2_pcgs.txt", data.table=F)
lncs_medians4 <- fread("meds_greaterthan4.txt", data.table=F)

#These are gene expression files  
lnc <- readRDS("liver_jp_lncRNA_expression_6028.rds")
all <- readRDS("liver_jp_pcg_expression.rds")
#change feature column with gene names so that they are the rownames
rownames(lnc) <- lnc[,1] ; lnc <- lnc[,-1] #13479 lncRNAs as defined by ensembl "lincRNA", "antisense", "sense_intronic", "sense_overlapping"
rownames(all) <- all[,1] ; all <- all[,-1] #18039 PCGs as defined by ensembl "protein_coding"

#clinical data
clin <- fread("300wgs_patient_data.txt", data.table=F)
#subset to include 68 patients that have gene expression data 
pats <- colnames(lnc)
z <- which(clin$ID %in% pats)
clin <- clin[z,]
clin[1,23] <- 1

#ucsc gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#[1]------------------------------------------------------
#SUBSET LNC, ALL, TO INCLUDE ONLY GENES IN EITHER
#TIER1 OR TIER2 RESPECTIVELY 
#---LNC----------------------------------
lncs <- c(tier1_lncs[,1], tier2_lncs[,1])
z <- which(rownames(lnc) %in% lncs)
lnc <- lnc[z,]
#---PCG----------------------------------
pcgs <- c(tier1_pcgs[,1], tier2_pcgs[,1])
z <- which(rownames(all) %in% pcgs)
all <- all[z,]

#[2]------------------------------------------------------
#ADD NEW COLUMN INDICATING WHICH TIER 
#GENE BELONGS TO
#---LNC----------------------------------
lnc$tier <- ""
for(i in 1:nrow(lnc)){
	gene <- rownames(lnc)[i]
	z <- which(tier1_lncs[,1] %in% gene)
	if(length(z)==0){
		lnc$tier[i] <- "TIER2"
	}
	if(!(length(z)==0)){
		lnc$tier[i] <- "TIER1"
	}
}


#---PCG----------------------------------
all$tier <- ""
for(i in 1:nrow(all)){
	gene <- rownames(all)[i]
	z <- which(tier1_pcgs[,1] %in% gene)
	if(length(z)==0){
		all$tier[i] <- "TIER2"
	}
	if(!(length(z)==0)){
		all$tier[i] <- "TIER1"
	}
}

#[3]------------------------------------------------------
#FOR LNCRNAS, ADD WHAT KIND OF LNCRNA THEY ARE 
#USING UCSC FILE
lnc$lnc_type <- ""
for(i in 1:nrow(lnc)){
	id <- rownames(lnc)[i]
	z <- which(ucsc[,6] %in% id)
	lnc$lnc_type[i] <- ucsc[z,7]
}


#---------------------------------------------------------
#Save Data 
#---------------------------------------------------------

saveRDS(lnc, file= "lnc_liver_expression_file.rds")
saveRDS(all, file= "pcg_liver_expression_file.rds")

#---------------------------------------------------------
#Visualize 
#---------------------------------------------------------

mypal = pal_npg("nrc", alpha = 0.7)(10)

#[1]------------------------------------------------------
#BARPLOT OF HOW MANY OF EACH LNCRNAS ARE PRESENT 
types <- as.data.frame(table(lnc$lnc_type))
colnames(types) <- c("Type", "Number")
#Script 2 Figure 1 
#+++++++++++++++++
p1 <- ggbarplot(types, "Type", "Number", fill="Type", color = "Type", palette = mypal, ylab="Number of Genes", label = TRUE, label.pos = "out", legend="none") + scale_x_discrete(labels=c("Antisense", "LincRNA", "Sense Intronic", "Sense Overlapping")) + 
theme(axis.ticks.x=element_blank()) 
#BARPLOT OF HOW MANY TIER/TYPE COMBINATIONS THERE ARE 
types <- as.data.frame(table(lnc$tier, lnc$lnc_type))
colnames(types) <- c("Tier", "Type", "Number")
p2 <- ggbarplot(types, "Type", "Number", fill="Tier", color = "Tier", palette = mypal, legend="right", ylab="Number of Genes", label = TRUE, position = position_dodge(0.65)) + scale_x_discrete(labels=c("Antisense", "LincRNA", "Sense Intronic", "Sense Overlapping")) + 
theme(axis.ticks.x=element_blank())
#COMBINE
#+++++++++++++++++
pdf("type_of_lncs_studied_inctier_new_tier1.pdf", pointsize=5, width=9, height=9)
grid.arrange(p1, p2, nrow = 2)
dev.off()

#[2]------------------------------------------------------
#HIEARCHIAL CLUSTERING OF PATIENTS BY LNCRNA EXPRESSION 
#TIER1 FIRST 

#first log expression values
lnc[,1:68] <- log1p(lnc[,1:68])
all[,1:68] <- log1p(all[,1:68])

#All lncRNAs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hmcol <- colorRampPalette(brewer.pal(9, "BrBG"))(100)
rv <- rowVars(lnc[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 
cols <- palette(mypal)[as.fumeric(clin[,12])]
head(cbind(colnames(lnc[,1:68]),cols))

pdf("lncs_liver_clustered_heatmap_all_lncs.pdf", pointsize=6)
heatmap.2(as.matrix(lnc[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          margins =c(12,9),
          col=hmcol,
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()

#Tier 1 lncRNAs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tier1_lncs <- lnc[lnc$tier=="TIER1",]
rv <- rowVars(tier1_lncs[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 

pdf("lncs_liver_clustered_heatmap_TIER1_lncs.pdf", pointsize=6)
heatmap.2(as.matrix(tier1_lncs[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          col=hmcol,
          margins =c(12,9),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()

#Tier 2 lncRNAs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tier2_lncs <- lnc[lnc$tier=="TIER2",]
rv <- rowVars(tier2_lncs[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 
hmcol <- colorRampPalette(brewer.pal(9, "BrBG"))(20)

pdf("lncs_liver_clustered_heatmap_TIER2_lncs.pdf", pointsize=6)
heatmap.2(as.matrix(tier2_lncs[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          col=hmcol,           
          margins =c(12,9),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()


#All PCGs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hmcol <- colorRampPalette(brewer.pal(9, "BrBG"))(100)
rv <- rowVars(all[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 
cols <- palette(mypal)[as.fumeric(clin[,12])]
head(cbind(colnames(all[,1:68]),cols))

pdf("all_liver_clustered_heatmap_all_PCGs.pdf", pointsize=6)
heatmap.2(as.matrix(all[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          margins =c(12,9),
          col=hmcol,
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()

#Tier 1 PCGs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tier1 <- all[all$tier=="TIER1",]
rv <- rowVars(tier1[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 

pdf("PCG_liver_clustered_heatmap_TIER1.pdf", pointsize=6)
heatmap.2(as.matrix(tier1[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          col=hmcol,
          margins =c(12,9),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()

#Tier 2 PCGs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tier2 <- all[all$tier=="TIER2",]
rv <- rowVars(tier2[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 
hmcol <- colorRampPalette(brewer.pal(9, "BrBG"))(20)

pdf("PCGs_liver_clustered_heatmap_TIER2.pdf", pointsize=6)
heatmap.2(as.matrix(tier2[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          col=hmcol,           
          margins =c(12,9),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          distfun = function(x) dist(x,method = 'manhattan'),
          scale = c("row"))

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()


#Combine lncs and PCGs - get most variable genes 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
all$lnc_type <- "protein_coding" 
all_genes <- rbind(lnc, all)
#scale 
all_genes[,1:68] <- apply(all_genes[,1:68], 1, function(x) scale(x))
#column bar colours
hmcol <- colorRampPalette(brewer.pal(9, "BrBG"))(100)
cols <- palette(mypal)[as.fumeric(clin[,12])]
head(cbind(colnames(all_genes[,1:68]),cols))
rv <- rowVars(all_genes[,1:68]) #variances 
idx <- order(-rv)[1:100] #100 most variable genes 

pdf("PCGs_and_LNCs_liver_clustered_heatmap.pdf", pointsize=6)
heatmap.2(as.matrix(all_genes[idx,1:68]), labCol=clin[,12],
          trace="none", 
          ColSideColors=cols, 
          col=hmcol,           
          margins =c(12,9),
          )

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("HCC", "cHCC/CC", "ICC"), # category labels
    col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
    )

dev.off()
