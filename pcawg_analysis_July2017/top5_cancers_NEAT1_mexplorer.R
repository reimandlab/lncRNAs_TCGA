#top5_cancers_NEAT1_mexplorer.R

#Karina Isaev
#August 11th, 2017

#Purpose: compare list of PCGs found through co-expression
#with Neat1 versus that obtained from binding site intersection

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
library(GeneOverlap)
library(VennDiagram)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

exp <- fread("3132_pcgs_coexpressedWneat1.txt", data.table=F)
load("neat1_bound_genes.rsav")
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]
#background would be the list of protein coding genes 
back <- ucsc[ucsc$hg19.ensemblSource.source == "protein_coding",]

#Processing-------------------------------------------------

neat1_bound_genes <- neat1_bound_genes[-(which(duplicated(neat1_bound_genes$gene_symbol))), ]
colnames(exp)[1] <- "gene_symbol"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Using both positive and negative co-expression relationships 

#Fisher's test----------------------------------------------
#GeneOverlap
obj <- newGeneOverlap(neat1_bound_genes$gene_symbol, exp$gene_symbol, genome.size=21225) #21225 unique PCGs in UCSC that both Drew and I started with but not all these were actually used in co-expression s
test <- testGeneOverlap(obj)

#Venn Diagram-----------------------------------------------
pdf("venn_all_coexpressed_bounds224incommon.pdf", pointsize=10)
draw.pairwise.venn(area1=2239, area2=900, cross.area=224, category=c("Bound", "Co-Expressed"), lty=rep("blank", 2), 
	fill = c("light blue", "pink"), alpha=rep(0.5, 2), cat.pos=c(0, 0), cat.dist=rep(0.025,2))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Using just positive co-expression relationships 
exp_pos <- exp[exp$cor == "pos", ]

#Fisher's test----------------------------------------------
#GeneOverlap
obj <- newGeneOverlap(neat1_bound_genes$gene_symbol, exp_pos$gene_symbol, genome.size=21225) #21225 unique PCGs in UCSC that both Drew and I started with but not all these were actually used in co-expression s
test <- testGeneOverlap(obj)

#Venn Diagram-----------------------------------------------
pdf("venn_pos_coexpressed_bounds73incommon.pdf", pointsize=10)
draw.pairwise.venn(area1=2239, area2=414, cross.area=73, category=c("Bound", "Positively Co-Expressed"), lty=rep("blank", 2), 
	fill = c("light blue", "pink"), alpha=rep(0.5, 2), cat.pos=c(0, 0), cat.dist=rep(0.025,2))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Using just negative co-expression relationships 
exp_neg <- exp[exp$cor == "neg", ]

#Fisher's test----------------------------------------------
#GeneOverlap
obj <- newGeneOverlap(neat1_bound_genes$gene_symbol, exp_neg$gene_symbol, genome.size=21225) #21225 unique PCGs in UCSC that both Drew and I started with but not all these were actually used in co-expression s
test <- testGeneOverlap(obj)

#Venn Diagram-----------------------------------------------
pdf("venn_neg_coexpressed_bounds151incommon.pdf", pointsize=10)
draw.pairwise.venn(area1=2239, area2=486, cross.area=151, category=c("Bound", "Negatively Co-Expressed"), lty=rep("blank", 2), 
	fill = c("light blue", "pink"), alpha=rep(0.5, 2), cat.pos=c(0, 0), cat.dist=rep(0.025,2))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Enrichment Map Setup---------------------------------------
#Both positive and negative coexpressed genes 
load("mExplorer_analysis_additional_tags_result.rsav")
paths <- result$resp_coefs
paths <- as.data.frame(paths)
paths$Pathway <- rownames(paths)
#need to covnvert "." to ":"
paths$Pathway <- gsub("\\.", ":", paths$Pathway)

#GMT file 
gmt <- fread("pathways_names.gmt", data.table=F, header=F)
colnames(gmt) <- c("Pathway", "Name")

#Output needs to look like
#GO.ID      {tab} Description                     {tab} p.Val
#GO:0000346 {tab} transcription export complex    {tab} 0.01
#GO:0030904 {tab} retromer complex                {tab} 0.05
#GO:0008623 {tab} chromatin accessibility complex {tab} 0.05
#GO:0046540 {tab} tri-snRNP complex               {tab} 0.01

paths <- merge(paths, gmt, by="Pathway")
#now need to add score value
scores <- result$scores
scores <- as.data.frame(scores)
scores$Pathway <- rownames(scores)
scores$Pathway <- gsub("\\.", ":", scores$Pathway)
paths <- merge(paths, scores, by="Pathway") #136 

#convert format to match required format above
paths_new <- paths[,c(1,7,8)]
colnames(paths_new)[1] <- "GO.ID"
colnames(paths_new)[2] <- "Description"
colnames(paths_new)[3] <- "p.Val"

write.table(paths_new, sep="\t", file="enrich_results_file.txt", quote=F, row.names=F)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Enrichment Map Setup---------------------------------------
#Just negative coexpressed genes 
load("mExplorer_analysis_neg_only_result.rsav")
paths <- result$resp_coefs
paths <- as.data.frame(paths)
paths$Pathway <- rownames(paths)
#need to covnvert "." to ":"
paths$Pathway <- gsub("\\.", ":", paths$Pathway)

#GMT file 
gmt <- fread("pathways_names.gmt", data.table=F, header=F)
colnames(gmt) <- c("Pathway", "Name")

#Output needs to look like
#GO.ID      {tab} Description                     {tab} p.Val
#GO:0000346 {tab} transcription export complex    {tab} 0.01
#GO:0030904 {tab} retromer complex                {tab} 0.05
#GO:0008623 {tab} chromatin accessibility complex {tab} 0.05
#GO:0046540 {tab} tri-snRNP complex               {tab} 0.01

paths <- merge(paths, gmt, by="Pathway")
#now need to add score value
scores <- result$scores
scores <- as.data.frame(scores)
scores$Pathway <- rownames(scores)
scores$Pathway <- gsub("\\.", ":", scores$Pathway)
paths <- merge(paths, scores, by="Pathway") #126 

#convert format to match required format above
paths_new <- paths[,c(1,5,6)]
colnames(paths_new)[1] <- "GO.ID"
colnames(paths_new)[2] <- "Description"
colnames(paths_new)[3] <- "p.Val"

write.table(paths_new, sep= "\t", file="NegativeCoGenesenrich_results_file.txt", quote=F, row.names=F)






















