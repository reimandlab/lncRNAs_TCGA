#---------------------------------------------------------------------------
#pancanatlas_files_65lncs_analysis_script1.R
#---------------------------------------------------------------------------

#Data: August 8th

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 
#[CODE] 
#makes plot with number of lncRNAs above different 
#median cutoffs per cancer type 

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

#RNA-Seq data for patients from the equivalent top5 cancer types obtained 
#using PCAWG + 65 lncRNAs 
exp_fantom_lncs <- readRDS("panAtlas_65lncs_RNA-Seq_data.rds")
exp_fantom_lncs <- as.data.frame(exp_fantom_lncs)

#Patient ids and centres, data n=2814 in RNA-Seq file above 
patients <- readRDS("TOP5_PCAWGcancertypes_tcga_rnaseqfile_patients_cancertypes_conversion.rds")

#Clinical data 
clin <- fread("clinical_PANCANatlas_patient_with_followup.tsv")

#---------------------------------------------------------------------------
#Processing 
#---------------------------------------------------------------------------

#[1] add cancer type data to expression file for easy stratification 
exp_fantom_lncs <- t(exp_fantom_lncs)
exp_fantom_lncs <- as.data.frame(exp_fantom_lncs)
exp_fantom_lncs$canc <- ""

for(i in 1:nrow(exp_fantom_lncs)){
	z <- which(patients$id %in% rownames(exp_fantom_lncs)[i])
	exp_fantom_lncs$canc[i] <- patients$cancer[z]
}

#quick plot of number of patients per cancer type 
pats <- as.data.table(table(exp_fantom_lncs$canc))
pats <- pats[order(N)]
pdf("top5_cancers_Pancanatlas_patient_distribution.pdf")
p <- ggbarplot(pats, "V1", "N", fill="V1", xlab="Cancer Type", palette=mypal) 
p <- ggpar(p,
 font.tickslab = c(6,"bold", "black"),
 xtickslab.rt = 70,
 legend="none", yticks.by = 100)
p
dev.off()


#---------------------------------------------------------
#Analysis - plot medians cutoffs
#---------------------------------------------------------

#For each subtype:
#Generate (1) plot showing how many lncRNAs left after different 
#median cutoffs --> then decide on a cutoff based on how many lncRNAs 
#there are

plots <- list()

for(i in 1:length(unique(exp_fantom_lncs$canc))){
	tis <- unique(exp_fantom_lncs$canc)[i]
	z <- which(exp_fantom_lncs$canc %in% tis)
	tis_exp <- exp_fantom_lncs[z,]
	#measure medians of genes and save them 
	meds <- as.data.frame(matrix(nrow=dim(tis_exp)[2]-1,ncol=2))
	colnames(meds) <- c("Gene", "MedianE")
	meds$Gene <- colnames(tis_exp[,1:dim(tis_exp)[2]-1])
	meds$MedianE <- apply(tis_exp[,1:dim(tis_exp)[2]-1], 2, median)
	meds$check1 <- ""
	meds$check2 <- ""
	meds$check3 <- ""
	meds$check4 <- ""
	meds$check5 <- ""
	meds$check10 <- ""
	
	for(y in 1:nrow(meds)){
		m <- meds$MedianE[y]
		if(m >=10){
			meds[y,3:8] <- 1
		}
		if(m >=5){
			meds[y,3:7] <- 1
		}
		if(m >=4){
			meds[y,3:6] <- 1
		}
		if(m >=3){
			meds[y,3:5] <- 1
		}
		if(m >=2){
			meds[y,3:4] <- 1
		}
		if(m >=1){
			meds[y,3] <- 1
		}
		if(m < 1){
			meds[y,3] <- 0
		}
	}
	#plot how many lncRNAs meet each median cutoff
	plot_meds <- as.data.frame(matrix(nrow=6,ncol=2))
	colnames(plot_meds) <- c("Median", "Number_lncRNAs")
	plot_meds[,1] <- c(1:5,10)
	plot_meds[6,2] <- length(which(meds$check10==1))
	plot_meds[5,2] <- length(which(meds$check5==1 & (!(meds$check10==1))))
	plot_meds[4,2] <- length(which(meds$check4==1 & (!(meds$check5==1))))
	plot_meds[3,2] <- length(which(meds$check3==1 & (!(meds$check4==1))))
	plot_meds[2,2] <- length(which(meds$check2==1 & (!(meds$check3==1))))
	plot_meds[1,2] <- length(which(meds$check1==1 & (!(meds$check2==1))))

	#save plot
	g <- ggbarplot(plot_meds, x="Median", y="Number_lncRNAs", palette=mypal, col="Median", fill="Median", label = TRUE, lab.pos = "in", lab.size = 2.6)
	g <- ggpar(g, legend="none")
	g <- g + labs(title = tis, y="Number of lncRNAs") + 
     theme(plot.title = element_text(hjust = 0.5))
    plots[[i]] <- g

    #save list of high expression lncRNAs 
    name_file <- paste(tis, "list_great5med_lncs.txt")
    save <- as.data.frame(meds[meds$MedianE >=10,1])
    colnames(save)[1] <- "gene"
    save$canc <- tis
    write.table(save, name_file, quote=F, row.names=F, sep="_")

} #end loop

g1 <- plots[[1]]
g2 <- plots[[2]]
g3 <- plots[[3]]
g4 <- plots[[4]]
g5 <- plots[[5]]
g6 <- plots[[6]]
g7 <- plots[[7]]

require(cowplot)
pdf("top5_cancers_lncRNAs_above_Pancanatlas_diffMedians.pdf", pointsize=5, height=13, width=14)
plot_grid(g1,g2,g3,g4,g5,g6,g7, labels = "AUTO", ncol = 2, align = 'v', label_size = 10, scale = 0.9)
dev.off()



