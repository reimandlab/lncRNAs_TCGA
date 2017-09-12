#top5_cancers_extraction_script1.R

#Karina Isaev
#July 24th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script2 - using the top 5 cancer types chosen
#PLOTS:

#1. For each cancer type, identify list of candidate lncRNAs 
#whose expression is specific to that cancer 

#2. Plot their violin plots with increasing medians 

#3. Plot to compare how cancer specific lncRNAs are expressed 
#in other cancer types 

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
library(cowplot)



mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,] #total unique genes left 32843

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
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

#Cancers to use 
tum_types <- fread("top5_cancers_andHISTO_to_keepJuly20.txt", data.table=F)

#---------------------------------------------------------
#Processing
#---------------------------------------------------------

###"NORMAL SAMPLES"
z <- which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))
norm_pats <- conversion$icgc_donor_id[z]

###"TUMOUR SAMPLES" - PROCESSING RNA FILE 
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z]

for(i in 1:ncol(rna)){
	z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna)[i])
	if(!(length(z)==0)){
		colnames(rna)[i] <- conversion$icgc_donor_id[z]
	}
}

extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1:2], 1, ensg)
z <- which(check==TRUE)
rna <- rna[z,]

#2. Check how many IDs match the lncRNA IDs
#none match while trancript number is present 
#remove ie, ENSG00000201285.1 --> ENSG00000201285
#in both rna file and lncs file

extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract3) ; 
fantom[,1] <- apply(fantom[,1:2], 1, extract3)

#Remove duplicate genes 
z <- which(duplicated(rna[,1]))
genes <- rna[z,1]
z <- which(rna[,1] %in% genes)
rna <- rna[-z,]

#Using UCSC keep only antisense, lincRNA and protein-coding genes
z <- which(rna[,1] %in% ucsc$hg19.ensGene.name2)
rna <- rna[z,]
 
#rows
rownames(rna) <- rna[,1]
rna <- rna[,-1]

##Keep only patient samples that also have clinical data
z <- which(colnames(rna) %in% clin$icgc_donor_id)
rna <- rna[,z]

#Divide RNA into lnc_RNA and pcg_RNA
z <- which(ucsc$hg19.ensemblSource.source == "protein_coding")
zz <- which(rownames(rna) %in% ucsc$hg19.ensGene.name2[z])

###################
###PCG_RNA#########
###################
pcg_rna <- rna[zz,]


###################
###LNC_RNA#########
###################
lnc_rna <- rna[-zz,]


###SUBSET CLINICAL AND EXPRESSION FILE TO ONLY THE TOP 5 CANCERS/HISTOS

##CLIN
clin_top5 <- subset(clin, clin$icgc_donor_id %in% tum_pats)
z <- which(duplicated(clin_top5$icgc_donor_id))
clin_top5 <- clin_top5[-z,]
clin_top5 <- subset(clin_top5, (clin_top5$histology_tier2 %in% tum_types$V1) & (clin_top5$histology_tier4 %in% tum_types$V2))
clin_top5$combined_tum_histo <- ""
clin_top5$combined_tum_histo <- paste(clin_top5[,15], clin_top5[,17])

#EXPRESSION - pcgs
z <- which(colnames(pcg_rna) %in% clin_top5$icgc_donor_id)
pcg_rna_top5 <- pcg_rna[,z] #20166
sums <- apply(pcg_rna_top5, 1, sum) 
s <- which(sums==0)
z <- which(rownames(pcg_rna_top5) %in% names(s))
pcg_rna_top5 <- pcg_rna_top5[-z,] #20022
pcg_rna_top5 <- t(pcg_rna_top5)
pcg_rna_top5 <- as.data.frame(pcg_rna_top5)
pcg_rna_top5$canc <- ""
#add patient tum type
for(i in 1:nrow(pcg_rna_top5)){
	pat <- rownames(pcg_rna_top5)[i]
	z <- which(clin_top5$icgc_donor_id %in% pat)
	hist <- clin_top5$combined_tum_histo[z]
	pcg_rna_top5$canc[i] <- hist
}


#EXPRESSION - lncs
z <- which(colnames(lnc_rna) %in% clin_top5$icgc_donor_id)
lnc_rna_top5 <- lnc_rna[,z] #12598
sums <- apply(lnc_rna_top5, 1, sum) 
s <- which(sums==0)
z <- which(rownames(lnc_rna_top5) %in% names(s))
lnc_rna_top5 <- lnc_rna_top5[-z,] #12543
lnc_rna_top5 <- t(lnc_rna_top5)
lnc_rna_top5 <- as.data.frame(lnc_rna_top5)
lnc_rna_top5$canc <- ""
#add patient tum type
for(i in 1:nrow(lnc_rna_top5)){
	pat <- rownames(lnc_rna_top5)[i]
	z <- which(clin_top5$icgc_donor_id %in% pat)
	hist <- clin_top5$combined_tum_histo[z]
	lnc_rna_top5$canc[i] <- hist
}


#---------------------------------------------------------
#Analysis - plot medians cutoffs
#---------------------------------------------------------

#For each subtype:
#Generate (1) plot showing how many lncRNAs left after different 
#median cutoffs --> then decide on a cutoff based on how many lncRNAs 
#there are

plots <- list()

for(i in 1:length(unique(lnc_rna_top5$canc))){
	tis <- unique(lnc_rna_top5$canc)[i]
	z <- which(lnc_rna_top5$canc %in% tis)
	tis_exp <- lnc_rna_top5[z,]
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
	plot_meds[5,2] <- length(which(meds$check5==1))
	plot_meds[4,2] <- length(which(meds$check4==1))  
	plot_meds[3,2] <- length(which(meds$check3==1))
	plot_meds[2,2] <- length(which(meds$check2==1)) 
	plot_meds[1,2] <- length(which(meds$check1==1))

	#save plot
	g <- ggbarplot(plot_meds, x="Median", y="Number_lncRNAs", palette=mypal, col=mypal[i], fill=mypal[i], label = TRUE, lab.pos = "in", lab.size = 3.5)
	g <- ggpar(g, legend="none")
	g <- g + labs(title = tis, y="Number of lncRNAs") + 
     theme(plot.title = element_text(hjust = 0.5))
    plots[[i]] <- g

    #save list of high expression lncRNAs 
    name_file <- paste(tis, "list_great5med_lncs.txt")
    save <- as.data.frame(meds[meds$MedianE >=5,1])
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

pdf("top5_cancers_lncRNAs_above_diffMedians.pdf", pointsize=8, height=14, width=14)
plot_grid(g1,g2,g3,g4,g5,g6,g7, labels = "AUTO", ncol = 2, align = 'v', label_size = 10, scale = 0.9)
dev.off()


#---------------------------------------------------------
#Analysis - how many lncRNAs in common? - 21 lncRNAs plot
#---------------------------------------------------------

#all_lncs_cancers.txt (obtain by cat of the 7 files produced in the above for-loop)

f <- fread("all_lncs_cancers.txt", data.table=F, sep="_")
#remove individual file headers
z <- which(f[,1] %in% "gene")
f <- f[-z,]

#subset to include only lncs covered by FANTOM
z <- which(f[,1] %in% fantom[,1])
f <- f[z,]

#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1]))
counts <- counts[order(N)]
#50/86 lncRNAs appear only once
#9/86 lncRNAs appear in all cancer types 
all <- counts$V1[counts$N ==7]
all <- as.data.frame(all)
#gene ids conversion ENSG-Hugo
all$gene <- ""
for(i in 1:nrow(all)){
	g <- all[i,1]
	z <- which(ucsc$hg19.ensGene.name2 %in% g)
	all$gene[i] <- ucsc$hg19.ensemblToGeneName.value[z]
}

#unique to a single cancer type
specific <- counts$V1[counts$N ==1]
specific <- as.data.frame(specific)
#gene ids conversion ENSG-Hugo
specific$gene <- ""
for(i in 1:nrow(specific)){
	g <- specific[i,1]
	z <- which(fantom$CAT_geneID %in% g)
	specific$gene[i] <- fantom$CAT_geneName[z]
}

#make violin plot for 9 lncRNAs that are highly expressed in all cancers

#num of rows = num_genes(9) * patients(497) = 4473

all_genes <- as.data.frame(matrix(nrow=4473, ncol=4))
colnames(all_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(all)){
	gene2 <- all[i,1]
	gene <- all[i,2]
	if(i == 1){
		all_genes[1:497,1] <- gene
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]	
		all_genes[1:497,2] <- dat$canc
		all_genes[1:497,3] <- rownames(dat)
		all_genes[1:497,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(all_genes[,1]))[1]
		all_genes[z:(z+496), 1] <- gene
		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		all_genes[z:(z+496),2] <- dat$canc
		all_genes[z:(z+496),3] <- rownames(dat)
		all_genes[z:(z+496),4] <- dat[,1]
	}
}

##Plot violin
all_genes_logged <- all_genes  
all_genes_logged$GeneE <- log1p(all_genes_logged$GeneE)

pdf("9_high_all_tissues_lncs_alsoFantom.pdf", pointsize=5, width=16, height=14)
g <- ggviolin(all_genes_logged, x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene")
g <- ggpar(g, font.legend = c(8, "plain", "black")) 
g <- g + labs(title = "9 lncRNAs Highly Expressed in all Cancers", y="log1p(FPKM)") + 
     theme(plot.title = element_text(hjust = 0.5))
g <- g + rremove("x.text")
g <- g + rremove("x.ticks")     
g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
g
dev.off()


#---------------------------------------------------------
#Analysis - can the lncRNAs seperate patients based on
#expression into cancer type?
#---------------------------------------------------------

##Task: make heatmap clustering patients at the top
#by cancer type based on lncRNA expression 
#need matrix
#using only 305 genes to start
genes_305 <- as.data.frame(counts[,1])
lncs_305 <- which(colnames(lnc_rna_top5) %in% genes_305[,1])
lncs_305 <- lnc_rna_top5[,c(lncs_305, 12544)]

df <- lncs_305[,1:213] #number of lncs 

desc_stats <- data.frame(
  Min = apply(df, 2, min), # minimum
  Med = apply(df, 2, median), # median
 Mean = apply(df, 2, mean), # mean
  SD = apply(df, 2, sd), # Standard deviation
  Max = apply(df, 2, max) # Maximum
  )

desc_stats <- round(desc_stats, 1)

library(factoextra)
#Observations are represented by points in the plot, using principal components
#PCA using logged values 
df <- lncs_305[,1:213] #number of lncs 
df <- log1p(df)
pdf("pca_using213MEDIANsof5MIN_toplncsLogged.pdf", pointsize=4, height=11, width=10)
fviz_cluster(list(data = df, cluster = as.factor(lncs_305$canc)), geom = "point", palette=mypal, ggtheme = theme_minimal(), ellipse.type="norm")
dev.off()

#---------------------------------------------------------
#Analysis - Cancer Specific lncRNAs 
#---------------------------------------------------------

##clustering using just 50 lncRNAs super specific
#make violin plot for 50 lncRNAs that are highly expressed in specific cancers to show it 

#num of rows = num_genes(50) * patients(497) = 24850
colnames(specific)[2] <- "gene2"

specific_genes <- as.data.frame(matrix(nrow=24850, ncol=4))
colnames(specific_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(specific)){
	gene2 <- specific[i,1]
	gene <- specific[i,2]
	if(i == 1){
		specific_genes[1:497,1] <- gene
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]	
		specific_genes[1:497,2] <- dat$canc
		specific_genes[1:497,3] <- rownames(dat)
		specific_genes[1:497,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(specific_genes[,1]))[1]
		specific_genes[z:(z+496), 1] <- gene
		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		specific_genes[z:(z+496),2] <- dat$canc
		specific_genes[z:(z+496),3] <- rownames(dat)
		specific_genes[z:(z+496),4] <- dat[,1]
	}
}

##Plot violin
specific_genes_logged <- specific_genes  
specific_genes_logged$GeneE <- log1p(specific_genes_logged$GeneE)
specific_genes_logged$Gene <- as.factor(specific_genes_logged$Gene)

pdf("50_high_specific_tissues_lncs.pdf", pointsize=8, width=14, height=13)

# Pairwise comparisons: Specify the comparisons you want
my_comparisons <- list()
k <- 1
for(i in 1:(length(unique(specific_genes_logged$Cancer))-1)){
	t <- unique(specific_genes_logged$Cancer)[i]
	for(j in i:(length(unique(specific_genes_logged$Cancer))-1)){
		add <- c(t, unique(specific_genes_logged$Cancer)[j+1])
		my_comparisons[[k]] <- add 
		k <- k+1
	}
}


#plot indiviaul violin plots for each of 50 genes coloured by cancer type 

for (i in seq(1, length(unique(specific_genes_logged$Gene)), 4)) {
    
    g <- ggviolin(specific_genes_logged[specific_genes_logged$Gene %in% levels(specific_genes_logged$Gene)[i:(i+3)], ], 
                  x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene",
                  order = unique(specific_genes_logged$Cancer)[c(1,3,7,2,4,5,6)])
   
    g <- ggpar(g, font.legend = c(6, "plain", "black")) 
	g <- g + labs(title = "50 lncRNAs Highly Expressed in Specific Cancers", y="log1p(FPKM)") + 
     	 theme(plot.title = element_text(hjust = 0.5))
	g <- g + rremove("x.text")
	g <- g + rremove("x.ticks")  
	g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
	print(g)
}

dev.off()


#---------------------------------------------------------
#Analysis - can the 163 lncRNAs seperate patients based on 
#cancer type?
#---------------------------------------------------------

##Task: make heatmap clustering patients at the top
#by cancer type based on lncRNA expression 
#need matrix
#using only 305 genes to start
lncs_163 <- which(colnames(lnc_rna_top5) %in% specific[,1])
lncs_163 <- lnc_rna_top5[,c(lncs_163, 12544)]

df <- lncs_163[,1:50]
df <- log1p(df)
pdf("pca_using50_toplncsLogged.pdf", pointsize=4, height=11, width=10)
fviz_cluster(list(data = df, cluster = as.factor(lncs_163$canc)), geom = "point", palette=mypal, ggtheme = theme_minimal(), ellipse.type="norm")
dev.off()



##RANDOMIZATION#-------------------------------------------
z <- which(colnames(lnc_rna_top5) %in% specific[,1])
random_lncs <- lnc_rna_top5[,-z]
random_lncs_sample <- random_lncs[,-12381]

random_list <- c()

for(i in 1:500){
s <- sample(1:ncol(random_lncs_sample), 117)
df <- random_lncs_sample[,s]
df <- log1p(df)
#pdf("pca_using163random_lncsLogged.pdf", pointsize=4, height=11, width=10)
f <- fviz_cluster(list(data = df, cluster = as.factor(random_lncs$canc)), geom = "point", palette=mypal, ggtheme = theme_minimal(), ellipse.type="norm")
random_list[[i]] <- f
}

ml<-marrangeGrob(random_list,nrow=2,ncol=2)
ggsave("pca_using_117_random_lncsLogged.pdf",ml)