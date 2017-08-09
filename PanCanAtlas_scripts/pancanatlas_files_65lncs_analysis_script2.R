#---------------------------------------------------------------------------
#pancanatlas_files_65lncs_analysis_script2.R
#---------------------------------------------------------------------------

#Data: August 8th

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 
#[CODE] 
#[1] identifies cancer specific and nonspecific lncRNAs 
#[2] PCA to cluster patients by lncRNA expression 

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

#---------------------------------------------------------
#Analysis 
#---------------------------------------------------------

f <- fread("all_lncs_cancers.txt", data.table=F, sep="_")
#remove individual file headers
z <- which(f[,1] %in% "gene")
f <- f[-z,]

#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1]))
counts <- counts[order(N)]
#7/49 lncRNAs appear only once
#25/49 lncRNAs appear in all cancer types 

all <- counts$V1[counts$N ==7]
all <- as.data.frame(all)

#unique to a single cancer type
specific <- counts$V1[counts$N ==1]
specific <- as.data.frame(specific)


#make violin plot for 25 lncRNAs that are highly expressed in all cancers
#num of rows = num_genes(25) * patients(2546) = 63650

all_genes <- as.data.frame(matrix(nrow=63650, ncol=4))
colnames(all_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(all)){
	gene2 <- all[i,1]
	if(i == 1){
		all_genes[1:2546,1] <- gene2
		z <- which(colnames(exp_fantom_lncs) %in% gene2)
		dat <- exp_fantom_lncs[,c(z,66)]	
		all_genes[1:2546,2] <- dat$canc
		all_genes[1:2546,3] <- rownames(dat)
		all_genes[1:2546,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(all_genes[,1]))[1]
		all_genes[z:(z+2545), 1] <- gene2
		z2 <- which(colnames(exp_fantom_lncs) %in% gene2)
		dat <- exp_fantom_lncs[,c(z2,66)]	
		all_genes[z:(z+2545),2] <- dat$canc
		all_genes[z:(z+2545),3] <- rownames(dat)
		all_genes[z:(z+2545),4] <- dat[,1]
	}
}

##Plot violin
all_genes_logged <- all_genes  
all_genes_logged$GeneE <- log1p(all_genes_logged$GeneE)

pdf("25_high_all_tissues_lncs_alsoFantom.pdf", pointsize=5, width=16, height=14)
g <- ggviolin(all_genes_logged, x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene")
g <- ggpar(g, font.legend = c(8, "plain", "black")) 
g <- g + labs(title = "25 lncRNAs Highly Expressed in all Cancers", y="log1p(RPKM)") + 
     theme(plot.title = element_text(hjust = 0.5))
g <- g + rremove("x.text")
g <- g + rremove("x.ticks")     
g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
g
dev.off()

#---------------------------------------------------------
#PCA Analysis - commonly expressed lncRNAs 
#---------------------------------------------------------

#PCA using logged values of 49 lncRNAs highly expressed in at least one cancer type
df <- which(colnames(exp_fantom_lncs) %in% unique(f$gene))
df <- exp_fantom_lncs[,df]
df <- log1p(df)
pdf("pca_using49_pancanatlas_toplncsLogged.pdf", pointsize=4, height=11, width=10)
fviz_cluster(list(data = df, cluster = as.factor(exp_fantom_lncs$canc)), geom = "point", palette=mypal, ggtheme = theme_minimal(), ellipse.type="norm")
dev.off()

#PCA using logged values of 49 lncRNAs highly expressed in at least one cancer type
df <- which(colnames(exp_fantom_lncs) %in% all[,1])
df <- exp_fantom_lncs[,df]
df <- log1p(df)
pdf("pca_using25_pancanatlas_toplncsLogged.pdf", pointsize=4, height=11, width=10)
fviz_cluster(list(data = df, cluster = as.factor(exp_fantom_lncs$canc)), geom = "point", palette=mypal, ggtheme = theme_minimal(), ellipse.type="norm")
dev.off()

#---------------------------------------------------------
#Analysis - Cancer Specific lncRNAs 
#---------------------------------------------------------

#make violin plot for 7 lncRNAs that are highly expressed in specific cancers to show it 

#num of rows = num_genes(7) * patients(2546) = 17822
specific_genes <- as.data.frame(matrix(nrow=17822, ncol=4))
colnames(specific_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(specific)){
	gene2 <- specific[i,1]
	if(i == 1){
		specific_genes[1:2546,1] <- gene2
		z <- which(colnames(exp_fantom_lncs) %in% gene2)
		dat <- exp_fantom_lncs[,c(z,66)]	
		specific_genes[1:2546,2] <- dat$canc
		specific_genes[1:2546,3] <- rownames(dat)
		specific_genes[1:2546,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(specific_genes[,1]))[1]
		specific_genes[z:(z+2545), 1] <- gene2
		z2 <- which(colnames(exp_fantom_lncs) %in% gene2)
		dat <- exp_fantom_lncs[,c(z2,66)]	
		specific_genes[z:(z+2545),2] <- dat$canc
		specific_genes[z:(z+2545),3] <- rownames(dat)
		specific_genes[z:(z+2545),4] <- dat[,1]
	}
}

##Plot violin
specific_genes_logged <- specific_genes  
specific_genes_logged$GeneE <- log1p(specific_genes_logged$GeneE)
specific_genes_logged$Gene <- as.factor(specific_genes_logged$Gene)

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

#plot indiviaul violin plots for each of 9 genes coloured by cancer type 
pdf("9_high_specific_tissues_lncs.pdf", pointsize=8, width=14, height=13)

for (i in seq(1, length(unique(specific_genes_logged$Gene)), 4)) {
    
    g <- ggviolin(specific_genes_logged[specific_genes_logged$Gene %in% levels(specific_genes_logged$Gene)[i:(i+3)], ], 
                  x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene",
                  order = unique(specific_genes_logged$Cancer)[c(1,3,7,2,4,5,6)])
   
    g <- ggpar(g, font.legend = c(6, "plain", "black")) 
	g <- g + labs(title = "9 lncRNAs Highly Expressed in Specific Cancers", y="log1p(RPKM)") + 
     	 theme(plot.title = element_text(hjust = 0.5))
	g <- g + rremove("x.text")
	g <- g + rremove("x.ticks")  
	g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
	print(g)
}

dev.off()

