#top5_cancers_survival_results_analysis.R

#Karina Isaev
#August 21th, 2017

#Purpose: Analyze further the lncRNAs obtained as signficantly
#associated with survival 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Libraries
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Data
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

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
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 

#Survival results 
result <- fread("results_coxAug28_median5fpkmMin.txt", sep=";")

#---------------------------------------------------------
##Processing
#---------------------------------------------------------

#Write function that takes a dataframe, calculates medians and 
#output list of genes with median greater than that 

check_medians <- function(column){
  med <- median(column)
  if(med >=5){
    return(med)
  } 
}

#save results
high_lncs <- as.data.frame(matrix(ncol=3))
colnames(high_lncs) <- c("median", "gene", "canc")

#apply to dataframe 

for(i in 1:length(unique(lnc_rna$canc))){

#subset RNA-dataset to one cancer type
df <- subset(lnc_rna, lnc_rna$canc %in% unique(lnc_rna$canc)[i])

#apply function
res <- apply(df[,1:5607], 2, check_medians)
res <- Filter(Negate(is.null), res)  
res <- data.frame(median=(matrix(unlist(res), nrow=length(res), byrow=T)), gene = names(res), canc=unique(lnc_rna$canc)[i])
high_lncs <- rbind(high_lncs, res)

}#end loop

high_lncs <- high_lncs[-1,]
write.table(high_lncs, file="high_lncsmed4top5cancersPCAWG.txt", sep=";", quote=F, row.names=F)

#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------

lnc_rna <- lnc_rna[,c((which(colnames(lnc_rna) %in% high_lncs$gene)), 5608,5609)] #215 lncRNAs remain 

#For each patient add survival status and days since last seen 
lnc_rna$status <- ""
lnc_rna$time <- ""
lnc_rna$sex <- ""
lnc_rna$tumour_stage <- ""
lnc_rna$tumour_grade <- ""

#lncs
for(i in 1:nrow(lnc_rna)){
  pat <- rownames(lnc_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  lnc_rna$status[i] <- clin$donor_vital_status[z]
  lnc_rna$sex[i] <- clin$donor_sex[z]
  lnc_rna$tumour_stage[i] <- clin$tumour_stage[z]
  lnc_rna$tumour_grade[i] <- clin$tumour_grade[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Analysis
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sig <- filter(result, pval < 0.05)
sig <- as.data.frame(sig)
write.table(sig, "42sig_lncRNACancerAssociations.txt", quote=F, row.names=F, sep=";")

#Plot expression distribution of high expressing lncRNAs versus significant survival associated ones
#and compare to GTEX pattern 
fdr_sig <- sig[sig$fdr <0.1,]

lnc_rna[,1:215] <- log1p(lnc_rna[,1:215])

cancers_list <- unique(high_lncs$canc)
#only care about pancreas, ovary, liver and clear cell 
cancers_list <- cancers_list[c(1,4,5,6)]

#Function 1
#input: list of cancers who had high expressing genes 
#output: list of cancer-genes associated with it because they had high expression
get_genes <- function(cancer){
	genes <- high_lncs$gene[high_lncs$canc==cancer]
	return(list(cancer, genes))
}

canc_genes_list <- lapply(cancers_list, get_genes)

#Function 2
#input: list of cancer-genes associated with it 
#output: list of dataframes with only samples associated with that tissue and rna expression of only those genes 

get_rna_data <- function(canc_genes){
	tissue <- canc_genes[[1]]
	genes <- canc_genes[[2]]
	z2 <- which(lnc_rna$canc == tissue)  
	z <- which(colnames(lnc_rna) %in% genes)
	df <- lnc_rna[z2,c(z, ncol(lnc_rna)-6)]
	return(df)
}

canc_gene_expression <- lapply(canc_genes_list, get_rna_data)


#Function 3
#input: list of datframes
#output: dataframe prepared for plotting violin plots side by side for each gene 

for(i in 1:length(canc_gene_expression)){
	dataframe <- canc_gene_expression[[i]]
	#nrow = #genes * #patients
	rows <- (nrow(dataframe) * (ncol(dataframe)-1))
	df3 <- as.data.frame(matrix(nrow=rows, ncol=4))
	colnames(df3) <- c("gene", "patient", "exp", "PCAWG")

	genes <- colnames(dataframe)[1:(ncol(dataframe)-1)]
	
	meds <- list()
	sig_genes <- sig$gene[which(sig$canc %in% dataframe$canc[1])]
	fdr_genes <- fdr_sig$gene[which(fdr_sig$canc %in% dataframe$canc[1])]

	for(j in 1:length(genes)){
		gene <- genes[j]
		z <- which(is.na(df3[["gene"]]))[1]       
		gene_exp <- which(colnames(dataframe) %in% gene)
		gene_exp <- dataframe[,gene_exp]
		h <- length(gene_exp)
		df3[z:(z+h-1),1] <- gene
		df3[z:(z+h-1),2] <- rownames(dataframe) 
		df3[z:(z+h-1),3] <- gene_exp
		sig_tag <- length(which(sig_genes %in% gene))
		fdr_tag <- length(which(fdr_genes %in% gene))
		if(sig_tag == 0){
			sig_tag <- "Not Sig"
		}
		if(sig_tag == 1){
			sig_tag <- "Sig"
		}
		if(fdr_tag==1){
			sig_tag <- "FDR Sig"
		}
		df3[z:(z+h-1),4] <- sig_tag
		med <- median(gene_exp)
		meds[[j]] <- c(gene, med)
	}
	#get list of ordered genes
	meds <- as.data.table(matrix(unlist(meds), ncol=2, byrow=T))	
	meds$V2 <- as.numeric(meds$V2)
	meds <- meds[order(V2)]
	order <- meds$V1
	tissue <- dataframe$canc[1]
	pdf(paste(tissue, "high_lncs_expression.pdf", sep="_"), pointsize=9, width=17, height=13)
	g <- ggboxplot(df3, x= "gene", y="exp", add="boxplot", order=order, palette=mypal, fill="PCAWG")
	g <- ggpar(g, main= paste(tissue, "High Expressing lncRNAs in Cancer,", nrow(dataframe), "PCAWG Samples") ,xlab = "lncRNA", ylab = "log1p(FPKM)", xtickslab.rt = 65, font.tickslab = c(7,"plain", "black"))
	g <- g + geom_hline(aes(yintercept=log1p(4.5)), colour="#990000")
	print(g) 
	dev.off()
}





















#[1]
#For each cancer type, assess whether the lncRNAs significant 
#in that cancer are co-expressed with each other 

#write function divide sig df into cancer type
#using genes in that cancer type, run co-expression 
#so really just need function that taked n number of unique genes
#finds them in dataframe and runs lm 

#Function 1 
#input: cancer type 
#output: vector cancer type + list of genes 

cancers_list <- unique(sig$canc)
get_genes <- function(cancer){
	genes <- sig$gene[sig$canc==cancer]
	return(list(cancer, genes))
}

canc_genes_list <- lapply(cancers_list, get_genes)

#Function 2
#input: cancer type and list of genes
#output: cancer with list of median tag corresponding to patients 

gene_tag <- function(column){
  median <- median(as.numeric(column))
  #median2 <- median(df[,1])
  res <- c()
  for(y in 1:length(column)){
    genexp <- column[y]
    if(genexp >= median){
      res <- c(res, 1)
      }
    if(genexp < median){
      res <- c(res, 0)
      }
    }
    return(res)   
}

#Function 3
#input: cancer type and list of genes
#output: list of dataframes corresponding to each cancer (4 dataframes in total)
#with gene tags for each gene in that cancer 
multivariate_prep <- function(canc_genes){
	canc <- canc_genes[[1]]
	genes <- canc_genes[[2]]
	df <- lnc_rna[lnc_rna$canc==canc,]
	#add high or low tag to each genes 
	df <- df[,c(which(colnames(df) %in% genes),245:251)]
	#change the values in each gene's column to median tag values
	j <- ncol(df)-7
	tags <- apply(df[,1:j], 2, gene_tag)
	rownames(tags) <- rownames(df)
	j <- ncol(df)-6
	j2 <- ncol(df)
	df <- cbind(df[,j:j2], tags)
	return(df)
}

#apply to canc_genes_list
genes_to_test <- lapply(canc_genes_list, multivariate_prep) #list of dataframes 

#Function 4
#input: dataframe with 0/1 indicating high or low genes
#output: patients given a cluster number to be used further in survival analysis 

make_cluster <- function(dataframe){
	 km1 = kmeans(dataframe[,8:dim(dataframe)[2]], 2, nstart=100)
	 df <- dataframe
	 df$cluster <- km1$cluster
	 return(df)
}

clustered_patients <- lapply(genes_to_test, make_cluster) #list of dataframes 

#Function 5
#input: dataframe with 0/1 indicating high or low genes
#output: heatmap showing how genes tags are related 

make_heatmap <- function(dataframe){
	myColors <- brewer.pal(5,"Set1")[1:2]
	names(myColors) <- levels(as.factor(dataframe$cluster))
	dataframe$colour <- ""
	dataframe$colour[dataframe$cluster==1] <- myColors[1]
	dataframe$colour[dataframe$cluster==2] <- myColors[2]
	j <- dim(dataframe)[2] - 2
	
	heatmap.2(
	as.matrix(dataframe[,8:j]),
	main = paste(dataframe$canc[1], "lncRNA Tags"), # heat map title
	RowSideColors=dataframe$colour,
	cexRow = 0.75, cexCol=0.75,
  	density.info="none",  # turns off density plot inside color legend
  	trace="none", 
  	keysize =0.75,
  	#margins =c(12,9),       # turns off trace lines inside the heat map
  	col=mypal[c(10,1)], key.title=NA, key.xlab=NA) 

	par(lend = 1)           # square line ends for the color legend
	legend("topright",      # location of the legend on the heatmap plot
    legend = c("cluster1", "cluster2"), # category labels
    col = unique(dataframe$colour)[c(2,1)], lwd=1, cex=0.6)

}

pdf("canc_lncRNA_candidate_TAgs_heatmaps_clustered.pdf", pointsize=10, width=13, height=11)
lapply(clustered_patients, make_heatmap)
dev.off()


#Function 6
#input: dataframe wity clustered patients 
#output: univariate survival analysis based on cluster number 

uni_cluster_survival <- function(dataframe){
	#cox
    dataframe$status[dataframe$status=="alive"] <- 0
    dataframe$status[dataframe$status=="deceased"] <- 1
    dataframe$status <- as.numeric(dataframe$status)
    dataframe$time <- as.numeric(dataframe$time)
    dataframe$cluster <- as.numeric(dataframe$cluster)
  	#cox regression 
    res.cox <- coxph(Surv(time, status) ~ cluster, data = dataframe)
    row <- c(dataframe$canc[1], summary(res.cox)$coefficients[1,c(1,2,5)])
    return(row)
}

clustered_cox_results <- lapply(clustered_patients, uni_cluster_survival) 

#Function 7
#input: dataframe wity clustered patients 
#output: kaplan mier plots with log rank pvalue 

uni_cluster_plot <- function(dataframe){
  #cox
  dataframe$status[dataframe$status=="alive"] <- 0
  dataframe$status[dataframe$status=="deceased"] <- 1
  dataframe$status <- as.numeric(dataframe$status)
  dataframe$time <- as.numeric(dataframe$time)
      
          #plot survival plot
          fit <- survfit(Surv(time, status) ~ cluster, data = dataframe)
          s <- ggsurvplot(
          fit, 
          main = dataframe$canc[1],       
          legend.labs = c(paste("Cluster1", dataframe$canc[1]), paste("Cluster2", dataframe$canc[1])),             # survfit object with calculated statistics.
          data = dataframe,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,2000),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 500,     # break X axis in time intervals by 500.
          palette = colorRampPalette(mypal)(14), 
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s) 
}

pdf("clustered_patients_usingkmeans_survival_curves.pdf", pointsize=8, width=14, height=12)
lapply(clustered_patients, uni_cluster_plot)
dev.off()

#Function 8
#input: list of dataframes 
#output: multivaraite cox hazards ratio results 
multivariate <- function(dataframe){
  		#cox
        dataframe$status[dataframe$status=="alive"] <- 0
        dataframe$status[dataframe$status=="deceased"] <- 1
        dataframe$status <- as.numeric(dataframe$status)
        dataframe$time <- as.numeric(dataframe$time)
      	
      	covariates <- colnames(dataframe)[8:dim(dataframe)[2]]
      	covariates <- sub("-", ".", covariates, fixed=TRUE)
		univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status)~', x)))

		#univariate 
		colnames(dataframe)[8:dim(dataframe)[2]] <- sub("-", ".", colnames(dataframe)[8:dim(dataframe)[2]], fixed=TRUE)
		univ_models <- lapply(univ_formulas, function(x){coxph(x, data = dataframe)})

		# Extract data 
		univ_results <- lapply(univ_models,
                       	  function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                         })
		res <- t(as.data.frame(univ_results, check.names = FALSE))
		res <- as.data.frame(res)
		res$patient <- rownames(res)
		res$p.value <- as.numeric(res$p.value)
		res$fdr <- p.adjust(res$p.value, method="bonferroni")
		res <- as.data.table(res)
		res <- res[order(fdr)]
		
		pdf(paste(dataframe$canc[1], "univaraite_comparisons.pdf", sep="_"), pointsize=9, height=11)		
		grid.table(res)
		dev.off()

        }

lapply(clustered_patients, multivariate)

















