#42candidate_lncs_limma_diffCoexpression .R

#Karina Isaev
#September 5th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run co-expression analysis with cancer type and sex? as confounders
#to get list of co-expressed PCGs with NEAT1 to run through m:Explorer  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)


#Data#-------------------------------------------------

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

#split by cancer type 

z = which(str_detect(colnames(pcg), "ENSG"))	
pcg = as.data.frame(pcg)
pcg[,z] <- log1p(pcg[,z])

#2. Get lncRNA - median within each tissue type
tissues <- unique(pcg$type)
tissues = tissues[which(tissues %in% c("LGG", "GBM"))]

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------
#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- pcg[pcg$type==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  

#genes to test

get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}

genes_test = c("CATSPER1", "SCN9A", "AQP9", "KCNN4")
genes_test = unlist(llply(genes_test, get_ensg_pcg))

getPCGS <- function(df){
	#add RCC1 state column

	z = which(colnames(df) %in% c("patient", genes_test[[4]]))
	rcc1 = df[,z]
	med = median(as.numeric(rcc1[,2]))
	z = which(rcc1[,2] >= med)
	rcc1$rcc1_tag = ""
	rcc1$rcc1_tag[z] = "high"
	rcc1$rcc1_tag[-z] = "low"

	#merge back with orgiinal dataframe
	df = merge(df, rcc1, by=c("patient", genes_test[[4]]))

	#Remove PCGs with median E < 5 FPKM 	
	#get medians of all PCGs
	z = which(str_detect(colnames(df), "ENSG"))	

	meds <- apply(df[,z], 2, median)

	#names of pcgs with median <5 
	low <- names(meds[which(meds <(log1p(5)))]) 
	df <- df[,-(which(colnames(df) %in% low))] 

	return(df)
}

dividedWpcgs <- llply(tissues_data, getPCGS, .progress = "text")
print("pass2")

#FUNCTION3 - limma differential expression between high and low lncRNA groups
diffE <- function(d){
	design <- model.matrix(~ 0 + factor(d$rcc1_tag))
	colnames(design) <- c("high", "low")
	rownames(d) <- d$patient
	#remove rcc1 from matrix 

	z = which(colnames(d) == genes_test[[4]])
	d = d[,-z]

	z = which(str_detect(colnames(d), "ENSG"))	
	expression <- t(d[,z])

	fit <- lmFit(expression, design)
	cont.matrix <- makeContrasts(LowvsHigh=low-high, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	ps <- fit2$p.value
	ps <- p.adjust(ps, method="fdr")
	#numGenes <- length(which(ps <= 0.05))

	numGenes <- length(ps)

	genes=rownames(expression)
    t <- topTable(fit2,coef=1,adjust.method="fdr",n=numGenes,p.value=1,genelist=genes)

    if(dim(t)[1] > 10){

    #rank list of genes before making heatmap
    t <- as.data.table(t)
    #first by adj p values then by decreasing FC
    #t <- t[order(adj.P.Val)]
    t <- t[order(-abs(as.numeric(logFC)))]
    #t = filter(t, abs(logFC) >=2)
    if(dim(t)[1] >= 5){

    #save top gene names 
    top <- c(paste(genes_test[[4]], d$type[1]), t$ID)
    
    print(paste("t before subset", dim(t)[1]))
    
    #generate volcano plot
    point <- quantile(as.numeric(-log10(t$P.Value)),0.95)
   
	}

	t = as.data.frame(t) 
	t$cancer = d$type[1]

    if(dim(t)[1] <= 1){
    	top <- c(paste(colnames(d)[6], d$canc[1]), "none")
   	}
    return(t)   	
}
}


diffEresults = llply(dividedWpcgs, diffE, .progress="text")
diffEresults1 = ldply(diffEresults, data.frame)

saveRDS(diffEresults1, file=paste(genes_test[[4]], "diff_exp_results.rds", sep="_"))


#FC > 1 means higher in low RCC1 group 


















