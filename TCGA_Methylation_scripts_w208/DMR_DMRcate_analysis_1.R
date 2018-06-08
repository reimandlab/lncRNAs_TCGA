options(stringsAsFactors=F)
#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
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

# load packages required for analysis
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)


mypal = pal_npg("nrc", alpha = 0.7)(10)

#MAIN PACKAGE
library(DMRcate)

#all Methylation files - already processed (NAs and 0 probes removed)
#file converted into matrix 

list.files("files_generated_by_part2_transform_firehose_methylation_files_23cancers_june6/", pattern="methylation_matrix.rds")


##--------Methylation File----------------------------------------------------------------------------------------------

MESO = readRDS("MESO_methylation_matrix.rds")

#Before doing DMR analysis need to assign each lncRNA candidate high risk & low risk patients 

##--------RNA expressionfile--------------------------------------------------------------------------------------------

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

##--------lncRNA prognostic candidates----------------------------------------------------------------------------------

#Data---
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
cands = filter(cands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, 166 unique lncRNAs, 23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
cands$Combo = NULL
cands = cands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
cands = cands[!duplicated(cands), ]
cands_dups = unique(cands$gene[which(duplicated(cands$gene))])

cands$canc[cands$canc=="Ovarian serous cystadenocarcinoma"] = "OV"
cands$canc[cands$canc=="Liver hepatocellular carcinoma"] = "LIHC"
cands$canc[cands$canc=="Kidney renal clear cell carcinoma"] = "KIRC"
cands$canc[cands$canc=="Lung squamous cell carcinoma"] = "LUSC"
cands$canc[cands$canc=="Breast invasive carcinoma"] = "BRCA"
cands$canc[cands$canc=="Lung adenocarcinoma"] = "LUAD"
cands$canc[cands$canc=="Pancreatic adenocarcinoma"] = "PAAD"
cands$canc[cands$canc=="Adrenocortical carcinoma"] = "ACC"
cands$canc[cands$canc=="Bladder Urothelial Carcinoma"] = "BLCA"
cands$canc[cands$canc=="Stomach adenocarcinoma"] = "STAD"
cands$canc[cands$canc=="Head and Neck squamous cell carcinoma"] = "HNSC"
cands$canc[cands$canc=="Brain Lower Grade Glioma"] = "LGG"
cands$canc[cands$canc=="Sarcoma"] = "SARC"
cands$canc[cands$canc=="Kidney renal papillary cell carcinoma"] = "KIRP"
cands$canc[cands$canc=="Mesothelioma"] = "MESO"
cands$canc[cands$canc=="Uterine Corpus Endometrial Carcinoma"] = "UCEC"
cands$canc[cands$canc=="Uveal Melanoma"] = "UVM"
cands$canc[cands$canc=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] = "CESC"
cands$canc[cands$canc=="Colon adenocarcinoma"] = "COAD"
cands$canc[cands$canc=="Rectum adenocarcinoma"] = "READ"
cands$canc[cands$canc=="Thyroid carcinoma"] = "THCA"
cands$canc[cands$canc=="Glioblastoma multiforme"] = "GBM"
cands$canc[cands$canc=="Esophageal carcinoma"] = "ESCA"

#Subset candidates to cancer type in methylation file
rna <- rna[rna$type=="MESO",]
cancer = rna$type[1]


##--------Analaysis start--------------------------------------------------------------------------------------------------

#1. label patients for each lncRNA candidate as high expression/low expression
#2. run DMR between high and low expression patients
#3. save coordinates 


##--------lncRNA prognostic candidates label--------------------------------------------------------------------------------

z = which((cands$canc == cancer))
lncs = as.list((as.character(cands$gene[z])))

get_lnc_cor = function(lnc){
		canc = rna$type[1]
		z = which(colnames(rna) == lnc)
		lnc_data = rna[,c(z, 1, 5922)]
		#add lnc risk group 
		median2 = median(as.numeric(lnc_data[,1]))
		lnc_data$lncRNA_exp = ""
		if(median2 == 0){
			lnc_data$lncRNA_exp[lnc_data[,1] == 0] = "Low"
			lnc_data$lncRNA_exp[lnc_data[,1] > 0] = "High"
		}
		if(!(median2 == 0)){
			lnc_data$lncRNA_exp[lnc_data[,1] < median2] = "Low"
			lnc_data$lncRNA_exp[lnc_data[,1] >= median2] = "High"
		}

		z = which((cands$canc == cancer) & (cands$gene == lnc))
		HR = as.numeric(cands$HR[z])

		if(HR > 1){
			risk = "High"
		}

		if(HR < 1){
			risk = "Low"
		}

		lnc_data$risk_type = risk

		return(lnc_data)		

		}#end get_lnc_cor


lnc_data = llply(lncs, get_lnc_cor, .progress="text")


##--------lncRNA prognostic candidates label--------------------------------------------------------------------------------

#1. invert patients = columns , probes=rows
meth_matrix = t(MESO)

print("surccess transposing matrix")

#2. convert NAs to 0s to can feed to the following functions
meth_matrix[is.na(meth_matrix)] <- 0
print(is.matrix(meth_matrix))
mybetas = apply(meth_matrix, 2, as.numeric)

print("surccess matrix now numeric")

#3. convert to M values
myMs = logit2(mybetas)
print(is.matrix(mybetas))
rownames(mybetas) = rownames(meth_matrix)

rownames(myMs) = rownames(meth_matrix)

print("surccess M values")
print(paste(nrow(myMs), "probes prior to SNP filtering"))

#4. remove probes that are 0-2 nucleotides away from snp with MAF > 0.05

myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
print(paste(nrow(myMs.noSNPs), "probes after SNP filtering"))


##--------run DMR on each candidate------------------------------------------------------------------------------------------

#First subset to only patients that have both RNA-Seq data and metylation data
z = which(colnames(myMs.noSNPs) %in% lnc_data[[1]]$patient)
myMs.noSNPs = myMs.noSNPs[,z]

get_DMR_lnc = function(dtt){

	#need to somehow tell DMR which patients have high expression and whihc patients have low expression
	#change columns of methylation matrix so it's like
	#patientID_highexp or patientID_lowexp

	#write function that changes column names
	
	dim(dtt)	

	change_col = function(colname){
		z = which(dtt$patient == colname)
		new_pat = paste(colname, dtt$lncRNA_exp[z], sep="_")
		return(new_pat)
	}

	lnc_meth = myMs.noSNPs
	colnames(lnc_meth) = llply(colnames(lnc_meth), change_col, .progress="text")

	patient <- factor(sub("_.*", "", colnames(lnc_meth)))
	type <- (sub(".*_", "", colnames(lnc_meth)))
	design <- model.matrix(~type)

	myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=2)

	dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
	results.ranges <- extractRanges(dmrcoutput, genome = "hg19")

	groups <- c(High="magenta", Low="forestgreen")
	cols <- groups[as.character(type)]

	samps <- c(1:2)


	 DMR.plot(ranges=results.ranges, dmr=1, CpGs=mybetas, what="Beta", arraytype = "450K",
         phen.col=cols, genome="hg19", samps=samps)

}


llply(lnc_data, get_DMR_lnc, .progress="text")









































