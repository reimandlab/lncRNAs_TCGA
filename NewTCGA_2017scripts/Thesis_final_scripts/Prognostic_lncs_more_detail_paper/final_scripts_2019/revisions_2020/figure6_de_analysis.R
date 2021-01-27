set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
colnames(allCands)[7] = "Cancer"
allCands = merge(allCands, canc_conv, by="Cancer")
allCands$combo = paste(allCands$gene, allCands$type)

#------------------------------------------------------------------------------
#Take each lncRNA candidate and conduct LIMMA diff exp between lncRNA
#high and low risk groups
#------------------------------------------------------------------------------

library(edgeR)
library(limma)
library(corrplot)

#COSMIC cancer gene census
census = read.csv("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

pcg_counts = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/counts_19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
table(pcg_counts$type)

#add in counts from LAML and SKCM
laml=readRDS("/u/kisaev/TCGA_LAML_count_data/LAML_count_data.rds")
skcm=readRDS("/u/kisaev/TCGA_SKCM_count_data/SKCM_count_data.rds")

#clean up pcg counts
cols = colnames(pcg_counts)[which(str_detect(colnames(pcg_counts), "ENSG"))]
z = which(colnames(pcg_counts) %in% c(cols, "patient", "type"))
pcg_counts = pcg_counts[,z]

#need to bind with new data
laml<-laml[names(pcg_counts)]
skcm<-skcm[names(pcg_counts)]
pcg_counts = rbind(pcg_counts, laml, skcm)

#keep only patients that are the in the main rna matrix to keep same list of patients
z = which(pcg_counts$patient %in% rna$patient)
pcg_counts = pcg_counts[z,]
table(pcg_counts$type)

#------FEATURES-----------------------------------------------------

#1. Get lncRNA - median within each tissue type
allCands$combo = paste(allCands$gene, allCands$Cancer, sep="_")
combos = unique(allCands$combo)

#canc type to cancer conversion
canc_conv = unique(rna[,c("type", "Cancer")])

get_name_pcg = function(pcg){
        z = which(hg38$ensgene == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(hg38$symbol[z])
}

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group

#---------------------------------------------------------
#Function 1
#for each lnc-cancer, label patient as lncRNA-risk or non-risk
#---------------------------------------------------------

#z = which(str_detect(combos, "Brain Lower Grade Glioma"))
#combos = combos[z]
rna=as.data.frame(rna)

get_lnc_canc = function(comb){
	lnc = unlist(strsplit(comb, "_"))[1]
	canc = unlist(strsplit(comb, "_"))[2]
	canc_type = canc_conv$type[which(canc_conv$Cancer == canc)]
	canc_data = subset(pcg_counts, type == canc_type)

	z = which(colnames(rna) %in% c("patient", lnc, "type"))
	lnc_dat = rna[,z]
	z = which(lnc_dat$type == canc_type)
	lnc_dat = lnc_dat[z,]
	med = median(as.numeric(lnc_dat[,which(colnames(lnc_dat)==lnc)]))
	if(med ==0){
		z = which(lnc_dat[,which(colnames(lnc_dat)==lnc)] > 0)
		lnc_dat$lnc_tag = ""
		lnc_dat$lnc_tag[z] = "high"
		lnc_dat$lnc_tag[-z] = "low"
	}

	if(!(med ==0)){
		z = which(lnc_dat[,which(colnames(lnc_dat)==lnc)] >= med)
		lnc_dat$lnc_tag = ""
		lnc_dat$lnc_tag[z] = "high"
		lnc_dat$lnc_tag[-z] = "low"
	}

	#merge back with orgiinal dataframe
	canc_data = merge(canc_data, lnc_dat, by=c("patient", "type"))

	#keep only PCGs not lncRNAs
	pcgs_id = unique(colnames(pcg))
	z = which(colnames(canc_data) %in% c(pcgs_id, colnames(lnc_dat)))
	canc_data = canc_data[,z]

	#Remove PCGs with median E < 5 FPKM
	#get medians of all PCGs
	z = which(str_detect(colnames(canc_data), "ENSG"))

  		#cox ph
  		z = which(allCands$combo == comb)
  		HR = as.numeric(allCands$HR[z])

  		if(HR <1){
  			risk = "low"
  			canc_data$risk = ""
  			canc_data$risk[canc_data$median=="high"] ="noRISK"
  			canc_data$risk[canc_data$median=="low"] ="RISK"
  		}
  		if(HR >1){
  			risk = "high"
  			canc_data$risk = ""
  			canc_data$risk[canc_data$median=="high"] ="RISK"
  			canc_data$risk[canc_data$median=="low"] ="noRISK"
  		}

  		canc_data$lnc = lnc
  		canc_data$canc = canc
  		colnames(canc_data)[which(colnames(canc_data)==lnc)] = "lncRNA"

  	return(canc_data)
	}#end function evaluate_each_lnc

all_canc_lnc_data = llply(combos, get_lnc_canc, .progress="text")

#---------------------------------------------------------
#Function 2
#wtihin each cancer
#calculate for each lncRNAs differentially expressed PCGs
#---------------------------------------------------------

#for lgg add IDH mutation
lgg_dat = readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")

diffE <- function(d){

  print(d$lnc[1])
	if(d$type[1] == "LGG"){
		d = merge(d, lgg_dat, by="patient")
		d$IDH.status = as.character(d$IDH.status)
		z = which(is.na(d$IDH.status))
		d = d[-z,]
		design <- model.matrix(~ 0 + factor(d$lnc_tag) + factor(d$IDH.status))
		colnames(design) <- c("high", "low", "IDH_WT")
	}

	z = which(str_detect(colnames(d), "ENSG"))
	rownames(d) = d$patient

	design <- model.matrix(~ 0 + factor(d$lnc_tag))
	colnames(design) <- c("high", "low")
	rownames(d) <- d$patient

	expression <- t(d[,z])

	# Obtain CPMs
	myCPM <- cpm(expression)
	# Have a look at the output
	head(myCPM)

	# Which values in myCPM are greater than 0.5?
	thresh <- myCPM > 0.5
	# This produces a logical matrix with TRUEs and FALSEs
	head(thresh)

	# Summary of how many TRUEs there are in each row
	table(rowSums(thresh))

	# we would like to keep genes that have at least 10 TRUES in each row of thresh
	keep <- rowSums(thresh) >= 10
	# Subset the rows of countdata to keep the more highly expressed genes
	myCPM = myCPM[keep,]

	counts.keep <- expression[keep,]
	summary(keep)
	dim(counts.keep)

	y <- DGEList(counts.keep)
	# have a look at y
	#y

	# Library size information is stored in the samples slot
	y$samples
	# Get log2 counts per million
	logcounts <- cpm(y,log=TRUE)

	#TMM normalization to eliminate composition biases between libraries
	# Apply normalisation to DGEList object
	y <- calcNormFactors(y)

	#apply voom normalization
	v <- voom(y,design,plot = TRUE)

	# Fit the linear model
	fit <- lmFit(v)
	names(fit)

	cont.matrix <- makeContrasts(LowvsHigh=high-low, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	ps <- fit2$p.value
	ps <- p.adjust(ps, method="fdr")
	#numGenes <- length(which(ps <= 0.05))
  t <- topTable(fit2, coef=1, n="Inf")
  t$ID = rownames(t)
  t$cancer = d$type[1]

	#dist of fold change
	summary(t$logFC)
	t$gene_name = llply(t$ID, get_name_pcg)
	t$lnc = d$lnc[1]
	t$gene_name = as.character(t$gene_name)
    #plot fold changes
    #p <- ggplot(t, aes(logFC, -log10(adj.P.Val)))
	#print(p + geom_point(alpha = 0.55, color="lightcyan4") +
	#geom_vline(xintercept=log(2), linetype="dashed", color = "red")+
	#geom_vline(xintercept=log(0.5), linetype="dashed", color = "red")+
	#geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray40")+
	#geom_text_repel(
    #data = filter(t, -log10(adj.P.Val) > -log10(0.05) , abs(logFC) >= 2),
    #aes(label=gene_name), size =2,
    #segment.size  = 0.2,
    #segment.color = "grey50")+
	#ggtitle(paste(d$Cancer[1], d$lnc[1], get_name(d$lnc[1]))))

	#only include those with FC >2 or <2
  	t = as.data.table(filter(t, adj.P.Val <=0.05 & ((logFC <= log(0.5)) | (logFC >= log(2)))))

    #if(dim(t)[1] <= 1){
    #	t <- c(paste(d$lnc[1], d$canc[1]), "none")
   	#}

    return(t)
}

#pdf("LGG_HOXA10_volcano_plot_Aug21.pdf")
#diffE(all_canc_lnc_data[[2]])
#dev.off()

#pdf("volcano_plots_diffE_lncRNA_risks.pdf")
#diffEresults = llply(all_canc_lnc_data, diffE, .progress="text")
diffEresults = llply(all_canc_lnc_data, diffE, .progress="text")
#dev.off()

diffEresults1 = ldply(diffEresults, data.frame)
diffEresults1 = as.data.table(diffEresults1)
saveRDS(diffEresults1, file="diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
saveRDS(diffEresults1, file="diff_expressed_PCGs_lncRNA_risk_groups_lgg_nov30.rds")

##########
###DONE###
##########
