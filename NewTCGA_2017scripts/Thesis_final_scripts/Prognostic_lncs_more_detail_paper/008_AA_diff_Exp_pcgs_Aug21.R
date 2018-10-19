#008_AA_diff_Exp_pcgs_Aug21.R

#------------------------------------------------------------------------------
#Take each lncRNA candidate and conduct LIMMA diff exp between lncRNA
#high and low risk groups
#------------------------------------------------------------------------------

library(edgeR)
library(limma)
#library(Glimma)
#library(gplots)
#library(org.Mm.eg.db)
#library(RColorBrewer)

#source code
source("check_lnc_exp_cancers.R")
library(corrplot)

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

#HiC data 
load("hic_data.rsav")
#remove rownames
rownames(hic_data) = c(1:nrow(hic_data))

pcg_counts = readRDS("counts_19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
z = which(str_detect(colnames(all), "ENSG"))	
all[,z] <- log1p(all[,z])

#2. Get lncRNA - median within each tissue type
allCands$combo = paste(allCands$gene, allCands$cancer, sep="_")
combos = unique(allCands$combo)

#canc type to cancer conversion
canc_conv = unique(rna[,c("type", "Cancer")])

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group
#---------------------------------------------------------
#Function 1
#for each lnc-cancer, label patient as lncRNA-risk or non-risk 
#---------------------------------------------------------

get_lnc_canc = function(comb){
	lnc = unlist(strsplit(comb, "_"))[1]
	canc = unlist(strsplit(comb, "_"))[2]
	canc_type = canc_conv$type[which(canc_conv$Cancer == canc)]
	canc_data = subset(pcg_counts, type == canc_type)
	
	z = which(colnames(rna) %in% c("patient", lnc, "type"))
	lnc_dat = rna[,z]
	z = which(lnc_dat$type == canc_type)
	lnc_dat = lnc_dat[z,]
	med = median(as.numeric(lnc_dat[,2]))
	if(med ==0){
		z = which(lnc_dat[,2] > 0)
		lnc_dat$lnc_tag = ""
		lnc_dat$lnc_tag[z] = "high"
		lnc_dat$lnc_tag[-z] = "low"
	}
	
	if(!(med ==0)){
		z = which(lnc_dat[,2] >= med)
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

all_canc_lnc_data = mclapply(combos, get_lnc_canc, mc.cores = 3L)

#---------------------------------------------------------
#Function 2
#wtihin each cancer 
#calculate for each lncRNAs differentially expressed PCGs
#---------------------------------------------------------

diffE <- function(d){
	
	z = which(str_detect(colnames(d), "ENSG"))	
	rownames(d) = d$patient
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

	design <- model.matrix(~ 0 + factor(d$lnc_tag))
	colnames(design) <- c("high", "low")
	rownames(d) <- d$patient

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
diffEresults = mclapply(all_canc_lnc_data, diffE, mc.cores = 3L)
#dev.off()

diffEresults1 = ldply(diffEresults, data.frame)
diffEresults1 = as.data.table(diffEresults1)
saveRDS(diffEresults1, file="diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")

##########
###DONE###
##########

#--------------------------------------------------------------------------------------------------
#-----------------------------------old analysis below---------------------------------------------
#using linear regression with binary predictor 
#--------------------------------------------------------------------------------------------------

get_pcg_enrichment = function(dat){
	lncs = unique(dat$lnc)
	
	get_pcgs_high = function(lncrna){
		newdat = subset(dat, lnc == lncrna)
		#which PCGs have higher expression in the high risk group
		z = which(str_detect(colnames(newdat), "ENSG"))
		meds = apply(newdat[,z], 2, median)
		z = which(meds <= log1p(5))
		rm = names(meds)[z]
		z = which(colnames(newdat) %in% rm)
		newdat = newdat[,-z]

		pcgs = colnames(newdat)[which(str_detect(colnames(newdat), "ENSG"))]
		pcgs_exps = newdat[,c(pcgs)]
		#medians = apply(pcgs_exps, 2, median)
		#pcgs = names(medians[medians > 2])
		print(length(pcgs))
		#lnc_pcg_results = as.data.frame(matrix(ncol=5)) ; colnames(lnc_pcg_results) = c("lnc", "pcg", "canc", "mean_diff", "pvalue")
		#pcgs = pcgs[1:10]
		get_correlation = function(pcg){
			p = pcg
			z = which(colnames(newdat) %in% p)
			lncpcg = newdat[,c(z, 1, (ncol(newdat)-4):ncol(newdat))]
			colnames(lncpcg)[1] = "pcgExp"
			order = c("RISK", "noRISK")
			lncpcg$risk <- factor(lncpcg$risk, levels = order)
			
			fit <- lm(pcgExp ~ risk, data=lncpcg)
			#get p-value and generate boxplot with wilcoxon p-value 
			fit_pval = summary(fit)$coefficients[2,4]
			#which group is it higher in? 
			mean_diff = mean(lncpcg$pcgExp[lncpcg$risk == "RISK"]) - mean(lncpcg$pcgExp[lncpcg$risk == "noRISK"])
			#if higher than 0 --> more expressed in risk group, less than 0 --> more expressed in low risk group
			#g = ggboxplot(lncpcg, x = "risk", y="pcgExp", color="median", title=paste(lncpcg$lnc[1], p, lncpcg$canc[1]))
			#g = g + stat_compare_means()
			#print(g)
			row = c(lncpcg$lnc[1], p, lncpcg$canc[1], mean_diff, fit_pval)
			return(row)
			#names(row) = colnames(lnc_pcg_results)
			#lnc_pcg_results = rbind(lnc_pcg_results, row)
		
		}#end get_correlation function 
		
		pcgs_results = llply(pcgs, get_correlation, .progress="text")
		pcgs_results1 = as.data.frame(do.call("rbind", pcgs_results))
		colnames(pcgs_results1) = c("lnc", "pcg", "canc", "mean_diff", "pvalue")
		return(pcgs_results1)

	} #end get_pcgs_high function
	
	results_lncs = llply(lncs, get_pcgs_high, .progress="text")
	results_lncs1 = as.data.frame(do.call("rbind", results_lncs))
	#results for all lncRNA-PCG correlations in a single cancer type 
	return(results_lncs1)

}

#all_canc_lnc_data = all_canc_lnc_data[1:2] ###TEST CASE -------------------------------------------------------------
#all_canc_lnc_data = llply(all_canc_lnc_data, get_pcg_enrichment, .progress="text")

#all_canc_lnc_data1 = as.data.frame(do.call("rbind", all_canc_lnc_data))
#z = which(all_canc_lnc_data1$lnc %in% cands_dups)

#if(!(length(z))==0){
#all_canc_lnc_data1$lnc[z] = paste(all_canc_lnc_data1$lnc[z], all_canc_lnc_data1$canc[z], sep="_")
#}

#saveRDS(all_canc_lnc_data1, file="all_results_for_each_cancer_from_coexpression_analysis_july19_allCands.rds")

#For each cancer type, for each lncRNA ... 
#Summarize #of PCGs upregulated in risk group and #of PCGs upregulated in non-risk group 

#divide into high risk and low risk lncRNAs
#high_risk = subset(all_canc_lnc_data1, mean_diff >=1.5) #should set higher mean difference threshold?
#low_risk = subset(all_canc_lnc_data1, mean_diff <=0.75) #should set higher mean difference threshold? 

#library(reshape2)

#pcgs enriched in high risk lncRNAs 
#high_risk_matrix = acast(high_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
#      value.var = 'pvalue', fill = 'na')

#pcgs enriched in low risk lncRNAS
#low_risk_matrix = acast(low_risk, pcg ~ lnc, function(x) {sort(as.character(x))[1]},
#      value.var = 'pvalue', fill = 'na')

#columns are lncRNAs and rows are PCGs

#saveRDS(high_risk_matrix, file="high_risk_matrix_lncRNA_candidates_June6.rds")
#saveRDS(low_risk_matrix, file="low_risk_matrix_lncRNA_candidates_June6.rds")







