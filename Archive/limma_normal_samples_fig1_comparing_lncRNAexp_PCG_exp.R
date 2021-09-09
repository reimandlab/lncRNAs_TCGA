source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)
library(plyr)
library(doParallel)

rna = as.data.frame(rna)
saveRDS(rna, file="TCGA_rna_expression_Data_forgtex_analysis.rds")

norm = as.data.frame(norm)
rownames(norm) = norm$patient

#only look at cancers with at least 50 patients in them 
cancers_dist = filter(as.data.table(table(rna$Cancer)), N >=50)
rna = subset(rna, rna$Cancer %in% cancers_dist$V1)
norm = subset(norm, norm$Cancer %in% rna$Cancer)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#add location to lncRNAs 

#how many lncRNAs  --> 5,785

rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))	
z2 = which(colnames(rna) %in% "Cancer")
rna = rna[,c(z1,z2)]

#normal lncRNAs
z = which(colnames(norm) %in% colnames(rna))
norm = norm[,z]
#make sure looking at same lncRNAs 
z = which(colnames(rna) %in% colnames(norm))
rna = rna[,z]
#only keep normal tissues with at least 10 patient samples
pats = as.data.table(table(norm$Cancer))
pats = pats[order(N)]
pats = filter(pats, N >=10)
norm = subset(norm, norm$Cancer %in% pats$V1)

#subset rna file to just cancers in normal file
rna = subset(rna, rna$Cancer %in% norm$Cancer)

#how many detectable lncRNAs in normal tissues versus tumours 

cancers = unique(rna$Cancer)

get_canc_data = function(cancer){
	norm_data = subset(norm, Cancer == cancer)
	norm_data$type = "normal" 
	canc_data = subset(rna, Cancer == cancer)
	canc_data$type = "cancer"
	all_data = rbind(norm_data, canc_data)
	return(all_data)
}

canc_datas = llply(cancers, get_canc_data)

########calculate number of lncRNAs up/downregulated for each cancer type 

#for each cancer types, get number of lncRNAs "differentailly expressed"

###---------------------------------------------------------------------------------------------------------------------------
###----------------------------diff-Expression-Function-----------------------------------------------------------------------
###---------------------------------------------------------------------------------------------------------------------------

get_diff_exp = function(dtt){
	library(stringr)
	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	lncs = as.list(colnames(dtt)[z1]) #5781 
	dtt[,z1] = log1p(dtt[,z1])

	#for each lncRNA calculate mean difference in expression
	#plot expression
	pdf(file=paste(dtt$Cancer[1], "diff_expressed_lncRNAs.pdf", sep="_"))

	#use limma to get differenitally expressed lncRNA between tumour and normals per cancer type
	design <- model.matrix(~ 0 + factor(dtt$type))
	colnames(design) <- c("Cancer", "Normal")

	expression <- t(dtt[,z1])
	fit <- lmFit(expression, design)
	cont.matrix <- makeContrasts(CancervsNormal=Cancer-Normal, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	ps <- fit2$p.value
	ps <- p.adjust(ps, method="fdr")
	numGenes <- length(which(ps <= 0.05))


	genes=rownames(expression)
    t <- topTable(fit2,coef=1,adjust.method="fdr",n=numGenes,p.value=0.05,genelist=genes)

	if(dim(t)[1] > 10){

    #rank list of genes before making heatmap
    t <- as.data.table(t)
    #first by adj p values then by decreasing FC
    #t <- t[order(adj.P.Val)]
    t <- t[order(-abs(as.numeric(logFC)))]
    t = filter(t, abs(logFC) >=2)

    if(dim(t)[1] >= 5){

    	top <- c(paste(dtt$Cancer[1]), t$ID)

    	if(dim(t)[1] >500){
    		tt = t[1:500,]
			}

		if(dim(t)[1] <=500){
			tt = t
		}

    	#generate heatmap 
        heat <- expression[which(rownames(expression) %in% tt$ID),]
    	tags <- dtt$type
		color.map <- function(tags) { if (tags=="normal") "#FF0000" else "#0000FF" }
    	patientcolors <- unlist(lapply(tags, color.map))

		# cluster on correlation

		hc <- hclust(as.dist(1 - cor(t(heat))), method="ward.D2")
		# draw a heatmap
		my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 100)
		heatmap.2(as.matrix(heat), col=my_palette, ColSideColors= patientcolors, cexRow=0.5, cexCol=0.6, Rowv=as.dendrogram(hc), trace="none", scale="row")


    	}#dim(t)[1] >=5
	}#dim(t)[1]>10	

	colnames(t)[1] = "lnc"
	t$Cancer = dtt$Cancer[1]
	dev.off()
	return(t)
}

lncs_summary_all = llply(canc_datas, get_diff_exp, .progress="text")
lncs_summary_all1 = ldply(lncs_summary_all, data.frame)

saveRDS(lncs_summary_all1, "cancer_normals_lncRNA_diff_expression_analysis_June5.rds")

lncs_summary_all1 = readRDS("cancer_normals_lncRNA_diff_expression_analysis_June5.rds")

###---------------------------------------------------------------------------------------------------------------------------
###----------------------------summarize-results------------------------------------------------------------------------------
###---------------------------------------------------------------------------------------------------------------------------

#how many lncRNAs upregulated vs downregulated each cancer

lncs_summary_all1$reg[lncs_summary_all1$logFC >0] = "Upregulated_Cancer"
lncs_summary_all1$reg[lncs_summary_all1$logFC <0] = "Downregulated_Cancer"

sum_genes = as.data.table(table(lncs_summary_all1$Cancer, lncs_summary_all1$reg))
sum_genes = sum_genes[order(-V2, N)]

order = unique(sum_genes$V1)

pdf("sum_diff_expressed_genes_cancers.pdf", width=9)
g = ggbarplot(sum_genes, x="V1", y="N", order=order, color = "grey", fill = "V2", palette = "Paired",
  position = position_dodge(0.9),label.pos = "out", xlab="Cancer", ylab="Number of Genes",  label = TRUE) + theme_light() 

ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65, legend.title="Expression")

dev.off()


###---------------------------------------------------------------------------------------------------------------------------
###----------------------------check-candidates-------------------------------------------------------------------------------
###---------------------------------------------------------------------------------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#which of the candidates are differentially expressed? 

results_cands = as.data.frame(matrix(ncol=ncol(lncs_summary_all1)))
colnames(results_cands) = colnames(lncs_summary_all1)

for(i in 1:nrow(allCands)){
	lnc = as.character(allCands$gene[i])
	canc = allCands$cancer[i]
	z = which((lncs_summary_all1$lnc == lnc) & (lncs_summary_all1$Cancer == canc))
	print(z)
	row = lncs_summary_all1[z,]
	names(row) = colnames(results_cands)
	results_cands = rbind(results_cands, row)
}

results_cands = results_cands[-1,]

#9 lncRNAs from candidates coming from cancer types that have available 
#matched normal tissues are significantly differetnially expressed 
#across 5 cancer types 

colnames(allCands)[c(1, 5)] = c("lnc", "Cancer")

allCands = merge(allCands, results_cands, by=c("lnc", "Cancer"))


#check if biological match 
#if HR > 1 & significantly upregulated in cancer --> potential OG? 
#if HR < 1 & significantly downregulated in cacner --> potential TS? 

allCands$biological_match = ""
z = which((allCands$HR >1) & (allCands$reg == "Upregulated_Cancer"))
allCands$biological_match[z] = "PredictedOG"

z = which((allCands$HR <1) & (allCands$reg == "Downregulated_Cancer"))
allCands$biological_match[z] = "PredictedTS"

#plot boxplots 
#expression for patients wtih hihg expression in cancer
#patients wtih low expression in cancer
#patient with normal tissue expression 

lncs = as.list(as.character(allCands$lnc))

compare_exp_boxplots = function(lnc){

	#get tumour expression and seperate by high and low expression
	z = which(allCands$lnc == lnc)
	canc = allCands$Cancer[z]
	#subset tumour exp
	canc_exp = subset(rna, Cancer==canc)
	canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(lnc, "Cancer"))]
	#sep by median
	median2 = median(as.numeric(canc_exp[,1]))
	canc_exp$tag = ""
	if(median2 == 0){
		z = which(canc_exp[,1] ==0)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}
	if(!(median2==0)){
		z = which(canc_exp[,1] < median2)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}

	#risk 
	risk = allCands$HR[allCands$lnc == lnc]
	if(as.numeric(risk) > 1){
		risk_exp = "High_exp"
	}
	if(as.numeric(risk) < 1){
		risk_exp = "Low_exp"
	}

	canc_exp$risk_type = risk_exp
	canc_exp$exp_type = "Tumour"

	#get normal expression
	#subset tumour exp
	norm_exp = subset(norm, Cancer==canc)
	norm_exp = norm_exp[,which(colnames(norm_exp) %in% c(lnc, "Cancer"))]
	norm_exp$tag = "norm"
	norm_exp$risk_type = "norm"
	norm_exp$exp_type = "normal"

	#all expression data needed for boxplot
	all_exp = rbind(norm_exp, canc_exp)
	all_exp[,1] = log1p(all_exp[,1])

	colnames(all_exp)[1] = "lncRNA"

	#boxplot
	ggboxplot(all_exp, ylab="log1p(FPKM)", x="tag", y="lncRNA", palette = mypal[c(3,2,1)], add = "jitter", fill = "tag", order=c("norm", "Low", "High"), 
		title= paste(lnc, canc, "risk=", risk_exp))+ 
		stat_compare_means(label = "p.signif", 
                     ref.group = "norm")  

}

pdf("lncRNA_expression_tumours_TCGA_matched_normals_cancds.pdf", width=10)
llply(lncs, compare_exp_boxplots)
dev.off()


###---------------------------------------------------------------------------------------------------------------------------
###----------------------------check-pairwise-comparisone-normal-vs-tumour----------------------------------------------------
###---------------------------------------------------------------------------------------------------------------------------

compare_exp_boxplots_pairwise = function(lnc){

	#get tumour expression and seperate by high and low expression
	z = which(allCands$lnc == lnc)
	canc = allCands$Cancer[z]
	#subset tumour exp
	canc_exp = subset(rna, Cancer==canc)
	canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(lnc, "Cancer"))]
	canc_exp$patient = rownames(canc_exp)
	#sep by median
	median2 = median(as.numeric(canc_exp[,1]))
	canc_exp$tag = ""
	if(median2 == 0){
		z = which(canc_exp[,1] ==0)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}
	if(!(median2==0)){
		z = which(canc_exp[,1] < median2)
		canc_exp$tag[z] = "Low"
		canc_exp$tag[-z] = "High"
	}

	#risk 
	risk = allCands$HR[allCands$lnc == lnc]
	if(as.numeric(risk) > 1){
		risk_exp = "High_exp"
	}
	if(as.numeric(risk) < 1){
		risk_exp = "Low_exp"
	}

	canc_exp$risk_type = risk_exp
	canc_exp$exp_type = "Tumour"

	#get normal expression
	#subset tumour exp
	norm_exp = subset(norm, Cancer==canc)
	norm_exp = norm_exp[,which(colnames(norm_exp) %in% c(lnc, "Cancer"))]
	norm_exp$patient = rownames(norm_exp)
	norm_exp$tag = "norm"
	norm_exp$risk_type = "norm"
	norm_exp$exp_type = "normal"

	#keep patients in both datafiles
	z = which(norm_exp$patient %in% canc_exp$patient)
	pats_common = norm_exp$patient[z]

	norm_exp = subset(norm_exp, patient %in% pats_common)
	canc_exp = subset(canc_exp, patient %in% pats_common)

	all_exp = rbind(norm_exp, canc_exp)
	all_exp[,1] = log1p(all_exp[,1])

	colnames(all_exp)[1] = "lncRNA"

	#boxplot

	ggpaired(all_exp, x = "exp_type", y = "lncRNA",
         color = "exp_type", line.color = "gray", line.size = 0.4,
         palette = mypal[c(3,1)], add = "jitter", ylab="log1p(FPKM)", title= paste(lnc, canc, "risk=", risk_exp))+
		 stat_compare_means(paired = TRUE)

}


pdf("lncRNA_pairwise_expression_tumours_TCGA_matched_normals_cancds.pdf", width=10)
llply(lncs, compare_exp_boxplots_pairwise)
dev.off()


#now do the same thing for GTEx


























