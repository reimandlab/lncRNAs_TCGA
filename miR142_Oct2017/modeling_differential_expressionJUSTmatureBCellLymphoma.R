###----------------------------------
###Modeling_differential_expression.R
###----------------------------------

###
###Libraries––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-
###

source("source_file.R")
library(GenomicRanges)
library(stringr)

###
###Data-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––-––
###

#1. ppl with mutations in miR142
load("mirna142_muts_for_Lina_250117.rsav")
#20 mutations, 17 people in total 
muts_ppl = mirna142_muts$mirna_pre 

#2. identified mir142 targets comrpise set of genes tests with annotations 
genes = fread("file_for_network_w_target_type_and_cancer_genes.txt")

#3. UCSC coordinates 
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

#CNA 172 patients
cna = readRDS("172_Lymph_CNA_plus_MutStatus_Histo.rds")

#RNA-Seq 172 patients 
exp = readRDS("172_Lymph_RNASeq.rds")

#AID sampels 
aid = c("DO52717" ,"DO52692" , "DO52682",  "DO221124" ,"DO52675",  "DO52672")
z <- which(colnames(cna) %in% aid)
cna[960,z] = "Nomutation"

#Only keep Mature B-cell lymphoma 
z <- which(cna[961,] == "Mature B-cell lymphoma")
cna = cna[,z]

###
###Analysis-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
###


#LM plotting function 
#Run linear regression and obtain p-values and generate plot
linear_regression_wplot <- function(lm){
	#Make plot 
	print(ggplot(lm$model, aes_string(x = names(lm$model)[2], y = names(lm$model)[1])) + 
  	geom_point() +
  	theme_hc() + 
  	stat_smooth(method = "lm", col = "red") +
  	labs(title = paste(colnames(lm$model)[1], "Adj R2 = ",round(signif(summary(lm)$adj.r.squared, 5),digits=4),
                     "Intercept =",round(signif(lm$coef[[1]],5), digits=4),
                     " Slope =",round(signif(lm$coef[[2]], 5), digits=4),
                     " P =",round(signif(summary(lm)$coef[2,4], 5),digits=4))))
	}

###
###1 - Linear logged models-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
###

##Apply linear model with varying cofactors to list of genes 

gene_list = rownames(exp)

get_data = function(gene){
	#1. make matrix for that gene 
	z <- which(rownames(cna) == gene)
	gexp = t(cna[c(z,960:961),])
	gexp = as.data.frame(gexp)
	gexp$patient = rownames(gexp)
	colnames(gexp) = c(paste(gene, "CNA", sep="_"), "MutationStatus", "Subtype", "patient")
	z2 <- which(rownames(exp) == gene)
	expdat = t(exp[z2,])
	expdat = as.data.frame(expdat)
	expdat$patient = rownames(expdat)
	gexp = merge(gexp, expdat, by = "patient")
	return(gexp)
}

gene_data = llply(gene_list, get_data, .progress = "text")

wilcoxon_testing = function(data){
	res <- wilcox.test(data[,5] ~ MutationStatus, data = data,
                   exact = FALSE)
	p = res$p.value
	results = as.data.frame(matrix(ncol=2, nrow=1))
	colnames(results) = c("gene", "wilcoxon_pval")
	results$gene = colnames(data)[5]
	results$wilcoxon_pval = p

	#for plotting 
	#data$geneE = data[,5]
	#g = ggboxplot(data, x = "MutationStatus", y = "geneE", add="jitter", palette=mypal[c(2,1)], col="MutationStatus", main=paste(colnames(data)[5]))
	#g = g + stat_compare_means()
	#print(g)

	return(results)
}

wilcoxon_results = llply(gene_data, wilcoxon_testing, .progress = "text")
#pdf("Wilcoxon_tests_959genes_mutationVS_Expression.pdf", width=9)
wilcoxon_results = ldply(wilcoxon_results, data.frame)
wilcoxon_results$fdr = p.adjust(wilcoxon_results$wilcoxon_pval, method="fdr")
wilcoxon_results = as.data.table(wilcoxon_results)
wilcoxon_results = wilcoxon_results[order(fdr)]


run_linear_model = function(data){
	#1. log1p 
	data[,5] = log1p(as.numeric(data[,5]))
	#2. print median 
	med = median(data[,5])
	#3. build models 
	m0 = lm(data[,5] ~ 1) 
	m1 = lm(data[,5] ~ data[,3])
	colnames(m1$model)[1] = gene
	#plot
	print(linear_regression_wplot(m1))
	m5 = lm(data[,5] ~ data[,2])
	print(linear_regression_wplot(m5))

	mut_m1_coef = summary(m1)$coefficients[2,1]
	mut_m1_pval = summary(m1)$coefficients[2,4]

	m2 = lm(data[,5] ~ data[,3] + as.numeric(data[,2]))
	colnames(m2$model)[1] = gene

	mut_m2_coef = summary(m2)$coefficients[2,1]
	mut_m2_pval = summary(m2)$coefficients[2,4]

	gene = colnames(data)[5]

	#4. check if gene CNA has NAs 
	z <- length(which(is.na(data[,2])))
	#5. set up anovas 
	just_mut = anova(m0, m1)[2,6]
	if(!(z >0)){
	mut_cna = anova(m1, m2)[2,6]
	}

	if(z >0){
	mut_cna = "not_avail"
	}
	#6. summarize results and print 
	results = as.data.frame(matrix(ncol=6, nrow=2))
	colnames(results) = c("gene", "median", "test", "ANOVApval", "coef", "coef_pval")
	results$gene = colnames(data)[5]
	results$median = med
	results$test = c("just_mut", "mut_cna")
	results$ANOVApval = c(just_mut, mut_cna)
	results$coef = c(mut_m1_coef, mut_m2_coef)
	results$coef_pval = c(mut_m1_pval, mut_m2_pval)

	return(results)

}

pdf("Linear_Model_mutation_CNA.pdf", width=10, height=8)
linear_model_results = llply(gene_data, run_linear_model, .progress = "text")
dev.off()

###
###2 - Negative binomial rounded models-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
###

##Apply negative binomial model with varying cofactors to list of genes 

run_NB_model = function(data){
	#1. floor 
	data[,5] = floor(as.numeric(data[,5]))
	#2. print median 
	med = median(data[,5])
	#make sure that not all expression values are 0 after flooring
	check = dim(as.table(table(data[,5])))
	if(!(check ==1)){
	#3. build models 
	m0 = try(glm.nb(data[,5] ~ 1), TRUE) 
	m1 = try(glm.nb(data[,5] ~ data[,3]), TRUE) 

	if(!(inherits(m1,"try-error"))){
		mut_m1_coef = summary(m1)$coefficients[2,1]
		mut_m1_pval = summary(m1)$coefficients[2,4]
	}

	if(inherits(m1,"try-error")){
		mut_m1_coef = "not_avail"
		mut_m1_pval = "not_avail"
	}

	m2 = try(glm.nb(data[,5] ~ data[,3] + as.numeric(data[,2])), TRUE) 
	if(!(inherits(m2,"try-error"))){
		mut_m2_coef = summary(m2)$coefficients[2,1]
		mut_m2_pval = summary(m2)$coefficients[2,4]
	}

	if(inherits(m2,"try-error")){
		mut_m2_coef = "not_avail"
		mut_m2_pval = "not_avail"
	}

	#4. check if gene CNA has NAs 
	z <- length(which(is.na(data[,2])))
	
	#5. set up anovas 
	
	if((!(inherits(m0,"try-error")))&(!(inherits(m1,"try-error")))){
	just_mut = anova(m0, m1)[2,8]
	}	

	if(!(z >0)){

	if((!(inherits(m1,"try-error")))&(!(inherits(m2,"try-error")))){
	mut_cna = anova(m1, m2)[2,8]
	}

	if(!((!(inherits(m1,"try-error")))&(!(inherits(m2,"try-error"))))){
		mut_cna = "not_avail"
	}
	}

	if(z >0){
	mut_cna = "not_avail"
	}
	}

	if(check ==1){
		mut_m1_coef = "not_avail"
		mut_m2_coef = "not_avail"
		mut_m1_pval = "not_avail"
		mut_m2_pval = "not_avail"

		just_mut = "not_avail"
		mut_cna = "not_avail"
	}	


	#6. summarize results and print 
	results = as.data.frame(matrix(ncol=6, nrow=2))
	colnames(results) = c("gene", "median", "test", "ANOVApval", "coef", "coef_pval")
	results$gene = colnames(data)[5]
	results$median = med
	results$test = c("just_mut", "mut_cna")
	results$ANOVApval = c(just_mut, mut_cna)
	results$coef = c(mut_m1_coef, mut_m2_coef)
	results$coef_pval = c(mut_m1_pval, mut_m2_pval)
	
	return(results)
}

neg_model_results = llply(gene_data, run_NB_model, .progress = "text")

###
###3 - Analyze Results–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––---
###

lm_results = ldply(linear_model_results, data.frame)
lm_results$model = "lm"

neg_results = ldply(neg_model_results, data.frame)
neg_results$model = "NB"

results = rbind(lm_results, neg_results)

#1. Which genes significant when just mutation is predictor 
results = as.data.table(results)

#Correct for multiple testing using FDR
#within each test/model combintation 

#lm - mut only
lm_mut_only = filter(results, test == "just_mut", model=="lm")
lm_mut_only$anova_fdr = p.adjust(lm_mut_only$ANOVApval, method="fdr")
lm_mut_only$coef_fdr = p.adjust(lm_mut_only$coef_pval, method="fdr")

#lm - mut + cna
lm_mut_cna_only = filter(results, test == "mut_cna", model=="lm")
lm_mut_cna_only$anova_fdr = p.adjust(lm_mut_cna_only$ANOVApval, method="fdr")
lm_mut_cna_only$coef_fdr = p.adjust(lm_mut_cna_only$coef_pval, method="fdr")

#------NB
#NB - mut only
NB_mut_only = filter(results, test == "just_mut", model=="NB")
NB_mut_only$anova_fdr = p.adjust(NB_mut_only$ANOVApval, method="fdr")
NB_mut_only$coef_fdr = p.adjust(NB_mut_only$coef_pval, method="fdr")

#NB - mut + cna
NB_mut_cna_only = filter(results, test == "mut_cna", model=="NB")
NB_mut_cna_only$anova_fdr = p.adjust(NB_mut_cna_only$ANOVApval, method="fdr")
NB_mut_cna_only$coef_fdr = p.adjust(NB_mut_cna_only$coef_pval, method="fdr")

results = rbind(lm_mut_only, lm_mut_cna_only,
	NB_mut_only, NB_mut_cna_only)

colnames(genes)[2] = "gene"
results = merge(results, genes, by = "gene")

results = as.data.table(results)
results = results[order(coef_fdr)]

results$ANOVApval = round(as.numeric(results$ANOVApval), digits=4)
results$coef = round(as.numeric(results$coef), digits=4)
results$coef_pval = round(as.numeric(results$coef_pval), digits=4)
results$anova_fdr = round(as.numeric(results$anova_fdr), digits=4)
results$coef_fdr = round(as.numeric(results$coef_fdr), digits=4)

#first which genes had sig coefficient before adding new cofactors
mut_only = as.data.table(filter(results, test == "just_mut", coef_fdr <=0.05))
mut_cna = as.data.table(filter(results, test == "mut_cna", coef_fdr <=0.05))














