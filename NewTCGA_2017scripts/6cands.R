#list of 6 candidates 

genes = c("AC006126.4", "ADORA2A-AS1", "NEAT1", "LINC00665", "ZNF503-AS2", "GS1-251I9.4")

###
###PCG Processing---------------------------------------------------------------------------
###

#1. shorten ids 
#Change patient ids to shorted id
change = function(rowname){
  new = canc_conversion$id[which(canc_conversion$TCGA_id %in% rowname)]
  return(new)  
}

rownames(pcg) = sapply(rownames(pcg), change)
#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) #600 to remove 
z <- which(rownames(pcg) %in% ids_remove) #666 PCAWG samples in this TCGA RNA file
pcg = pcg[-z,]

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(pcg) %in% clin$bcr_patient_barcode)
pcg = pcg[z,] #all have clinical data - 7387 patients 

#Add survival info to rna file
pcg = as.data.frame(pcg)
pcg$canc = "" 
pcg$time = ""
pcg$status = ""
pcg$sex = ""

for(i in 1:dim(pcg)[1]){
	pat = rownames(pcg)[i]
	z <- which(clin$bcr_patient_barcode == pat)
	if(length(z) >1){z = z[1]}
	status = clin$vital_status[z]
	time = clin$days_to_death[z]
	if(is.na(time)){
		time = clin$days_to_last_followup[z]
	}
	pcg$status[i] = status
	pcg$time[i] = time
	pcg$sex[i] = clin$gender[z]
}

for(i in 1:dim(pcg)[1]){
	pat = rownames(pcg)[i]
	z = which(canc_conversion$id == pat)
	canc = canc_conversion$Cancer[z]
	pcg$canc[i] = canc
}

#2. remove 0s
sums = apply(pcg[,1:54564], 2, sum)
sums = as.numeric(sums)
z <- which(sums == 0) #134
#remove - MAYBE INSTEAD OF MEDIAN = 0 , SHOULD REMOVE THE ONES 
#THAT HAVE SUM OF 0 MEANING IT HAS 0 EXPRESSION IN EVERY SINGLE PATIENT 
pcg = pcg[,-z]

#still have 52,766 genes 
#within each cancer get the ones with highest variance 
z <- which(pcg$canc == "")
pcg = pcg[-z,]

cancers = unique(pcg$canc)

#only need liver, kidney and ovary 
cancers = cancers[c(6, 9, 14)]

#get list of PCGs
get_pgs = function(cancer){
	dat = subset(pcg, pcg$canc %in% cancer)	
	#calculate variance of each gene 
	vars = apply(dat[,1:52762], 2, var)
	#keep top 1000 genes with highest variance 
	vars = as.data.frame(vars)
	vars$gene = rownames(vars)
	vars = as.data.table(vars)
	vars = vars[order(-vars)]
	vars = vars[1:1000,]
	colnames(vars)[1] = cancer
	return(vars)
}

pcg_lists = llply(cancers, get_pgs)

###
###Analysis-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
###

###
###1 - Linear logged models-––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
###

##Apply linear model with varying cofactors to list of genes 

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

#plot wilcoxon results so that they are ordered
wilcoxon_results = filter(wilcoxon_results, fdr <=0.1)

pdf("172_patients_pvales0.1wilcoxon_boxplots.pdf", width=9)
for(i in 1:nrow(wilcoxon_results)){
	gene = wilcoxon_results$gene[i]
	z <- which(gene_list %in% gene)
	data = gene_data[[z]]
	res <- wilcox.test(data[,5] ~ MutationStatus, data = data,
                   exact = FALSE)
	p = res$p.value
	#for plotting 
	data$geneE = data[,5]
	g = ggboxplot(data, x = "MutationStatus", y = "geneE", add="jitter", palette=mypal[c(2,1)], col="MutationStatus", main=paste(colnames(data)[5], "FDR=", round(wilcoxon_results[i,3], digits=4)))
	g = g + stat_compare_means()
	print(g)
}
dev.off()


run_linear_model = function(data){
	#1. log1p 
	data[,5] = log1p(as.numeric(data[,5]))
	#2. print median 
	med = median(data[,5])
	#3. build models 
	m0 = lm(data[,5] ~ 1) 
	m1 = lm(data[,5] ~ data[,3])
	colnames(m1$model)[1] = gene

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

linear_model_results = llply(gene_data, run_linear_model, .progress = "text")



#3. Co-Expression analysis 

