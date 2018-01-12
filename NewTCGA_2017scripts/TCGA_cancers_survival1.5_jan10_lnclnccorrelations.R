###---------------------------------------------------------------
###TCGA_cancers_survival1.5_jan10_lnclnccorrelations.R
###---------------------------------------------------------------

###December 20th, 2017
###establish a list of lncRNAs to evaluate using survival analysis
###by looking at mean versus variance density for each lncRNA 
###within each cancer type 

#identify correlations between lncRNAs to get list of those
#that represent other ones since they are correlated 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. lncRNA expression in different cancers 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

lncrnas = colnames(rna)[1:(ncol(rna)-5)]

rna$patient = rownames(rna) ; 

#matched normal 
norm = readRDS("5919_lncs4matched_normal_tissues_TCGAnew.rds")

#3. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

#4. Fantom data 
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
colnames(fantom)[1] = "gene"

###---------------------------------------------------------------
###Get lncRNAs for each cancer 
###---------------------------------------------------------------

#1. Remove patients with unlabelled cancer type 
z <- which(rna$canc == "")
rna = rna[-z,]
#2. Remove patients without survival status 
z <- which(rna$status == "")
rna = rna[-z,]

cancers = as.list(unique(rna$canc)) #18 cancers wtih at least 90 patients in each cohort
cancers = cancers[c(6,14,9)]

#3. Remove any lncRNAs that are not expressed in any of the patients 
sums = apply(rna[,1:(ncol(rna)-5)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

#4. Now within each cancer get mean and variance for each gene 
rna = as.data.table(rna)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

get_correlated_lncs = function(cancer){
	canc_data = rna[canc == cancer]
	gene_data = canc_data[,1:(ncol(rna)-5)]
	#1. remove any genes that have 0 counts within cancer
	sums = apply(gene_data, 2, sum) #134 with 0 expression in ALL patients 
	zeroes = names(sums[which(sums ==0)]) #what are they?
	#in pcawg 
	z <- which(colnames(gene_data) %in% zeroes)
	gene_data = as.data.frame(gene_data)
	if(!(length(z)==0)){
		gene_data = gene_data[,-z]
	}
	#2. calculate correlation between each lncRNA  
	#first remove the very lowly expressed genes 
	medians = apply(gene_data, 2, median)
	means = apply(gene_data, 2, mean)  
	#remove the ones with median less than the 1st quantile 
	first = summary(means)[2]
	firsts = names(means[which(means <= first)]) #what are they?
	#in pcawg 
	z <- which(colnames(gene_data) %in% firsts)
	gene_data = as.data.frame(gene_data)
	if(!(length(z)==0)){
		gene_data = gene_data[,-z]
	}
	#3. log1p values to make them closer to each other 
	gene_data = log1p(gene_data)
	#4. calculate pairwise correlation between each lncRNA 
	res <- rcorr(as.matrix(gene_data), type = "spearman")
	# Extract the correlation coefficients
	res2= flattenCorrMatrix(res$r, res$P)
	res2$fdr = p.adjust(res2$p, method="fdr")
	res2 = subset(res2, fdr <=0.01)
	res2 = as.data.table(res2)
	res2 = res2[order(fdr)]
	#5. only keep correlations greater than 0.4
	res2 = as.data.frame(res2)
	res2 = subset(res2, abs(res2$cor) >=0.4)
	#6. which lncRNAs are connected to the most nodes? 
	rows = as.data.frame(table(res2$row))
	cols = as.data.frame(table(res2$column))
	all = merge(rows, cols, by = "Var1")
	all$total = all$Freq.x + all$Freq.y
	all = as.data.table(all)
	all = all[order(-total)]
	mediant = summary(all$total)[3]
	all = subset(all, total >= mediant) #937 lncRNAs with the most signficant, high, correlations
	#with other lncRNAs 
	all$canc = cancer
	return(all)
}

correlated_lncs = llply(cancers, get_correlated_lncs, .progress="text")

surv_test = function(gene){
  print(gene)
  results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% corlncs$canc)
  z <- which(colnames(df) %in% gene)
  if(!(length(z)==0)){
  df = as.data.frame(df)
  df <- df[,c(z,5786:5790)]  
  df[,1] <- log1p(df[,1])

  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
  if(median2 ==0){
  	median2 = mean(as.numeric(df[,1]))
  }

  #median2 <- median(df[,1])
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
  gene <- colnames(df)[1]
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
      
  #cox regression 
  res.cox <- coxph(Surv(time, status) ~ median, data = df)
  row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)], df$canc[1])
  names(row) <- names(results_cox)
}
 return(row)
}

calculate_survival = function(corlncs){
	
	canc_data = rna[canc == cancer]
	gene_data = canc_data[,1:(ncol(rna)-5)]
	
	#1. remove any genes that have 0 counts within cancer
	sums = apply(gene_data, 2, sum) #134 with 0 expression in ALL patients 
	zeroes = names(sums[which(sums ==0)]) #what are they?
	#in pcawg 
	z <- which(colnames(gene_data) %in% zeroes)
	gene_data = as.data.frame(gene_data)
	if(!(length(z)==0)){
		gene_data = gene_data[,-z]
	}

	genes = as.list(colnames(gene_data))
	tests_survival2 = llply(genes, surv_test, .progress="text")
	tests_survival3 = ldply(tests_survival2, rbind)
	tests_survival3$pval = as.numeric(tests_survival3$pval)

	tests_survival3$sig = ""
	tests_survival3$sig[tests_survival3$pval <=0.05] = "yes"
	tests_survival3$sig[tests_survival3$pval > 0.05] = "no"

	tests_survival3$fdr = as.numeric(p.adjust(tests_survival3$pval, method="fdr"))
	tests_survival3$fdrsig = ""
	tests_survival3$fdrsig[tests_survival3$fdr <=0.1] = "yes"
	tests_survival3$fdrsig[tests_survival3$fdr > 0.1] = "no"

}



