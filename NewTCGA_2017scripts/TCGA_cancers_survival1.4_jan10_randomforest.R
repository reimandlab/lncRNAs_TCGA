###---------------------------------------------------------------
###TCGA_cancers_survival1.1_dec20.R
###---------------------------------------------------------------

###December 20th, 2017
###establish a list of lncRNAs to evaluate using survival analysis
###by looking at mean versus variance density for each lncRNA 
###within each cancer type 

#first get distribution of expression of PCGs in the first 5% percentile 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. lncRNA expression in different cancers 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")
pcg = readRDS("19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

lncrnas = colnames(rna)[1:(ncol(rna)-5)]
pcgs = colnames(pcg)[1:(ncol(pcg)-5)]

rna$patient = rownames(rna) ; pcg$patient = rownames(pcg)
rna = merge(rna, pcg, by = c("patient", "canc", "time", "status", "sex", "patient"))

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

#3. Remove any lncRNAs that are not expressed in any of the patients 
sums = apply(rna[,6:ncol(rna)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

#4. Now within each cancer get mean and variance for each gene 
rna = as.data.table(rna)

get_mean_variance = function(cancer){
	canc_data = rna[canc == cancer]
	gene_data = canc_data[,6:ncol(rna)]
	#1. remove any genes that have 0 counts within cancer
	sums = apply(gene_data, 2, sum) #134 with 0 expression in ALL patients 
	zeroes = names(sums[which(sums ==0)]) #what are they?
	#in pcawg 
	z <- which(colnames(gene_data) %in% zeroes)
	gene_data = as.data.frame(gene_data)
	if(!(length(z)==0)){
		gene_data = gene_data[,-z]
	}
	#2. calculate mean variance for each gene and save as a new data table 
	#gene_data = log1p(gene_data)
	means = as.data.frame(apply(gene_data, 2, mean)); colnames(means) = c("mean")
	vars = as.data.frame(apply(gene_data, 2, sd)) ; colnames(vars) = c("sd")
	means$gene = rownames(means) ; vars$gene = rownames(vars)
	gene_overview = merge(means, vars, by = "gene")
	#gene_overview = merge(gene_overview, fantom, by ="gene")
	gene_overview = as.data.frame(gene_overview)
	gene_overview$canc = cancer
	gene_overview$type = ""
	z <- which(gene_overview$gene %in% lncrnas) ; gene_overview$type[z] = "lncRNA"
	z <- which(gene_overview$gene %in% pcgs) ; gene_overview$type[z] = "pcgs"
	z <- which(gene_overview$type == "")
	if(!(length(z)==0)){
		gene_overview = gene_overview[-z,]
	}
	gene_overview = as.data.table(gene_overview)
	gene_overview = gene_overview[order(mean)]
	return(gene_overview)
}

cancers = cancers[c(1,7,14)]
gene_overviews = llply(cancers, get_mean_variance, .progress="text")
#top 5 cancers

#5. Calculate how many lncRNA fall within each mean/sd bin 
#10*10 matrix , PCGs and lncRNAs seperatley 

make_bins = function(dt){
	
	lncdt = subset(dt, type == "lncRNA")
	meandt = lncdt[,c(1:2,4:5)]
	sdft = lncdt[,c(1,3,4:5)]
	meandt = meandt[order(mean)]
	sdft = sdft[order(sd)]
	meandt$meanrank = 1:nrow(meandt)
	sdft$sdrank = 1:nrow(sdft)
	meandt$meanper = meandt$meanrank/nrow(meandt)
	sdft$sdper = sdft$sdrank/nrow(sdft)
	lncdt = merge(meandt, sdft, by = c("gene", "canc", "type"))
    lncdt = lncdt[order(meanper)]

	pcgdt = subset(dt, type == "pcgs")
	meandt = pcgdt[,c(1:2,4:5)]
	sdft = pcgdt[,c(1,3,4:5)]
	meandt = meandt[order(mean)]
	sdft = sdft[order(sd)]
	meandt$meanrank = 1:nrow(meandt)
	sdft$sdrank = 1:nrow(sdft)
	meandt$meanper = meandt$meanrank/nrow(meandt)
	sdft$sdper = sdft$sdrank/nrow(sdft)
	pcgdt = merge(meandt, sdft, by = c("gene", "canc", "type"))

	dt = rbind(lncdt, pcgdt)
	return(dt)
}

percentile_ranks = llply(gene_overviews, make_bins, .progress="text")

#6. Change percentile cutoff, obtain list of candidate and do surival 
#analysis 

#What is PCGs low 5percntile?
percentile_compare = function(dtt){
	#what are the genes at the lowest 5% of mean in PCGs
	pcg = subset(dtt, type == "pcgs")
	pcg = subset(pcg, meanper <=0.4)
	#what is the expression of these genes?
	df = subset(rna, rna$canc %in% dtt$canc[1])
	df = as.data.frame(df)
	z = which(colnames(df) %in% pcg$gene)
	df = df[,z]
	#plot summaries of mean and medians 
	means = apply(df, 2, mean) ; medians = apply(df, 2, median) 
	dff = data.frame(gene = names(means), mean = means, median = medians)
	dff$mean = as.numeric(dff$mean) ; dff$median = as.numeric(dff$median)
	mean = ggboxplot(dff, y = "mean") ; median = ggboxplot(dff, y = "median")
	cutoff = summary(dff$mean)[3]

	#what does lncRNAs look like?
	lnc = subset(dtt, type == "lncRNA")
	lnc = subset(lnc, mean >= cutoff)

	#what is the expression of these genes?
	df = subset(rna, rna$canc %in% dtt$canc[1])
	df = as.data.frame(df)
	z = which(colnames(df) %in% lnc$gene)
	df = df[,z]
	#plot summaries of mean and medians 
	means = apply(df, 2, mean) ; medians = apply(df, 2, median) 
	dff = data.frame(gene = names(means), mean = means, median = medians)
	dff$mean = as.numeric(dff$mean) ; dff$median = as.numeric(dff$median)
	mean2 = ggboxplot(dff, y = "mean") ; median2 = ggboxplot(dff, y = "median")
	plot_grid(mean, mean2, median, median2, labels = c("A", "B", "C", "D"), align = "hv")
	dev.off()
	return(as.list(lnc$gene))
}

random_forest = function(dtt){
	lncs = percentile_compare(dtt)
	df <- subset(rna, rna$canc %in% dtt$canc)
  	z <- which(colnames(df) %in% lncs)
  	
  	df = as.data.frame(df)
  	df <- df[,c(1:5,z)]  
  	df[,6:ncol(df)] = log1p(df[,6:ncol(df)])
  	df$status[df$status=="Alive"] <- 0
 	df$status[df$status=="Dead"] <- 1
  	df$status <- as.numeric(df$status)
    df$time <- as.numeric(df$time)
  	
  	library(ranger)
	#ranger model
  	r_fit <- ranger(Surv(time, status) ~ .,     
        data = df,
		mtry = 4,
        importance = "permutation",
        splitrule = "extratrees",
        respect.unordered.factors = "order",
        verbose = TRUE)
  	# Average the survival models
	death_times <- r_fit$unique.death.times 
	surv_prob <- data.frame(r_fit$survival)
	avg_prob <- sapply(surv_prob,mean)

	# Plot the survival models for each patient
	plot(r_fit$unique.death.times,r_fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Days",
     ylab = "survival",
     main = "Patient Survival Curves")

	#
	cols <- colors()
	for (n in sample(c(2:dim(df)[1]), 20)){
  		lines(r_fit$unique.death.times, r_fit$survival[n,], type = "l", col = cols[n])
		}
	lines(death_times, avg_prob, lwd = 2)
	legend(500, 0.7, legend = c('Average = black'))
	

	vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
	names(vi) <- "importance"
	head(vi)






}

#percentile_ranks















