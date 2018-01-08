###---------------------------------------------------------------
###TCGA_cancers_survival1.1_dec20.R
###---------------------------------------------------------------

###December 20th, 2017
###establish a list of lncRNAs to evaluate using survival analysis
###by looking at mean versus variance density for each lncRNA 
###within each cancer type 

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
#10*10 matrix

make_bins = function(dt){
	meandt = dt[,c(1:2,4:5)]
	sdft = dt[,c(1,3,4:5)]
	meandt = meandt[order(mean)]
	sdft = sdft[order(sd)]
	meandt$meanrank = 1:nrow(meandt)
	sdft$sdrank = 1:nrow(sdft)
	meandt$meanper = meandt$meanrank/nrow(meandt)
	sdft$sdper = sdft$sdrank/nrow(sdft)
	dt = merge(meandt, sdft, by = c("gene", "canc", "type"))
}

percentile_ranks = llply(gene_overviews, make_bins, .progress="text")

#6. Make density plots for each cancer type 
make_density_plot = function(gene_data){
	#gene_data$mean = log2(gene_data$mean)
	#gene_data$sd = log2(gene_data$sd)
	sp = ggplot(gene_data, aes(x=meanper, y=sdper)) + geom_point(alpha = 0.3, aes(colour = factor(type)))
	sp = sp + geom_density_2d() + theme_minimal(base_size = 18) + 
	labs(title= paste(gene_data$canc[1], "Mean vs SD", nrow(gene_data), "Total Genes"),
       x="mean percentiles", y = "SD percentiles")
	return(sp)
}

plots = llply(percentile_ranks, make_density_plot, .progress="text")

ml <- marrangeGrob(plots, nrow=1, ncol=1)
ggsave("percentiles_JAN_PCGS_vs_lncRNAs_mean_SD_plots_lncRNAS_TCGA.pdf", ml, width = 20, height = 18)
dev.off()

#7. Density plots for just lncRNAs 
make_bins_lncs = function(dt){
	dt = subset(dt, type == "lncRNA")
	meandt = dt[,c(1:2,4:5)]
	sdft = dt[,c(1,3,4:5)]
	meandt = meandt[order(mean)]
	sdft = sdft[order(sd)]
	meandt$meanrank = 1:nrow(meandt)
	sdft$sdrank = 1:nrow(sdft)
	meandt$meanper = meandt$meanrank/nrow(meandt)
	sdft$sdper = sdft$sdrank/nrow(sdft)
	dt = merge(meandt, sdft, by = c("gene", "canc", "type"))
}

lnc_percentile_ranks = llply(gene_overviews, make_bins_lncs, .progress="text")

plots = llply(lnc_percentile_ranks, make_density_plot, .progress="text")

ml <- marrangeGrob(plots, nrow=1, ncol=1)
ggsave("percentiles_JAN_lncRNAs_mean_SD_plots_lncRNAS_TCGA.pdf", ml, width = 20, height = 18)
dev.off()

#8. Change percentile cutoff, obtain list of candidate and do surival 
#analysis 

surv_test = function(gene){
  
  results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")

  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% dt$canc)
  z <- which(colnames(df) %in% gene)
  if(!(length(z)==0)){
  df = as.data.frame(df)
  df <- df[,c(1:5,z)]  
  df[,6] <- log1p(df[,6])

  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,6]), 0.5)
  if(median2 ==0){
  	median2 = mean(as.numeric(df[,6]))
  }

  #median2 <- median(df[,1])
  for(y in 1:nrow(df)){
    genexp <- df[y,6]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
  gene <- colnames(df)[6]
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

analyze_cands = function(dt){
	dt = subset(dt, meanper >=0.25 & sdper >=0.25)
	genes = as.list(unique(dt$gene))
	surv_genes = llply(genes, surv_test, .progress="text")
	surv_results = ldply(surv_genes, rbind)
}


###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------
































