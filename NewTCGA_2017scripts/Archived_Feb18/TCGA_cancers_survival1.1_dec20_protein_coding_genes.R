###---------------------------------------------------------------
###TCGA_cancers_survival1.1_dec20_protein_coding_genes.R
###---------------------------------------------------------------

###December 20th, 2017
###establish a list of lncRNAs to evaluate using survival analysis
###by looking at mean versus variance density for each lncRNA 
###within each cancer type 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)
library(scater)

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. lncRNA expression in different cancers 
rna = readRDS("19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

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
sums = apply(rna[,1:(ncol(rna)-4)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]
zeroes = fantom$CAT_geneName[which(fantom$CAT_geneID %in% zeroes)] #should check if these are also the ones that are lowly expressed

#4. Now within each cancer get mean and variance for each gene 
rna = as.data.table(rna)

get_mean_variance = function(cancer){
	canc_data = rna[canc == cancer]
	gene_data = canc_data[,1:(ncol(rna)-4)]
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
	gene_data = log1p(gene_data)
	means = as.data.frame(apply(gene_data, 2, mean)); colnames(means) = c("mean")
	vars = as.data.frame(apply(gene_data, 2, sd)) ; colnames(vars) = c("sd")
	means$gene = rownames(means) ; vars$gene = rownames(vars)
	gene_overview = merge(means, vars, by = "gene")
	gene_overview = merge(gene_overview, fantom, by ="gene")
	gene_overview = as.data.frame(gene_overview)
	gene_overview$canc = cancer
	gene_overview = as.data.table(gene_overview)
	gene_overview = gene_overview[order(mean)]
	return(gene_overview)
}

gene_overviews = llply(cancers, get_mean_variance, .progress="text")

#5. Make density plots for each cancer type 

make_density_plot = function(gene_data){
	#gene_data$mean = log2(gene_data$mean)
	#gene_data$var = log2(gene_data$var)
	sp = ggplot(gene_data, aes(x=mean, y=sd)) + geom_point(alpha = 0.2, aes(colour = factor(functional_evidence)))
	sp = sp + geom_density_2d() + theme_minimal() + 
	labs(title= paste(gene_data$canc[1], "Mean vs SD", nrow(gene_data), "lncRNAs"),
       x="mean log1p(FPKM)", y = "SD log1p(FPKM)")
	return(sp)
}

plots = llply(gene_overviews, make_density_plot, .progress="text")

ml <- marrangeGrob(plots, nrow=2, ncol=2)
ggsave("w_colour_logged1p_allexpression_first_mean_SD_plots_lncRNAS_TCGA.pdf", ml, width = 20, height = 15)
dev.off()

#6. 






























