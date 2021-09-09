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

rna = as.data.frame(rna)

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

#need to get median expression for each gene for each cancer type 
cancers = unique(rna$Cancer)

get_canc_data = function(cancer){
	canc_data = subset(rna, Cancer == cancer)
	return(canc_data)
}

canc_datas = llply(cancers, get_canc_data)

#get median value for each lncRNA
colnames(ucsc)[6] = "gene"
ucsc = ucsc[,c(6, 2, 4,5)]
colnames(ucsc) = c("gene", "chr", "start", "end")

#need to figure out % of people that have expression greater than 1 FPKM 
#in the cohort 

get_medians = function(dtt){
	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	#dtt[,z1] = log1p(dtt[,z1])
	#meds 
	lncs = colnames(dtt)[1:5785]
	calc_freq = function(lnc){
		print(lnc)
		newdat = dtt[,which(colnames(dtt) %in% lnc)]
		l = (length(which(newdat >=1)))/nrow(dtt)
		if(l > 0.1){
			stat = "detectable"
		}
		if(l <= 0.1){
			stat = "NOTdetectable"
		}
		return(stat)
	}

	summary_dat = as.data.frame(matrix(ncol = 3, nrow=length(lncs))) ; colnames(summary_dat) = c("cancer", "lncRNA", "status")
	summary_dat$cancer = dtt$Cancer[1]
	summary_dat$lncRNA = lncs
	summary_dat$status = unlist(llply(lncs, calc_freq))

	#meds = apply(dtt[,z1], 2, median)
	#meds = as.data.frame(meds)
	#meds$canc = dtt$Cancer[1]
	#meds$gene = rownames(meds)
	#meds = merge(meds, ucsc, by="gene")
	return(summary_dat)

}

meds_cancers = llply(canc_datas, get_medians, .progress="text")
meds_cancers1 = ldply(meds_cancers, data.frame)
meds_cancers1$chr =  str_replace(meds_cancers1$chr, "chr", "")
meds_cancers1$chr[meds_cancers1$chr == "X"] = 23
meds_cancers1$chr[meds_cancers1$chr == "Y"] = 24
meds_cancers1$chr = as.numeric(meds_cancers1$chr)
meds_cancers1 = as.data.table(meds_cancers1)
meds_cancers1 = meds_cancers1[order(chr, start, end)]
order = unique(meds_cancers1$gene)
meds_cancers1$gene <- factor(meds_cancers1$gene, levels = order)
meds_cancers1$canc = as.factor(meds_cancers1$canc)
#meds_cancers1 = filter(meds_cancers1, meds > 7)

#convert to heatmap matrix --> rows = cancers and columns = genes
heatmap_matrix = as.data.frame(matrix(ncol = length(unique(meds_cancers1$gene)), nrow=length(unique(meds_cancers1$canc))))
colnames(heatmap_matrix) = unique(meds_cancers1$gene)
rownames(heatmap_matrix) = unique(meds_cancers1$canc)

#learn how to do this using reshape
for(i in 1:nrow(meds_cancers1)){
	x = which(rownames(heatmap_matrix) == meds_cancers1$canc[i])
	y = which(colnames(heatmap_matrix) == meds_cancers1$gene[i])
	heatmap_matrix[x,y] = meds_cancers1$meds[i]
}

# Example: grouping from the first letter:
my_group = as.numeric(as.factor(substr(rownames(data), 1 , 1)))
my_col=brewer.pal(9, "Set1")[my_group]

heatmap(heatmap_matrix, Colv = NA, Rowv = NA, scale="column")

#remove those with 0s in all cancers
sums = apply(heatmap_matrix, 2, sum)
zeroes = sums[which(sums ==0)]
z = which(colnames(heatmap_matrix) %in% names(zeroes))
heatmap_matrix = heatmap_matrix[,-z]
heatmap_matrix = as.matrix(heatmap_matrix)

heatmap_matrix = t(heatmap_matrix)
hc <- hclust(as.dist(1 - cor(t(heatmap_matrix))), method="ward.D2")
# draw a heatmap
pdf("clustering_32_TCGA_cancers_lncRNA_medians_May4.pdf")
my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 100)
heatmap.2(as.matrix(heatmap_matrix), col=my_palette, cexRow=0.01, cexCol=0.5, Rowv=as.dendrogram(hc), trace="none", scale="row", srtCol=30)
dev.off()


