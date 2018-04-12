###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

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

#####################################################################
###PCA starting with all genes
#####################################################################

groups <- as.factor(rna$Cancer)

#get lncRNAs that have median exp >0 in at least one cancer 
cancers = unique(rna$Cancer)
get_meds_high = function(canc){
	exp = subset(rna, Cancer == canc)
	meds = apply(exp[,2:(ncol(exp)-34)], 2, median)
	meds = as.data.frame(meds)
	meds$gene = rownames(meds)
	z = which(meds$meds >=5)
	meds$med = ""
	meds$med[z] = "Yes" #if median greater than 1 
	meds$med[-z] = "No" #if median less than 1 
	meds$canc = canc
	return(meds)
}

meds_summaries = llply(cancers, get_meds_high)
meds_summaries <- ldply(meds_summaries, data.frame)

#keep only those detectable at least in one cancer type 
meds_summaries = as.data.table(meds_summaries)
meds_summaries = filter(meds_summaries, med == "Yes")

#summary of how many are expressed in each cancer 
summ = as.data.table(table(meds_summaries$gene))
summ = summ[order(N)]

z = which(colnames(rna) %in% meds_summaries$gene)
rna = rna[,c(z, 1, 5787:5820)]

####log1p####
rna[,1:(ncol(rna)-35)] = log1p(rna[,1:(ncol(rna)-35)])

####look at each class of lncRNAs seperatyley? 
z = which(fantom$gene %in% colnames(rna))
fantom = fantom[z,]
#both category and class
lncs = as.data.table(table(fantom$CAT_geneCategory, fantom$CAT_geneClass))
lncs = lncs[order(N)]
lncs = filter(lncs, N >0)
colnames(lncs) = c("CATgeneCategory", "CATgeneClass", "Frequency")
cats = unique(lncs$CATgeneCategory)

get_pca_plot = function(cat){
	lnc_genes = fantom[which(fantom$CAT_geneCategory == cat),1]
	z = which(colnames(rna) %in% lnc_genes)
	exp = rna[,c(z, 5188:5222)]
	groups = as.factor(exp$Cancer)
	res.pca <- prcomp(exp[,1:(ncol(exp)-35)],  scale = TRUE)
	p = fviz_pca_ind(res.pca, geom="point", label="none", habillage=exp$Cancer, addEllipses=TRUE, ellipse.level=0.95) +
  	labs(title = paste("TCGA", cat, "lncRNAs PCA")) #+ 
   	#xlim(-25, 15) + ylim (-10, 20)
	# Change group colors using RColorBrewer color palettes
	#print(p + scale_color_brewer(palette="Dark2") +
    # theme_minimal())
	print(p + theme_minimal())
	print(cat)
}

pdf("ALL_cancers_PCA_using_all_lncRNAs.pdf", width = 15, height=15)
llply(cats, get_pca_plot, .progress="text")
dev.off()

#####################################################################
###Expression distribution for each cancer type lncRNAs versus PCGs
#####################################################################

pdf("all_TCGA_lncs_categories_classes.pdf", width=8, height=8)
ggbarplot(lncs, "CATgeneCategory", "Frequency",
  fill = "CATgeneClass", color = "CATgeneClass", palette = "Dark2",
  label = TRUE,
  position = position_dodge(0.9))
dev.off()

####################################################
#2. get median FPKM for each gene 
#####PLOTTING#######################################

pdf("dist_lncRNA_medians_TCGA_ind_cancers.pdf", width=5, height=5)
for(i in 1:length(unique(rna$Cancer))){
canc = subset(rna, Cancer == unique(rna$Cancer)[i])
meds_lncs = apply(canc[,1:5187], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$gene = rownames(meds_lncs)
meds_lncs$type = ""
for(i in 1:nrow(meds_lncs)){
	z = which(fantom$gene == meds_lncs$gene[i])
	meds_lncs$type[i] = fantom$CAT_geneCategory[z]
}
colnames(meds_lncs)[1] = "median" #<---- logged already, also lncs with med <5 in individual cancers were removed 
meds_lncs = as.data.table(meds_lncs)

pcg_canc = pcg[which(pcg$patient %in% canc$patient),]
pcg_canc[,2:19348] = log1p(pcg_canc[,2:19348])
meds_pcgs = apply(pcg_canc[,2:19348], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "median"
meds_pcgs = as.data.table(meds_pcgs)

all_meds = rbind(meds_lncs, meds_pcgs)
#all_meds$median = floor(all_meds$median)
all_meds$type = as.factor(all_meds$type)

#remove super highly expressed outliers
z = which(all_meds$median >= 1000)
if(!(length(z)==0)){
all_meds = all_meds[-z,]
}

#zoomed in 
all_meds$gene = NULL
gg <- ggplot(all_meds) + labs(title=canc$Cancer[1], y = "density", x="median log1p(FPKM)")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(ggpar(gg, xlim=c(0,100)))

#zoomed out
gg <- ggplot(all_meds) + labs(title=canc$Cancer[1], y = "density", x="median log1p(FPKM)")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
#ggpar(p, xlim=c(0,250))
print(gg)

#3. log medians first 
#all_meds$median = log1p(all_meds$median)
#zoomed in 
#gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="log1p(median FPKM)")
#gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
#gg <- gg + theme_bw()
#print(gg)
#print(ggpar(gg, xlim=c(0,100)))
}
dev.off()




####
#pca with just lncRNAs that have median >5 in only one cancer type? 








































###PCA using lncRNAs that are univariatley associated with survival 



















