###---------------------------------------------------------------
###008_clustering_lncRNA_expression_Figure1a.R 
###---------------------------------------------------------------

#what? 
#Using lncRNA expression profiles run t-SNE to get unsupervised clustering 

#load all required data
source("universal_LASSO_survival_script.R")

#load libraries 
require(caTools)
library(survAUC)
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

date = Sys.Date()

library("FactoMineR")

mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = unlist(cancers[which(!(cancers %in% canc_rm))])

#rna = pcg
rna = subset(rna, Cancer %in% cancers)

dim(rna)
rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% "type")
rna = as.data.frame(rna)
rna = rna[,c(z1, z2)]

#1. remove those not expressed at all and get top 1000 variable lncRNAs 
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, var) #get 500 most variables genes 
sums = sums[order(-sums)]
#keep = sums[1:1000] 
#keep = names(keep)
#z = which(colnames(rna) %in% c("type", keep))
#rna = rna[,z]

#2. cluster cancer types, ie -> each dot should be cancer type 

z1 = which(str_detect(colnames(rna), "ENSG"))
logged_rna = rna
logged_rna[,z1] = log1p(logged_rna[,z1])

z1 = which(str_detect(colnames(logged_rna), "ENSG"))

#try t-SNE 
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(logged_rna[,z1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500) # Run TSNE
saveRDS(tsne_out, file="tsne_out_28_cancers_Oct1_lncRNAs.rds")

tsne_plot = readRDS("tsne_out_28_cancers_Oct1_lncRNAs.rds")
# Show the objects in the 2D tsne representation
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = logged_rna$type)

z = which(duplicated(tsne_plot$col))
tsne_plot$label[z] ="no"
tsne_plot$label[-z] = "yes"

z = (which(tsne_plot$label == "yes"))
for(i in 1:length(z)){
	canc = tsne_plot$col[z[i]]
	tsne_plot$label[z[i]] = canc
}

#assign color to variable ******
pdf("tSNE_28_cancers_april10_top_var_1000_lncRNAs.pdf", width=9)
ggplot(tsne_plot,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col)) + scale_colour_manual(values=mypal)+
geom_text_repel(data = subset(tsne_plot, !(label == "no")))
dev.off()
