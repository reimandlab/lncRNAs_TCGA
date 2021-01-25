source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#------FEATURES-----------------------------------------------------

print(colours_palette)

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

#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, var) #get 500 most variables genes
sums = sums[order(-sums)]
#keep = sums[1:1000] #1000 most variable lncRNAs
keep = sums #1000 most variable lncRNAs
keep = names(keep)
z = which(colnames(rna) %in% c("type", keep))
rna = rna[,z]
z1 = which(str_detect(colnames(rna), "ENSG"))
rna[,z1] = log1p(rna[,z1])

library(umap)

set.seed(100)

cancer.labels = rna$type
df = rna
df$type = NULL
embedding = umap(df)
head(embedding)

layout = as.data.table(embedding$layout) ; colnames(layout)=c("x", "y")
layout$col = cancer.labels

z = which(duplicated(layout$col))
layout$label[z] ="no"
layout$label[-z] = "yes"

z = (which(layout$label == "yes"))
for(i in 1:length(z)){
	canc = layout$col[z[i]]
	layout$label[z[i]] = canc
}

layout$cancer=layout$col
layout$cancer=factor(layout$cancer)

pdf("/u/kisaev/Jan2021/UMAP_29_cancer_types_logged_fpkmuq_all_lncs.pdf", width=8)
ggplot(layout,aes(x, y, label = label)) +
geom_point(aes(x=x, y=y, colour=cancer),alpha=0.4, stroke=0) +
#scale_colour_manual(values=mypal)+
colScale_full+theme_bw()+
geom_text_repel(data = subset(layout, !(label == "no")))
dev.off()


write.table(layout, file="/u/kisaev/Jan2021/umap_29_cancers_data_for_figure.txt", quote=F, row.names=F, col.names=T, sep="\t")
















#2. remove those with MAD < 0?
#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
#sums = apply(rna[,z1], 2, mad)
#z = which(sums <= 0)
#rna = rna[,-z]

#2. cluster cancer types, ie -> each dot should be cancer type
library("FactoMineR")
z1 = which(str_detect(colnames(rna), "ENSG"))
logged_rna = rna
logged_rna[,z1] = log1p(logged_rna[,z1])
#z = which(colnames(logged_rna) %in% rm)
#logged_rna = logged_rna[,-z]
z1 = which(str_detect(colnames(logged_rna), "ENSG"))

#logged - DONE!
#pdf("logged_nonMAD0lncRNAs_cands_PCA_plots_Aug27.pdf", width=12)
#autoplot(prcomp(logged_rna[,z1]), data = logged_rna, colour = 'type')+
#scale_colour_manual(values = mypal)
#dev.off()

#try t-SNE
library(Rtsne)
#iris_unique <- unique(iris) # Remove duplicates
#iris_matrix <- as.matrix(iris_unique[,1:4])
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
