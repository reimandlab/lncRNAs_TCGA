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
colScale_full+theme_classic()+
geom_text_repel(data = subset(layout, !(label == "no")))
dev.off()


write.table(layout, file="/u/kisaev/Jan2021/umap_29_cancers_data_for_figure.txt", quote=F, row.names=F, col.names=T, sep="\t")
