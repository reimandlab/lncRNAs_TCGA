source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

mypal = c("#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

mypal = "purple"

library(M3C)
library(umap)

lgg=readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")
lgg = unique(lgg[,c("patient", "IDH.status")])


lgg_rna = as.data.table(filter(rna, type =="LGG"))
lgg_pfi = unique(lgg_rna[,c("patient", "PFI", "PFI.time")])
lgg_pfi = merge(lgg, lgg_pfi, by="patient")
saveRDS(lgg_pfi, "/u/kisaev/TCGA_LGG_PFS_time.rds")


kirc = readRDS("KIRC_maf.rds")
kirc = unique(kirc[,c("Hugo_Symbol", "Tumor_Sample_Barcode")])
kirc$Tumor_Sample_Barcode = sapply(kirc$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
colnames(kirc)[2] = "patient"
kirc = as.data.table(filter(kirc, Hugo_Symbol == "BAP1"))

kirc_sub = readRDS("/u/kisaev/TCGA_kirc_wsubtype_info_biolinks.rds")
kirc_sub = unique(kirc_sub[,c("patient", "mRNA_cluster")])

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

dim(rna)
rownames(rna) = rna$patient

kirc_rna = as.data.table(filter(rna, type=="KIRC"))


#KIRC
	z1 = which(str_detect(colnames(kirc_rna), "ENSG"))
	z2 = which(colnames(kirc_rna) %in% c("type", "patient"))
	kirc_rna = as.data.frame(kirc_rna)
	kirc_rna = kirc_rna[,c(z1, z2)]
	kirc_rna = merge(kirc_rna, kirc_sub, by="patient")

	cancer.labels = as.character(kirc_rna$mRNA_cluster)
	kirc_rna$mRNA_cluster = NULL

	#1. remove those not expressed at all
	kirc_rna = as.data.table(kirc_rna)
	z1 = which(str_detect(colnames(kirc_rna), "ENSG"))
	kirc_rna[,z1] = log1p(kirc_rna[,..z1])
	
	set.seed(100)
	df = kirc_rna 
	df$type = NULL
	df$patient = NULL
	embedding = umap::umap(df)
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

pdf("/u/kisaev/KIRC_mRNA_subtypes_umap_plot.pdf")
	g = ggplot(layout,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col),alpha=0.5, stroke=0) + 
	scale_colour_manual(values=c("purple", "orange", "red", "green", "black","blue"))+
		geom_text_repel(data = subset(layout, !(label == "no")))
	print(g)
dev.off()

