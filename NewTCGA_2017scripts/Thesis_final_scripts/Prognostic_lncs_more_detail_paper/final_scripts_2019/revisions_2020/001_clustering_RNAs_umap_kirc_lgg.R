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

kirc = readRDS("KIRC_maf.rds")
kirc = unique(kirc[,c("Hugo_Symbol", "Tumor_Sample_Barcode")])
kirc$Tumor_Sample_Barcode = sapply(kirc$Tumor_Sample_Barcode, function(x){paste(unlist(strsplit(x, "-"))[1:3], collapse="-")})
colnames(kirc)[2] = "patient"
kirc = as.data.table(filter(kirc, Hugo_Symbol == "BAP1"))

kirc_sub = readRDS("/u/kisaev/TCGA_kirc_wsubtype_info_biolinks.rds")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

dim(rna)
rownames(rna) = rna$patient

kirc_rna = as.data.table(filter(rna, type=="KIRC"))
#get patient risk
surv = unique(kirc_rna[,c("patient", "OS", "OS.time", "histological_grade")])
surv = as.data.table(filter(surv, !(histological_grade %in% c("unknown", "GX"))))
risks = coxph(Surv(OS.time, OS) ~ histological_grade, data = surv)
relRisk <- predict(risks, surv, type="risk")   # relative risk
surv$risk = relRisk
risks = coxph(Surv(OS.time, OS) ~ risk, data = surv)
med = median(surv$risk)
surv$risk_med = ""
surv$risk_med[surv$risk >=med] = "high_risk"
surv$risk_med[surv$risk < med] = "low_risk"


z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% c("type", "patient", "OS", "OS.time"))
rna = as.data.frame(rna)
rna = rna[,c(z1, z2)]


#KIRC
	z1 = which(str_detect(colnames(kirc_rna), "ENSG"))
	z2 = which(colnames(kirc_rna) %in% c("type", "patient"))
	kirc_rna = as.data.frame(kirc_rna)
	kirc_rna = kirc_rna[,c(z1, z2)]
	kirc_rna = merge(kirc_rna, surv, by="patient")

	cancer.labels = kirc_rna$risk_med
	kirc_rna$risk_med = NULL
	kirc_rna$risk = NULL
	kirc_rna$histological_grade = NULL
	kirc_rna$OS = NULL
	kirc_rna$OS.time = NULL

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

pdf("/u/kisaev/KIRC_grade_risk_umap_plot.pdf")
	g = ggplot(layout,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col),alpha=0.5, stroke=0) +
	scale_colour_manual(values=c("purple", "orange", "red", "green", "black","blue"))+
		geom_text_repel(data = subset(layout, !(label == "no")))
	print(g)
dev.off()

get_kirc_umap = function(var){

	z1 = which(str_detect(colnames(kirc_rna), "ENSG"))
	z2 = which(colnames(kirc_rna) %in% c("type", "patient", var))
	kirc_rna = as.data.frame(kirc_rna)
	kirc_rna = kirc_rna[,c(z1, z2)]
	cancer.labels = kirc_rna[,which(colnames(kirc_rna)==var)]
	kirc_rna = kirc_rna[,-which(colnames(kirc_rna)==var)]
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

	g = ggplot(layout,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col),alpha=0.5, stroke=0) +
	scale_colour_manual(values=c("purple", "orange", "red", "green", "black","blue"))+
		geom_text_repel(data = subset(layout, !(label == "no")))
	print(g)


}

pdf("/u/kisaev/KIRC_grade_umap_plot.pdf")
get_kirc_umap("histological_grade")
dev.off()


get_umap = function(canc){

	canc_dat = as.data.table(filter(rna, type==canc))

	if(canc == "LGG"){
		canc_dat = as.data.table(merge(lgg, canc_dat, by="patient"))
		cancer.labels = canc_dat$IDH.status
		canc_dat$IDH.status = NULL
	}

	if(canc == "KIRC"){
		canc_dat$vhl = ""
		z = which(canc_dat$patient %in% kirc$patient)
		canc_dat$vhl[z] = "VHL_mut"
		canc_dat$vhl[-z] = "VHL_wt"
		cancer.labels = canc_dat$vhl
		canc_dat$vhl = NULL
	}

	#1. remove those not expressed at all
	z1 = which(str_detect(colnames(canc_dat), "ENSG"))
	canc_dat[,z1] = log1p(canc_dat[,..z1])

	set.seed(100)

	if(!(canc %in% c("LGG", "KIRC"))){
		cancer.labels = canc_dat$type
	}
	df = canc_dat
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

	g = ggplot(layout,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col),alpha=0.5, stroke=0) + scale_colour_manual(values=c("purple", "orange"))+
		geom_text_repel(data = subset(layout, !(label == "no")))
	print(g)
}

pdf("/u/kisaev/LGG_KIRC_umap_plots.pdf")
get_umap("KIRC")
get_umap("LGG")
dev.off()
