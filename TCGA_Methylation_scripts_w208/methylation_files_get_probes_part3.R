###now have methylation files for each cancer 
###with probes that overlap any of the lncRNA candidates (might not be the cancer specific candidate)

library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggthemes)

tss_codes = read.csv("TCGA_TissueSourceSite_Codes2017.csv")
source_codes = source = read.csv("TCGA_sample_codes.csv")

#1. cands 
cands = readRDS("36_unique_cands_4cancers_TCGA_Feb6.rds")
#candidate lncrnas 
cands = readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
colnames(cands)[7] = "canc"
#colnames(cands)[3] = "canc"

#2. lncRNAs that intersected with probes
probes = fread("fantom_lncrnas_mapped_to450_probes.bed")
colnames(probes) = c("probe_chr", "cpg_start", "cpg_end", "cgid", "rand", "cpgstrand", "lncchr", "lncstart", "lncend", "rand2", "ensg", "lncstrand", "lncname")

##add cancer type to which that gene belong to 
add_canc = function(probeid){
	lnc = probes$ensg[which(probes$cgid == probeid)]
	canc = cands$canc[which(cands$gene == lnc)][1]
	return(canc)
}
probes$canc = llply(probes$cgid, add_canc, .progress="text")

#3. methylation files
#KIRC
KIRC = readRDS("KIRC_methylation_data_lncs_cands.rds")
#LIHC
LIHC = readRDS("LIHC_methylation_data_lncs_cands.rds")
#OV
OV = readRDS("OV_methylation_data_lncs_cands.rds")
#PAAD
PAAD = readRDS("PAAD_methylation_data_lncs_cands.rds")
#LUAD
LUAD = readRDS("LUADmethylation_data_lncs_cands.rds")
#BRCA
BRCA = readRDS("BRCAmethylation_data_lncs_cands.rds")
#
ACC = readRDS(file="ACCmethylation_data_lncs_cands.rds")

ESCA = readRDS(file="ESCAmethylation_data_lncs_cands.rds")

GBM = readRDS(file="GBMmethylation_data_lncs_cands.rds")

HNSC = readRDS(file="HNSCmethylation_data_lncs_cands.rds")

LGG = readRDS(file="LGGmethylation_data_lncs_cands.rds")

MESO = readRDS(file="MESOmethylation_data_lncs_cands.rds")

READ = readRDS(file="READmethylation_data_lncs_cands.rds")

THCA = readRDS(file="THCAmethylation_data_lncs_cands.rds")

UVM = readRDS(file="UVMmethylation_data_lncs_cands.rds")

CESC = readRDS(file="CESCmethylation_data_lncs_cands.rds")

KIRP = readRDS(file="KIRPmethylation_data_lncs_cands.rds")

STAD = readRDS(file="STADmethylation_data_lncs_cands.rds")

LUSC = readRDS(file="LUSCmethylation_data_lncs_cands.rds")

BLCA =  readRDS(file="BLCAmethylation_data_lncs_cands.rds")

COAD = readRDS(file="COADmethylation_data_lncs_cands.rds")

UCEC = readRDS(file="UCECmethylation_data_lncs_cands.rds")

SARC = readRDS(file="SARCmethylation_data_lncs_cands.rds")

#combine all 
methylation_data = list(KIRC,
LIHC,
OV,
PAAD,
LUAD,
BRCA,
BRCA,
ACC,
ESCA,
GBM,
HNSC,
LGG,
MESO,
READ,
THCA,
UVM,
CESC,
KIRP,
STAD,
LUSC,
BLCA,
COAD,
UCEC,
SARC)

cancers_order = c("KIRC",
"LIHC",
"OV",
"PAAD",
"LUAD",
"BRCA",
"BRCA",
"ACC",
"ESCA",
"GBM",
"HNSC",
"LGG",
"MESO",
"READ",
"THCA",
"UVM",
"CESC",
"KIRP",
"STAD",
"LUSC",
"BLCA",
"COAD",
"UCEC",
"SARC")

cands$canc[cands$canc=="Ovarian serous cystadenocarcinoma"] = "OV"
cands$canc[cands$canc=="Liver hepatocellular carcinoma"] = "LIHC"
cands$canc[cands$canc=="Kidney renal clear cell carcinoma"] = "KIRC"
cands$canc[cands$canc=="Lung squamous cell carcinoma"] = "LUSC"
cands$canc[cands$canc=="Breast invasive carcinoma"] = "BRCA"
cands$canc[cands$canc=="Lung adenocarcinoma"] = "LUAD"
cands$canc[cands$canc=="Pancreatic adenocarcinoma"] = "PAAD"
cands$canc[cands$canc=="Adrenocortical carcinoma"] = "ACC"
cands$canc[cands$canc=="Bladder Urothelial Carcinoma"] = "BLCA"
cands$canc[cands$canc=="Stomach adenocarcinoma"] = "STAD"
cands$canc[cands$canc=="Head and Neck squamous cell carcinoma"] = "HNSC"
cands$canc[cands$canc=="Brain Lower Grade Glioma"] = "LGG"
cands$canc[cands$canc=="Sarcoma"] = "SARC"
cands$canc[cands$canc=="Kidney renal papillary cell carcinoma"] = "KIRP"
cands$canc[cands$canc=="Mesothelioma"] = "MESO"
cands$canc[cands$canc=="Uterine Corpus Endometrial Carcinoma"] = "UCEC"
cands$canc[cands$canc=="Uveal Melanoma"] = "UVM"
cands$canc[cands$canc=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] = "CESC"
cands$canc[cands$canc=="Colon adenocarcinoma"] = "COAD"
cands$canc[cands$canc=="Rectum adenocarcinoma"] = "READ"
cands$canc[cands$canc=="Thyroid carcinoma"] = "THCA"
cands$canc[cands$canc=="Glioblastoma multiforme"] = "GBM"
cands$canc[cands$canc=="Esophageal carcinoma"] = "ESCA"


probes$canc[probes$canc=="Ovarian serous cystadenocarcinoma"] = "OV"
probes$canc[probes$canc=="Liver hepatocellular carcinoma"] = "LIHC"
probes$canc[probes$canc=="Kidney renal clear cell carcinoma"] = "KIRC"
probes$canc[probes$canc=="Lung squamous cell carcinoma"] = "LUSC"
probes$canc[probes$canc=="Breast invasive carcinoma"] = "BRCA"
probes$canc[probes$canc=="Lung adenocarcinoma"] = "LUAD"
probes$canc[probes$canc=="Pancreatic adenocarcinoma"] = "PAAD"
probes$canc[probes$canc=="Adrenocortical carcinoma"] = "ACC"
probes$canc[probes$canc=="Bladder Urothelial Carcinoma"] = "BLCA"
probes$canc[probes$canc=="Stomach adenocarcinoma"] = "STAD"
probes$canc[probes$canc=="Head and Neck squamous cell carcinoma"] = "HNSC"
probes$canc[probes$canc=="Brain Lower Grade Glioma"] = "LGG"
probes$canc[probes$canc=="Sarcoma"] = "SARC"
probes$canc[probes$canc=="Kidney renal papillary cell carcinoma"] = "KIRP"
probes$canc[probes$canc=="Mesothelioma"] = "MESO"
probes$canc[probes$canc=="Uterine Corpus Endometrial Carcinoma"] = "UCEC"
probes$canc[probes$canc=="Uveal Melanoma"] = "UVM"
probes$canc[probes$canc=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] = "CESC"
probes$canc[probes$canc=="Colon adenocarcinoma"] = "COAD"
probes$canc[probes$canc=="Rectum adenocarcinoma"] = "READ"
probes$canc[probes$canc=="Thyroid carcinoma"] = "THCA"
probes$canc[probes$canc=="Glioblastoma multiforme"] = "GBM"
probes$canc[probes$canc=="Esophageal carcinoma"] = "ESCA"



#4. Expression data 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
	dat = dplyr::filter(probes, canc == cancer, ensg == lnc)
	if(!(dim(dat)[1]==0)){
	z = which(cancers_order == cancer)
	meth_dat = methylation_data[[z]]
	z = which(meth_dat$probe %in% dat$cgid)
	meth_dat = meth_dat[z,]

	#break up patients into individual rows
	#also only want tumours 
	pats = as.list(unique(colnames(meth_dat)[2:ncol(meth_dat)]))
	rearrange = function(pat){
		z = which(colnames(meth_dat)==pat)
		newdat = meth_dat[,c(1,z), with=FALSE]
		colnames(newdat) = c("probe", "beta", "gene", "post", "coord")
		newdat$patient = pat
		return(newdat)
	}
	newdats = llply(pats, rearrange, .progress="text")
	dat <- ldply(newdats, data.table)
	dat$gene = lnc

	dat$source = ""
	get_source = function(id){
		source = unlist(strsplit(as.character(id), '-'))[4]
		source = substr(source, 1, 2)
		return(source)
		}
	dat$source = llply(dat$patient, get_source, .progress="text")
	#only keep primary solid tumour samples 
	dat = dplyr::filter(dat, source =="01")
	clean_tcga_id = function(id){
		s1 = unlist(strsplit(id, '-'))[1]
		s2 = unlist(strsplit(id, '-'))[2]
		s3 = unlist(strsplit(id, '-'))[3]
		return(paste(s1, s2, s3, sep="-"))
		}	
	dat$patient = llply(dat$patient, clean_tcga_id, .progress="text")

	exp_data = rna[which(rna$patient %in% dat$patient), ]
	#assign high or low to each patient in expression file
	z <- which(colnames(exp_data) %in% lnc)
  	if(!(length(z)==0)){
  	df = as.data.frame(exp_data)
  	df <- df[,c(z,(ncol(exp_data)-4):ncol(exp_data))]  

	df$median <- ""
 	median2 <- quantile(as.numeric(df[,1]), 0.5)
  	if(median2 ==0){
    median2 = mean(as.numeric(df[,1]))
  	}

  	#median2 <- median(df[,1])
  	for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
  	gene <- colnames(df)[1]
  	df$status[df$status=="Alive"] <- 0
  	df$status[df$status=="Dead"] <- 1
  	df$status <- as.numeric(df$status)
  	df$time <- as.numeric(df$time)
  	df$median[df$median ==0] = "Low"
  	df$median[df$median==1] = "High"
  	df$median = factor(df$median, levels=c("Low", "High"))

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	df = merge(df, dat, by=c("patient"))
  	colnames(df)[2] = "gene"
  	#is copy number aberation associated with expression? 
  	df$gene = log1p(df$gene)
  	library("ggExtra")
  	name = probes$lncname[which(probes$ensg == lnc)][1]
	z = which(is.na(df$beta))
	if(length(z) >=1){
		df = df[-z,]
	}
	z =length(table(df$probe))
	if(z > 1){
	for(k in 1:z){
	new = subset(df, df$probe == unique(df$probe)[k])
	new$gene = as.numeric(new$gene)
	new$beta = as.numeric(new$beta)

	sp = ggscatter(new, font.x = c(15, "plain", "black"), main= paste(name, new$canc[1], unique(new$probe)[1]), 
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
		x = "beta", y = "gene", 
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   shape = "median", ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", xlab="Beta Value") + stat_cor() 
	
	print(sp)

	xplot = ggboxplot(new, main= paste(name, new$canc[1], "Methylation vs Exp", "n=", length(unique(new$patient))),
		x = "median", y = "beta", legend.title = "Expression Tag", font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
                   fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Beta Value")+rotate()
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5) 
	print(xplot)

	yplot <- ggboxplot(new, x = "median", y = "gene", font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
                   fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5) 
    yplot = yplot + rremove("legend")
    print(yplot)

    #sp = sp + rremove("legend")
    #plots <- align_plots(xplot, sp, align = 'v', axis = 'l')
    #bottom_row <- plot_grid(plots[[2]], yplot, labels = c('B', 'C'), align = 'h', rel_widths = c(4.3, 1))
	#p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    #print(p)
    print(lnc)
	}
	}

	if(z==1){
	sp = ggscatter(df, font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
		x = "beta", y = "gene", 
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   shape = "median", ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", xlab="Beta Value") + stat_cor() 
	
	
	xplot = ggboxplot(df, x="median", y="beta", orientation="horizontal", palette="jco", fill="median", main= paste(name, df$canc[1], "Methylation vs Exp", "n=", length(unique(df$patient)), df$probe[1]),
		legend.title = "Expression Tag", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),order=(c("Low", "High")),  ggtheme = theme_light(),
          xlab="Expression", ylab="Beta Value")
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5) 
	
	yplot <- ggboxplot(df, x = "median", y = "gene", font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
                   fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5) 
    yplot = yplot + rremove("legend")

    sp = sp + rremove("legend")
    plots <- align_plots(xplot, sp, align = 'v', axis = 'l')
    bottom_row <- plot_grid(plots[[2]], yplot, labels = c('B', 'C'), align = 'h', rel_widths = c(2.85, 1))
	p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    #p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, nrow =2,align = "hv", 
    #      rel_widths = c(2, 1), rel_heights = c(1, 2))
    print(p)
    print(lnc)
	}

    return(dat)
}
}
}

pdf("candidate_lncRNAs_methylation_versus_Expression_April3.pdf", height=5, width=6)
genes = as.list(unique(cands$gene[which(cands$gene %in% probes$ensg)])) #110/190 have methylation probes overlapping them 
lnc_meth_cancer_data = llply(genes, get_data, .progress="text")
dev.off()






