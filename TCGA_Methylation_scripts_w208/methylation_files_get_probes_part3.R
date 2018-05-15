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
	z = which(is.na(dat$beta))
	if(!(length(z)==0)){
		dat = dat[-z,]
	}

	if(!(dim(dat)[1] == 0)){

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

	#exp_data = rna[which(rna$patient %in% dat$patient), ]
	exp_data = rna[which(rna$type %in% cancer),]
	#assign high or low to each patient in expression file
	z <- which(colnames(exp_data) %in% lnc)
  	if(!(length(z)==0)){
  	df = as.data.frame(exp_data)
  	df <- df[,c(1, z,(ncol(exp_data)-30):ncol(exp_data))]  

	df$median <- ""
 	median2 <- quantile(as.numeric(df[,2]), 0.5)
  	#if(median2 ==0){
    #median2 = mean(as.numeric(df[,2]))
  	#}

  	 if(median2 ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(df[,2] > 0)
		    l2 = which(df[,2] ==0)
		    df$median[l1] = 1
		    df$median[l2] = 0
		    }

	  if(!(median2 ==0)){
		    l1 = which(df[,2] >= median2)
		    l2 = which(df[,2] < median2)
		    df$median[l1] = 1
		    df$median[l2] = 0
		}


  	gene <- colnames(df)[2]
  	df$OS <- as.numeric(df$OS)
  	df$OS.time <- as.numeric(df$OS.time)
  	df$median[df$median ==0] = "Low"
  	df$median[df$median==1] = "High"
  	df$median = factor(df$median, levels=c("Low", "High"))

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	df = merge(df, dat, by=c("patient"))
  	colnames(df)[2] = "geneExp"
  	#is copy number aberation associated with expression? 
  	df$geneExp = log1p(df$geneExp)
  	library("ggExtra")
  	name = probes$lncname[which(probes$ensg == lnc)][1]
	z = which(is.na(df$beta))
	if(length(z) >=1){
		df = df[-z,]
	}
	z =length(table(df$probe))
	
	dim1 = table(df$median)[1]
	dim2 = table(df$median)[2]
	both = (dim1>=1) & (dim2 >=1)

	if(both){

	if(z > 1){
		results_all_probes = as.data.frame(matrix(ncol = 9)) ; colnames(results_all_probes) = c("cancer", "gene", "num_patients", "mean_beta_low", 
			"mean_beta_high", "median_beta_low", "median_beta_high", "wilcoxon_pval", "probe")
	for(k in 1:z){
	new = subset(df, df$probe == unique(df$probe)[k])
	new$geneExp = as.numeric(new$geneExp)
	new$beta = as.numeric(new$beta)

		wilcoxon_pval = wilcox.test(beta ~ median, data =new)$p.value  
		mean_low = mean(new$beta[new$median=="Low"])
		median_low = median(new$beta[new$median=="Low"])
		mean_high = mean(new$beta[new$median=="High"])
		median_high = median(new$beta[new$median=="High"])
		probe = unique(new$probe)
		p =  ggdensity(new, x = "beta", color = "median", fill="median", alpha=0.25, xlab="DNA Methylation",
			title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(df$median=="High"),
				"\npats wLow Exp=" , length(new$median=="Low"), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5)))
		p = p + geom_vline(xintercept = c(0, 0.5, 1), linetype="dotted", 
                color = "blue")
		print(p)

    results = c(cancer, gene, length(unique(new$patient)), mean_low, mean_high, median_low, median_high, wilcoxon_pval, probe)
    names(results) = c("cancer", "gene", "num_patients", "mean_beta_low", "mean_beta_high", "median_beta_low", "median_beta_high", "wilcoxon_pval", "probe")
    results_all_probes = rbind(results_all_probes, results)


    print(lnc)
	}
	results_all_probes = results_all_probes[-1,]
	results = results_all_probes
	}

	if(z==1){

		#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =df)$p.value  
		mean_low = mean(df$beta[df$median=="Low"])
		median_low = median(df$beta[df$median=="Low"])
		mean_high = mean(df$beta[df$median=="High"])
		median_high = median(df$beta[df$median=="High"])
		probe = unique(df$probe)

    results = c(cancer, gene, length(unique(df$patient)), mean_low, mean_high, median_low, median_high, wilcoxon_pval, probe)
    names(results) = c("cancer", "gene", "num_patients", "mean_beta_low", "mean_beta_high", "median_beta_low", "median_beta_high", "wilcoxon_pval", "probe")

   		p =  ggdensity(df, x = "beta", color = "median", fill="median", alpha=0.25, xlab="DNA Methylation",
			title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(df$median=="High"),
				"\npats wLow Exp=" , length(df$median=="Low"), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5)))
   		p = p + geom_vline(xintercept = c(0, 0.5, 1), linetype="dotted", 
                color = "blue")
		print(p)
    print(lnc)
	}
	}

	if(!(both)){
		results = as.data.frame(matrix(ncol = 9)) ; colnames(results) = c("cancer", "gene", "num_patients", "mean_beta_low", 
			"mean_beta_high", "median_beta_low", "median_beta_high", "wilcoxon_pval", "probe")
	}

    return(results)
}
}
}
}

pdf("candidate_lncRNAs_methylation_versus_Expression_May15_only_NOFDR_candidates.pdf")
cands = filter(cands, AnalysisType == "noFDR")
genes = as.list(unique(cands$gene[which(cands$gene %in% probes$ensg)])) #110/190 have methylation probes overlapping them 
lnc_meth_cancer_data = llply(genes, get_data, .progress="text")
dev.off()


#########################################################################################################################
#SUMMARIZE#
#########################################################################################################################

lnc_meth_cancer_data2 = as.data.frame(do.call("rbind", lnc_meth_cancer_data))
z = which(is.na(lnc_meth_cancer_data2[,1]))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[-z,]

lnc_meth_cancer_data2 = as.data.table(lnc_meth_cancer_data2)

lnc_meth_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_meth_cancer_data2$wilcoxon_pval))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[order(wilcoxon_pval)]
lnc_meth_cancer_data2$fdr = p.adjust(lnc_meth_cancer_data2$wilcoxon_pval, method="fdr")

sig_diff = filter(lnc_meth_cancer_data2, fdr <=0.05)
sig_diff = as.data.frame(sig_diff)

sig_diff$mean_beta_high = as.numeric(sig_diff$mean_beta_high)
sig_diff$mean_beta_low = as.numeric(sig_diff$mean_beta_low)
sig_diff$median_beta_high = as.numeric(sig_diff$median_beta_high)
sig_diff$median_beta_low = as.numeric(sig_diff$median_beta_low)


sig_diff$mean_diff = sig_diff$mean_beta_low/sig_diff$mean_beta_high
sig_diff$median_diff = sig_diff$median_beta_low/sig_diff$median_beta_high

gene_ids = cands[,c(1,9)]
sig_diff = merge(sig_diff, gene_ids, by="gene")

#sig_diff = sig_diff[order(-median_diff)]

#check which probes have higher mean and median beta vlaue in the low expression group and if they are also
#in the promoter 
sig_diff$meth_status = ""

check_meth_exp = function(row){
	mean_diff = row[[11]]
	if(mean_diff > 1.2){
		stat = "LowExpHighMeth"
	}
	if(mean_diff < 0.8){
		stat = "HighExpHighMeth"
	}
	if((mean_diff >=0.8) & (mean_diff <= 1.2)){
		stat = "noDifferent"
	}
	return(stat)
}

sig_diff$meth_status = apply(sig_diff, 1, check_meth_exp)

#summary how many lncRNAs have methylation correlation 

sig_diff$mean_diff = log2(sig_diff$mean_diff)
sig_diff$median_diff = log2(sig_diff$median_diff)

sig_diff$combo = paste(sig_diff$gene, sig_diff$probe, sig_diff$cancer, sep="_")
sig_diff = sig_diff[!duplicated(sig_diff), ]

#is probe in promoter or exon? 
colnames(probes)[11] = "gene"
probes = probes[,c(1:4, 6:9, 11:14)]
colnames(sig_diff)[9] = "cgid"
sig_diff = merge(sig_diff, probes, by=c("gene", "cgid"))

check_promoter = function(row){
	strand = row$lncstrand
	lnc_prom_start = row$lncstart-2000
	lnc_prom_end = row$lncstart+2000
	cpg = row$cpg_start+1

	#promoter
	if((cpg >lnc_prom_start) & (cpg<lnc_prom_end)){
		region="promoter"
	}

	if(!((cpg >lnc_prom_start) & (cpg<lnc_prom_end))){
		region="body"
	}
	return(region)
}
sig_diff$region = apply(sig_diff, 1, check_promoter)

pdf("summary_methylation_of_candidates_noFDR_candidates_only.pdf", width=18)
ggbarplot(sig_diff, x = "combo", y = "mean_diff",
          #facet.by = "cancer",
          fill = "region",           # change fill color by mpg_level
          color = "region",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          #sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Mean Beta Value",
          xlab = FALSE,
          legend.title = "Region bound by CpG"
          ) + theme(axis.text.x = element_text(size=6, angle=90))

dev.off()

write.table(sig_diff, file="methylation_analysis_of_candidate_lncRNAs_wilcoxon_resutls_May15.txt", quote=F, row.names=F, sep="\t")





































###plot 
#z > 1
	sp = ggscatter(new, font.x = c(15, "plain", "black"), main= paste(name, new$canc[1], unique(new$probe)[1]), 
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
		x = "beta", y = "geneExp", 
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

	yplot <- ggboxplot(new, x = "median", y = "geneExp", font.x = c(15, "plain", "black"),
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


###previous plots
sp = ggscatter(df, font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), 
		x = "beta", y = "geneExp", 
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   shape = "median", ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", xlab="Beta Value") + stat_cor() 
	
	
	xplot = ggboxplot(df, x="median", y="beta", orientation="horizontal", palette="jco", fill="median", main= paste(name, cancer, "Methylation vs Exp", "n=", length(unique(df$patient)), df$probe[1]),
		legend.title = "Expression Tag", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),order=(c("Low", "High")),  ggtheme = theme_light(),
          xlab="Expression", ylab="Beta Value")
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5) 
	
	yplot <- ggboxplot(df, x = "median", y = "geneExp", font.x = c(15, "plain", "black"),
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




