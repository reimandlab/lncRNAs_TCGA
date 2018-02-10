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
colnames(lncswcnas)[4] = "gene"

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
kirc = readRDS("KIRC_methylation_data_lncs_cands.rds")
#LIHC
lihc = readRDS("LIHC_methylation_data_lncs_cands.rds")
#OV
ov = readRDS("OV_methylation_data_lncs_cands.rds")
#PAAD
paad = readRDS("PAAD_methylation_data_lncs_cands.rds")
methylation_data = list(lihc, ov, kirc, paad)

#4. Expression data 
liver = readRDS("LIHC_269_pats_RNASeq_data_Feb7.rds")
liver$canc = "liver"
ov = readRDS("OV_295_pats_RNASeq_data_Feb7.rds")
ov$canc = "ovary"
kidney = readRDS("KIRC_463_pats_RNASeq_data_Feb7.rds")
kidney$canc = "kidney"
pancreas = readRDS("PAAD_169_pats_RNASeq_data_Feb7.rds")
pancreas$canc = "pancreas"
expression_data = list(liver, ov, kidney, pancreas)
order_cancers = c("liver", "ovary", "kidney", "pancreas")

get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
	dat = dplyr::filter(probes, canc == cancer, ensg == lnc)
	if(!(dim(dat)[1]==0)){
	z = which(order_cancers == cancer)
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

	z = which(order_cancers == cancer)
	exp_data = expression_data[[z]]
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
	sp = ggscatter(df, 
		x = "beta", y = "gene", facet.by = "probe", 
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   shape = "median", ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", xlab="Beta Value") + stat_cor() 
	
	xplot = ggboxplot(df, main= paste(name, df$canc[1], "Methylation vs Exp", "n=", length(unique(df$patient))),
		x = "median", y = "beta", legend.title = "Expression Tag",facet.by = "probe",
                   color = "median", fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Beta Value")+rotate()
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5)
	yplot <- ggboxplot(df, x = "median", y = "gene", 
                   color = "median", fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    yplot = yplot + rremove("legend")
    sp = sp + rremove("legend")
    plots <- align_plots(xplot, sp, align = 'v', axis = 'l')
    bottom_row <- plot_grid(plots[[2]], yplot, labels = c('B', 'C'), align = 'h', rel_widths = c(4.3, 1))
	p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    #p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, nrow =2,align = "hv", 
    #      rel_widths = c(2, 1), rel_heights = c(1, 2))
    print(p)
    print(lnc)
	}

	if(z==1){
	sp = ggscatter(df, 
		x = "beta", y = "gene", 
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   shape = "median", ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", xlab="Beta Value") + stat_cor() 
	
	xplot = ggboxplot(df, main= paste(name, df$canc[1], "Methylation vs Exp", "n=", length(unique(df$patient))),
		x = "median", y = "beta", legend.title = "Expression Tag", 
                   color = "median", fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Beta Value")+rotate()
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5)
	xplot = ggpar(xplot, legend="top")
	yplot <- ggboxplot(df, x = "median", y = "gene", 
                   color = "median", fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    yplot = yplot + rremove("legend")
    sp = sp + rremove("legend")
    plots <- align_plots(xplot, sp, align = 'v', axis="l")
    bottom_row <- plot_grid(plots[[2]], yplot, labels = c('A', 'B'), align = 'h', axis="l", rel_widths = c(4,1))
	p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    
    print(p)
    print(lnc)	
	}

    return(dat)
}
}
}

pdf("candidate_lncRNAs_methylation_versus_Expression_Feb9.pdf", width=13, height =11)
genes = as.list(unique(cands$gene[which(cands$gene %in% probes$ensg)])) #23 have methylation probes overlapping them 
lnc_meth_cancer_data = llply(genes, get_data, .progress="text")
dev.off()




