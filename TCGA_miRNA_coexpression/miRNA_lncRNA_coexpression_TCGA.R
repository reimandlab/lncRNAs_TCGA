###miRNA_lncRNA_coexpression_TCGA.R

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:5790)]

library(glmnet)
library(survcomp)
library(caret)
library(stringr)

###---------------------------------------------------------------
###Cands
###---------------------------------------------------------------

#List of canddidates and cox results
cands <- readRDS("36_unique_cands_4cancers_TCGA_Feb6.rds")
z = which(duplicated(cands$gene))
cands = cands[-z,]

###---------------------------------------------------------------
###miRNA expression 
###---------------------------------------------------------------
#1. miRNA expression data 
#OV
ov_mirna = readRDS("OV_miRNA_Expression_data_Feb12.rds")
#PAAD
paad_mirna = readRDS("PAAD_miRNA_Expression_data_Feb12.rds")
#KIRC
kirc_mirna = readRDS("KIRC_miRNA_Expression_data_Feb12.rds")
#LIHC
lihc_mirna = readRDS("LIHC_miRNA_Expression_data_Feb12.rds")
miRNA_data = list(lihc_mirna, ov_mirna, kirc_mirna, paad_mirna)

#2. Expression data 
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
	z = which(order_cancers == cancer)
	mir_dat = miRNA_data[[z]]

	#break up patients into individual rows
	#also only want tumours 
	pats = as.list(unique(colnames(mir_dat)[2:ncol(mir_dat)]))
	rearrange = function(pat){
		z = which(colnames(mir_dat)==pat)
		newdat = mir_dat[,c(1,z), with=FALSE]
		colnames(newdat) = c("probe", "readcount", "RPKM", "post")
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
  	#are there lncRNAs correlated with gene expression?? 
  	df$gene = log1p(df$gene)
  	library("ggExtra")
  	name = fantom$CAT_geneName[which(fantom$gene == lnc)]
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

    return(dat)
}
}






