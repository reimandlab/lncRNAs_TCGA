library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggthemes)


#1. See if any candidates have CNAs
lncswcnas = fread("fantom_lncrnas_wTCGA_CNAs_4cancers.bed")
lncswcnas = as.data.frame(lncswcnas)

#2. cands 
cands = readRDS("36_unique_cands_4cancers_TCGA_Feb6.rds")
colnames(lncswcnas)[4] = "gene"

#keep gene-cancer combinations, don't really care right now if gene has CNA
#for a different cancer where it's not a candidate 
genes = as.list(unique(cands$gene[which(cands$gene %in% lncswcnas$gene)])) #29/33 have CNAs overlapping them 
colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "gene", "name", "Chromosome" , "Start" , 
	"End", "Num_Probes" , "Segment_Mean", "canc", "rm", "patient")
lncswcnas$rm = NULL

#3. Expression data 
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
	dat = dplyr::filter(lncswcnas, canc == cancer, gene == lnc)
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

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	df = merge(df, dat, by=c("patient", "canc"))
  	colnames(df)[3] = "gene"
  	#is copy number aberation associated with expression? 
  	df$gene = log1p(df$gene)
  	df$median = factor(df$median, levels=c("Low", "High"))
  	library("ggExtra")
	sp = ggscatter(df, 
		x = "Segment_Mean", y = "gene",
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), xlab="Segment Mean SCNA") + stat_cor() 
	
	xplot = ggboxplot(df, main= paste(df$name[1], df$canc[1], "CNA vs Exp", "n=", length(unique(df$patient))),
		x = "median", y = "Segment_Mean", legend.title = "Expression Tag", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
                  fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Segment Mean SCNA")+rotate()
	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5)
	yplot <- ggboxplot(df, x = "median", y = "gene", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    yplot = yplot + rremove("legend")
    sp = sp + rremove("legend")
    #p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, nrow =2,align = "hv", 
    #      rel_widths = c(2, 1), rel_heights = c(1, 2))
   
    plots <- align_plots(xplot, sp, align = 'v', axis = 'l')
    bottom_row <- plot_grid(plots[[2]], yplot, labels = c('B', 'C'), align = 'h', rel_widths = c(2.85, 1))
    p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    print(p)
    print(lnc)
    return(dat)
}
}
pdf("candidate_lncRNAs_CNA_versus_Expression_Feb7.pdf", height=5.5, width=8.3)
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
dev.off()

#how to make sense of this? 
#plot segment? and where within it lies lncRNA? 

