library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggthemes)


#1. See if any candidates have CNAs
lncswcnas = fread("fantom_lncrnas_wTCGA_CNAs_23cancers.bed")
lncswcnas = as.data.frame(lncswcnas)

#2. cands 
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
colnames(cands)[7] = "canc"
colnames(lncswcnas)[4] = "gene"
colnames(lncswcnas)[11] = "canc"

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

lncswcnas$canc[lncswcnas$canc=="ovary"] = "OV"
lncswcnas$canc[lncswcnas$canc=="liver"] = "LIHC"
lncswcnas$canc[lncswcnas$canc=="kidney"] = "KIRC"
lncswcnas$canc[lncswcnas$canc=="brca"] = "BRCA"
lncswcnas$canc[lncswcnas$canc=="luad"] = "LUAD"
lncswcnas$canc[lncswcnas$canc=="pancreas"] = "PAAD"

#turn cancers to capitals 
lncswcnas = lncswcnas %>% 
      mutate(canc = toupper(canc))

#keep gene-cancer combinations, don't really care right now if gene has CNA
#for a different cancer where it's not a candidate 
genes = as.list(unique(as.character(cands$gene[which(cands$gene %in% lncswcnas$gene)]))) #170/190 have CNAs overlapping them 
colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "gene", "name", "Chromosome" , "Start" , 
	"End", "Num_Probes" , "Segment_Mean", "canc", "rm", "patient")
lncswcnas$rm = NULL

  ##3. Expression data 
  #lihc = readRDS("LIHC_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
  #lihc$canc = "lihc"
  #ov = readRDS("OV_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
  #ov$canc = "ovary"
  #kirc = readRDS("KIRC_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
  #kirc$canc = "kirc"
  #paad = readRDS("PAAD_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
  #paad$canc = "paad"
  #expression_data = list(lihc, ov, kirc, paad)
  #order_cancers = c("lihc", "ov", "kirc", "paad")


rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
	dat = dplyr::filter(lncswcnas, canc == cancer, gene == lnc)
	dat$canc = NULL

	exp_data = subset(rna, type == cancer)
	#assign high or low to each patient in expression file
	z <- which(colnames(exp_data) %in% lnc)
  	if(!(length(z)==0)){
  	df = as.data.frame(exp_data)
  	df <- df[,c(1, z,(ncol(exp_data)-32):ncol(exp_data))]  

	df$median <- ""
 	median2 <- quantile(as.numeric(df[,2]), 0.5)
  	if(median2 ==0){
    median2 = mean(as.numeric(df[,2]))
  	}

  	#median2 <- median(df[,1])
  	for(y in 1:nrow(df)){
    genexp <- df[y,2]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
  	gene <- colnames(df)[2]
  	df$OS <- as.numeric(df$OS)
  	df$OS.time <- as.numeric(df$OS.time)
  	df$median[df$median ==0] = "Low"
  	df$median[df$median==1] = "High"

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	df = merge(df, dat, by=c("patient"))
  	colnames(df)[2] = "geneexp"
  	#is copy number aberation associated with expression? 
  	df$geneexp = log1p(df$geneexp)
  	df$median = factor(df$median, levels=c("Low", "High"))

    df = df[,c(1:2, 3, 36:46)]
    df$V1 = NULL
    df$cna_status = ""
    get_cna_stat = function(seg_mean){
       dup = seg_mean >=0.1
       del = seg_mean < (-0.1)
       if(dup){
        return("DUP")
       }        
       if(del){
        return("DEL")
       } 
       if((!(dup)) & (!(del))){
        return("unreliable")
       }
    }
    df$cna_status = as.character(llply(df$Segment_Mean, get_cna_stat))

    df$median = as.character(df$median)
    df$exp_cna_match = ""
    get_cna_exp = function(row){
      exp = row[[3]]
      cna = row[[14]]
      if(((exp == "High") & (cna == "DUP"))){
        return("HighExp_Duplicated")
      }
        if(((exp == "Low") & (cna == "DEL"))){
        return("LowExp_Deleted")
      }
        else{
          return("DontMatch")
        }
    }
    df$exp_cna_match = as.character(apply(df, 1, get_cna_exp))
    #get_wilcoxon_pval 

    wilcoxon_pval = wilcox.test(Segment_Mean ~ median, data =df)$p.value  

    results = c(cancer, unique(df$gene), unique(df$name), length(unique(df$patient)), length(which(df$exp_cna_match == "HighExp_Duplicated")), 
      length(which(df$exp_cna_match == "LowExp_Deleted")), wilcoxon_pval)
    names(results) = c("cancer", "gene", "name", "num_patients", "numHighCNAmatch", "numLowCNAmatch", "wilcoxon_pval")


  	library("ggExtra")
	    sp = ggscatter(df, main = paste(df$gene[1], df$type[1]), 
	   	x = "Segment_Mean", y = "geneexp",
               color = "median", palette = "jco",
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), xlab="Segment Mean SCNA") + stat_cor() 
    	print(sp)
	     xplot = ggboxplot(df, main= paste(df$name[1], df$canc[1], "CNA vs Exp", "n=", length(unique(df$patient))),
		    x = "median", y = "Segment_Mean", legend.title = "Expression Tag", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
                  fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Segment Mean SCNA")+rotate()
    	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    	print(xplot)
       yplot <- ggboxplot(df, x = "median", y = "geneexp", main = df$name[1], font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), fill = "median", palette = "jco", order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	     yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    
      print(yplot)
      print(lnc)
    return(results)
}
}
pdf("candidate_lncRNAs_CNA_versus_Expression_May9.pdf", height=5.5, width=8.3)
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
dev.off()


lnc_cna_cancer_data2 = as.data.frame(do.call("rbind", lnc_cna_cancer_data))
lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)
lnc_cna_cancer_data2 = lnc_cna_cancer_data2[order(wilcoxon_pval, -numHighCNAmatch, -numLowCNAmatch)]
lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_cna_cancer_data2$wilcoxon_pval))
sig_diff = filter(lnc_cna_cancer_data2, wilcoxon_pval <=0.05)
sig_diff = as.data.frame(sig_diff)
sig_diff$no_match = ""
get_nomatch = function(row){
  nomatch = as.numeric(as.character(row[[4]])) #all patients #
  high = as.numeric(as.character(row[[5]])) #of pats high match
  low = as.numeric(as.character(row[[6]]))
  nomatch = nomatch-high-low
}
sig_diff$no_match = apply(sig_diff, 1, get_nomatch)

match_h = sig_diff[,c(1:5)]
colnames(match_h)[5] = "num_cna_exp_match"
match_h$type = "HighExpDup"
match_h$num_cna_exp_match = as.numeric(as.character(match_h$num_cna_exp_match))

match_l = sig_diff[,c(1:4, 6)]
colnames(match_l)[5] = "num_cna_exp_match"
match_l$type = "LowExpDel"
match_l$num_cna_exp_match = as.numeric(as.character(match_l$num_cna_exp_match))

match_no = sig_diff[,c(1:4, 8)]
colnames(match_no)[5] = "num_cna_exp_match"
match_no$type = "NoExpCNAMatch"

matched_sig = rbind(match_h, match_l, match_no)
matched_sig = as.data.table(matched_sig)
matched_sig = matched_sig[order(num_cna_exp_match, -type)]
matched_sig = as.data.frame(matched_sig)
order = as.character(unique(matched_sig$gene))
matched_sig$gene <- factor(matched_sig$gene, levels = x$name[order(x$val)])
matched_sig$gene  # notice the changed order of factor levels

# Stacked bar plots, add labels inside bars
pdf("num_matching_CNA_exp_tags_23cancers_57lncs.pdf", width= 10)
p = ggbarplot(matched_sig, x = "name", y = "num_cna_exp_match",
  fill = "type", color = "type", 
  palette = "jco",
  label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.5)
ggpar(p,
 font.tickslab = c(5,"plain", "black"),
 xtickslab.rt = 45)

dev.off()

write.table(lnc_cna_cancer_data2, file="wilcoxon_CNA_expression_Results_all_cancers_TCGA_cands_May10.txt", quote=F, sep="\t", row.names=F)
#how to make sense of this? 
#plot segment? and where within it lies lncRNA? 

    #yplot = yplot + rremove("legend")
    #sp = sp + rremove("legend")
    #p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, nrow =2,align = "hv", 
    #      rel_widths = c(2, 1), rel_heights = c(1, 2))
   
    #plots <- align_plots(xplot, sp, align = 'v', axis = 'l')
    #bottom_row <- plot_grid(plots[[2]], yplot, labels = c('B', 'C'), align = 'h', rel_widths = c(2.85, 1))
    #p = plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.5))
    #print(p)