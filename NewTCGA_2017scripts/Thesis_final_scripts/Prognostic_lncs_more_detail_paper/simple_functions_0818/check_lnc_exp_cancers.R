###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
rna = rna[-which(rna$Cancer %in% canc_rm),]

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer", "type", "OS", "OS.time"))
#all = all[,1:25170]

#lncRNA candidates, n = 166, n=173 combos 
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#Co-expression results of PCGs 
coexp = readRDS("coexpression_results_processed_july24.rds")
coexp$combo2 = paste(coexp$lnc, coexp$canc, coexp$pcg, sep="_")

#PCG lncRNA results
pcg_lnc = readRDS("summary_pcg_analysis_wHRs_jul2y24.rds") #all these have at least 1, 50-pcg signature 
pcg_lnc = pcg_lnc[order(-NumPCGs)]
pcg_lnc$HR = as.numeric(pcg_lnc$HR)
pcg_lnc$lnc_stat = ""
pcg_lnc$lnc_stat[which(pcg_lnc$HR < 0)] = "Favourable"
pcg_lnc$lnc_stat[which(pcg_lnc$HR > 0)] = "Unfavourable"


#------FUNCTIONS----------------------------------------------------

#######
##[1]##-------------------------------------------------------------
#######

#get gene ID from name 
#input: GeneName
#output: GeneID 

get_ensg = function(lnc){
	z = which(fantom$CAT_geneName == lnc)
	return(fantom$CAT_geneID[z])
}

get_name = function(ensg){
    z = which(ucsc$hg19.ensGene.name2 == ensg)
    return(ucsc$hg19.ensemblToGeneName.value[z][1])
}

#######
##[2]##-------------------------------------------------------------
#######

#get lncRNA expression across cancers plot 
#input: lncRNA id
#output: boxplot of expression across cancer types

#test case: pca3 - ENSG00000225937

get_exp_plots = function(lnc){
	dat = rna[,which(colnames(rna) %in% c("type", lnc))]
	colnames(dat)[1] = "lncRNAExp"
	dat$lncRNAExp = log1p(dat$lncRN)
	dat = as.data.table(dat)
	dat = dat[order(lncRNAExp)]

	#get order of cancer types by median 
	sum = as.data.table(dat %>% dplyr::group_by(type) %>% 
		dplyr::summarize(median=median(lncRNAExp)))
	sum = sum[order(median)]

	dat$type = factor(dat$type, levels=unique(sum$type))
	#p = ggboxplot(dat, x="type", y = "lncRNAExp", title=lnc, error.plot = "linerange")+
	#rotate_x_text(45)
	
	#try density
	p = ggdensity(dat, x = "lncRNAExp", title=paste(lnc, get_name(lnc)), color="type")+
  xlab("log1p(FPKM)")
	print(p)
	print("done plot")
}

pdf("pca3_expression_across_cancers.pdf")
get_exp_plots("ENSG00000225937")
dev.off()

lnc = get_ensg("MALAT1")
pdf("malat1_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

lnc = get_ensg("NEAT1")
pdf("neat1_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

lnc = get_ensg("HOTAIR")
pdf("hotair_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

lnc = get_ensg("TINCR")
pdf("tincr_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

lnc = get_ensg("XIST")
pdf("xist_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

lnc = get_ensg("RP5-1158E12.3")
pdf("RP5-1158E12.3_expression_across_cancers.pdf")
get_exp_plots(lnc)
dev.off()

#plot all lncRNAs 
lncs = unique(colnames(rna))
z = which(str_detect(lncs, "ENSG"))
lncs = lncs[z]

#pdf("all_lncs_dist_plots.pdf")
#llply(lncs, get_exp_plots, .progress="text")
#dev.off()

#######
##[3]##-------------------------------------------------------------
#######

#check how many patients have expression value X FPKM 
#(because we want at least 15 patients with 100FPKM expression)

get_num_pats = function(lnc, canc, exp_cut){
	dat = rna[,which(colnames(rna) %in% c("Cancer", lnc))]
	dat = dat[which(dat$Cancer == canc),]
	colnames(dat)[1] = "lncRNAExp"
	num_pats = length(which(dat[,1] > exp_cut))
	if(num_pats >= 15){
		return("great success")
	}
	if(num_pats < 15){
		return("problem")
	}
}

#######
##[4]##-------------------------------------------------------------
#######

#check expression of PCG how it's different between a lncRNA's risk group
#input: lncRNA, pcg, cancer 

get_pcg_enrich = function(lnc, pcg, canc){
	dat = all[,which(colnames(all) %in% c("Cancer", lnc, pcg, "patient"))]
	dat = dat[which(dat$Cancer == canc),]
	colnames(dat)[which(colnames(dat)==lnc)] = "lncRNAExp"
	colnames(dat)[which(colnames(dat)==pcg)] = "pcgExp"
	
	#add risk group
	lnc_hr = as.numeric(allCands$HR[which(allCands$combo == paste(lnc, canc, sep="_"))])

	#get med
	med = median(dat$lncRNAExp)
	dat$lnc_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat$lncRNAExp > 0)
        l2 = which(dat$lncRNAExp ==0)
        dat$lnc_median[l1] = 1
        dat$lnc_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(dat$lncRNAExp >= med)
        l2 = which(dat$lncRNAExp < med)
        dat$lnc_median[l1] = 1
        dat$lnc_median[l2] = 0
    }

   if(lnc_hr > 1){
   	dat$lnc_median[dat$lnc_median == 1] = "risk"
   	dat$lnc_median[dat$lnc_median == 0] = "less_risk"
   } 

    if(lnc_hr < 1){
   	dat$lnc_median[dat$lnc_median == 0] = "risk"
   	dat$lnc_median[dat$lnc_median == 1] = "less_risk"
   } 

   mean_diff = round(as.numeric(coexp$mean_diff[which(coexp$combo2 == paste(lnc, canc, pcg, sep="_"))]), digits=4)
   dat$pcgExp = log1p(dat$pcgExp)
   g = ggboxplot(dat, x = "lnc_median", y="pcgExp", title=paste(lnc, pcg, mean_diff)) +
   	stat_compare_means()
   print(g)
}

pdf("pcg_diff_exp_example.pdf")
get_pcg_enrich(get_ensg("MAPT-AS1"), get_ensg("MAPT"), "Breast invasive carcinoma"))
dev.off()


#######
##[5]##-------------------------------------------------------------
#######

gene = "ENSG00000225937" #PCA3
cancer = "PAAD"

###EASY WAY TO MAKE KM PLOT
get_km_plot = function(gene, cancer){
  dat = all[,c(which(colnames(all) %in% c("type", gene, "OS", "OS.time")))]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients 
  med = median(dat$gene)
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)
  dat$OS.time = dat$OS.time/365
  dat$gene = factor(dat$gene, levels = c(0,1))
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(get_name(gene), dat$type[1]),
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}
}	

pdf("pca3_km_plot_prad.pdf")
get_km_plot("ENSG00000225937", "PRAD")
dev.off()

pdf("hotair_km_plot_prad.pdf")
gene = get_ensg("HOTAIR")
get_km_plot(gene, "BRCA")
dev.off()

pdf("malat1_esca_km_plot_prad.pdf")
gene = get_ensg("MALAT1")
get_km_plot(gene, "ESCA")
dev.off()

pdf("malat1_stad_km_plot_prad.pdf")
gene = get_ensg("MALAT1")
get_km_plot(gene, "STAD")
dev.off()

pdf("malat1_ov_km_plot_prad.pdf")
gene = get_ensg("MALAT1")
get_km_plot(gene, "OV")
dev.off()

pdf("xist_brca_km_plot_prad.pdf")
gene = get_ensg("XIST")
get_km_plot(gene, "BRCA")
dev.off()


#######
##[7]##-------------------------------------------------------------
#######

###EAST WAY TO FOREST PLOT (1 variable)
get_forest_plot = function(gene, cancer){
  dat = all[,c(which(colnames(all) %in% c("type", gene, "OS", "OS.time")))]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients 
  med = median(dat$gene)
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ gene, data = dat)

  #compare coefficeints
  print(ggforest(cox_lnc, data=exp_dat))
  print("done")
}
}


#######
##[7]##-------------------------------------------------------------
#######

###EAST WAY TO BUILD AND TEST MULTIVARIATE MODELS

#*****
ggstatsplot::ggcoefstats(x = stats::lm(formula = mpg ~ am * cyl,
                                       data = datasets::mtcars)) 


ggstatsplot::ggcoefstats(
  x = stats::lm(formula = mpg ~ am * cyl,
                data = datasets::mtcars),
  point.color = "red",
  vline.color = "#CC79A7",
  vline.linetype = "dotdash",
  stats.label.size = 3.5,
  stats.label.color = c("#0072B2", "#D55E00", "darkgreen"),
  title = "Car performance predicted by transmission and cylinder count",
  subtitle = "Source: 1974 Motor Trend US magazine"
) +                                    
  # further modification with the ggplot2 commands
  # note the order in which the labels are entered
  ggplot2::scale_y_discrete(labels = c("transmission", "cylinders", "interaction")) +
  ggplot2::labs(x = "regression coefficient",
                y = NULL)

  


get_surv_model = function(){

}

#model with how many features would have the best performance? 

#######
##[8]##-------------------------------------------------------------
#######

###ROC SURVIVAL CURVES 
#https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/









