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
library(mixtools)

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

#cindices
r = readRDS("lncs_cross_validations_Results_nov20.rds")

#Co-expression results of PCGs 
#coexp = readRDS("coexpression_results_processed_july24.rds")
#coexp$combo2 = paste(coexp$lnc, coexp$canc, coexp$pcg, sep="_")

#-most recent differential expression analysis results from running LIMMA
#using counts and VOOM 

all_de_results = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
all_de_results = as.data.table(all_de_results)

#PCG lncRNA results
pcg_lnc = readRDS("summary_pcg_analysis_wHRs_jul2y24.rds") #all these have at least 1, 50-pcg signature 
pcg_lnc = pcg_lnc[order(-NumPCGs)]
pcg_lnc$HR = as.numeric(pcg_lnc$HR)
pcg_lnc$lnc_stat = ""
pcg_lnc$lnc_stat[which(pcg_lnc$HR < 0)] = "Favourable"
pcg_lnc$lnc_stat[which(pcg_lnc$HR > 0)] = "Unfavourable"

mypal5 = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#------FUNCTIONS----------------------------------------------------

#######
##[1]##-------------------------------------------------------------
#######

#get gene ID from name 
#input: GeneName
#output: GeneID 

get_name_pcg = function(pcg){
	z = which(ucsc$hg19.ensGene.name2 == pcg)
  if(length(z)>1){
    z = z[1]
  }
	return(ucsc$hg19.ensemblToGeneName.value[z])
}

get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}

get_ensg = function(lnc){
  z = which(fantom$CAT_geneName == lnc)
  return(fantom$CAT_geneID[z])
}

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

#######
##[2]##-------------------------------------------------------------
#######

#get lncRNA expression across cancers plot 
#input: lncRNA id
#output: boxplot of expression across cancer types

#test case: pca3 - ENSG00000225937
#mypal = readRDS("best_pal.rds")

get_exp_plots = function(lnc){
	dat = rna[,which(colnames(rna) %in% c("type", lnc))]
	colnames(dat)[1] = "lncRNAExp"
	dat$lncRNAExp = log1p(dat$lncRN)
	dat = as.data.table(dat)
	dat = dat[order(lncRNAExp)]

	#get order of cancer types by median 
	sum = as.data.table(dat %>% dplyr::group_by(type) %>% 
		dplyr::summarize(median=quantile(lncRNAExp, 0.75)))
	sum = sum[order(median)]

	dat$type = factor(dat$type, levels=unique(sum$type))
	#p = ggboxplot(dat, x="type", y = "lncRNAExp", title=lnc, error.plot = "linerange")+
	#rotate_x_text(45)
	
	#try density plot
	p = ggdensity(dat, x = "lncRNAExp", title=paste(lnc, get_name(lnc)), color="type", palette=mypal5)+
  xlab("log1p(FPKM)")
	print(p)
	print("done plot")

}


get_num_peaks_den = function(lnc){
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
  
  #how many peaks within each cancer?
  types = as.character(unique(dat$type))
  get_peaks = function(canc){
     D = dat$lncRNAExp[dat$type == canc]
     
     firstquantile = summary(D)[2]
     thirdquantile = summary(D)[5]
     med = median(D)

     #some data
     d <- density(D)

     #make it a time series
     ts_y<-ts(d$y)

     #calculate turning points (extrema)
     require(pastecs)
     tp<-turnpoints(ts_y)
     print(canc)
     row = c(lnc, canc, firstquantile, med, thirdquantile, length(d$x[tp$tppos]))
     names(row) = c("lnc", "cancer", "firstQ", "median", "thirdQ", "numPeaks")
     return(row)
  }
  canc_peaks = llply(types, get_peaks)
  canc_peaks = ldply(canc_peaks)
  canc_peaks = as.data.table(canc_peaks)
  canc_peaks$firstQ = as.numeric(canc_peaks$firstQ)
  canc_peaks$median = as.numeric(canc_peaks$median)
  canc_peaks$thirdQ = as.numeric(canc_peaks$thirdQ)
  canc_peaks = canc_peaks[order(-firstQ, -median, -thirdQ)]
  canc_peaks$med_count[canc_peaks$median == 0] = 0
  canc_peaks$med_count[canc_peaks$median > 0] = 1
  canc_peaks$firstQ_count[canc_peaks$firstQ == 0] = 0
  canc_peaks$firstQ_count[canc_peaks$firstQ > 0] = 1
  canc_peaks$thirdQ_count[canc_peaks$thirdQ == 0] = 0
  canc_peaks$thirdQ_count[canc_peaks$thirdQ > 0] = 1
  s = as.data.table(table(canc_peaks$firstQ_count, canc_peaks$med_count, canc_peaks$thirdQ_count))
  s = as.data.table(filter(s, N >0))
  if(dim(s)[1] > 1){
    lnc_type = "bimodal"
  }
   if(dim(s)[1] == 1){
    lnc_type = "unimodal"
  }
  return(c(lnc, lnc_type))
}

pdf("pca3_expression_across_cancers.pdf")
get_exp_plots("ENSG00000225937")
get_num_peaks_den("ENSG00000225937")
dev.off()


#pdf("all_lncs_cands_dist_plots.pdf")
#lncs = unique(allCands$gene)
#llply(lncs, get_exp_plots, .progress="text")
#dev.off()

#lnc_types = llply(lncs, get_num_peaks_den, .progress="text")


#######
##[x]##-------------------------------------------------------------
#######

#identify which lncRNAs are bimodal and which ones are 
#highly expressed in all cancers 

#get peaks 
#density(data$V2)$x[which.max(density(data$V2)$y)]



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
   g = ggboxplot(dat, x = "lnc_median", y="pcgExp", title=paste(lnc, pcg, mean_diff), fill="lnc_median") +
   	stat_compare_means()
   print(g)
}

#######
##[5]##-------------------------------------------------------------
#######

gene = "ENSG00000225937" #PCA3
cancer = "PAAD"

###EASY WAY TO MAKE KM PLOT
get_km_plot = function(gene, cancer){
  all_g = all
  all_g = as.data.frame(all_g)
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
  gene_name = get_name(gene)
  if(is.na(gene_name)){
    gene_name = get_name_pcg(gene)
  }
  cox_mod = coxph(Surv(OS.time, OS) ~ gene, data = dat)
  print(glance(cox_mod)$concordance)
  conc = round(glance(cox_mod)$concordance, digits=2)
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene_name, dat$type[1], "\nConcordance=", conc),
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


#model with how many features would have the best performance? 

#######
##[8]##-------------------------------------------------------------
#######

###ROC SURVIVAL CURVES 
#https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/


#Compare expression of lncRNAs to mRNAs 
rownames(rna) = rna$patient
#keep cancer type 
rna = rna[,c(2:5786)]
#rna = as.data.frame(rna)
rna = t(rna)
rna = as.data.frame(rna)

#get median for each gene 
rna$median = apply(rna[,1:8503], 1, median)
rna$CAT_geneID = rownames(rna)

rnaa = merge(rna, fantom, by="CAT_geneID")



#Compare expression of lncRNAs to mRNAs 
rownames(pcg) = pcg$patient
#keep cancer type 
pcg = pcg[,c(2:19351)]
#rna = as.data.frame(rna)
pcg = t(pcg)
pcg = as.data.frame(pcg)

#get median for each gene 
pcg$median = apply(pcg, 1, median)
pcg$hg19.ensGene.name2 = rownames(pcg)
pcga = merge(pcg, ucsc, by="hg19.ensGene.name2")

#what i need 
#gene type and median and gene name 

rnaa = rnaa[,c("CAT_geneID", "median", "CAT_geneClass")]
pcga = pcga[,c("hg19.ensGene.name2", "median", "hg19.ensemblSource.source")]

colnames(rnaa) = c("gene", "median", "type")
rnaa$typesec = "lncRNA"
colnames(pcga) = c("gene", "median", "type")
pcga = subset(pcga, type == "protein_coding")
pcga$typesec = "pcg" #17196 genes 


all = rbind(pcga, rnaa) #22981 genes in total
all$medlog = log1p(all$median)

all$type = factor(all$type, levels = c("lncRNA_intergenic", "lncRNA_sense_intronic", "lncRNA_antisense", 
  "lncRNA_divergent", "protein_coding"))

all = as.data.table(all)
all = all[order(medlog)]
all$rank = 1:nrow(all)
all$score = all$rank/nrow(all)

pdf("all_genes_plot_median_exp_all_cancers.pdf")

#boxplot
p = ggplot(all, aes(type, medlog))
p + geom_boxplot(outlier.alpha = 0.1) + stat_n_text()


dev.off()















