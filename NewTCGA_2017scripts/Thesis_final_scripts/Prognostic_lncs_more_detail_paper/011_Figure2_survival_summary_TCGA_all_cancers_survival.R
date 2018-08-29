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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)
library(patchwork)
library(tidyr)

rna = as.data.frame(rna)
dim(rna)
dim(pcg)
dim(norm)
dim(met)

all_genes = as.data.frame(unique(c(colnames(rna), colnames(pcg))))
z = which(str_detect(all_genes[,1], "ENSG"))
all_genes = all_genes[z,]
all_genes = as.data.frame(all_genes)
colnames(all_genes)[1] = "gene"

all_genes$type = ""
z = which(all_genes$gene %in% colnames(rna))
all_genes$type[z] = "lncRNA"
z = which(all_genes$gene %in% colnames(pcg))
all_genes$type[z] = "pcg"

saveRDS(all_genes, file="all_genes_used_in_TCGA_april17.rds")

###---------------------------------------------------------------

#function that tests each lncRNA's survival 

#1. remove discrepancy 
#z = which(rna$vital_status == "[Discrepancy]")
#rna = rna[-z,]

#2. list of cancers to apply function to 
cancers = as.list(unique(rna$Cancer))

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = cancers[which(!(cancers %in% canc_rm))]

#3. function that splits data into cancers 
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. function that calculates survival for each gene 
#det_lncs = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
#det_lncs =filter(det_lncs, status =="detectable")

canc_survival_genes = function(dato){
	#look at all lncRNAs that are expressed in at least some patients 
  z = which(str_detect(colnames(dato), "ENSG"))
  sums = apply(dato[,z], 2, sum)
  rm = names(sums[which(sums == 0)])
  z = which(colnames(dato) %in% rm)
  dato = dato[,-z]

  print(dato$type[1])

  z = which(str_detect(colnames(dato), "ENSG"))
  genes = unique(colnames(dato)[z])	

  #TEST------------------------------------------------------------------------------------------
  #genes = genes[1:100]
	canc_data_genes_analyze = dato 
	
	get_survival = function(gene){
	  print(gene)
  	results_cox <- as.data.frame(matrix(ncol=8)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients")
  	z = which(colnames(canc_data_genes_analyze) == gene)
  	dat = canc_data_genes_analyze[,c(1,z,(ncol(canc_data_genes_analyze)-33):ncol(canc_data_genes_analyze))]
  	dat$OS.time = as.numeric(dat$OS.time)
  	dat$OS = as.numeric(dat$OS)
  	#remove NAs
  	z = which(is.na(dat$OS.time))
  	if(!(length(z) ==0)){
  	dat = dat[-z,]}
	  med_gene = median(as.numeric(dat[,2]))  	
	  dat$med = ""
	  if(med_gene ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat[,2] > 0)
        l2 = which(dat[,2] ==0)
        dat$med[l1] = 1
        dat$med[l2] = 0
        }

    if(!(med_gene ==0)){
        l1 = which(dat[,2] >= med_gene)
        l2 = which(dat[,2] < med_gene)
        dat$med[l1] = 1
        dat$med[l2] = 0
    }

    check1 = table(dat$med)[1] >= 10
    check2 = table(dat$med)[2] >= 10
    if(check1 & check2){
      if(dim(table(dat$med)) ==2){
  	  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = dat)
    	hr = summary(res.cox)$coefficients[1,c(2)]
      num_pat = nrow(dat)
      if(hr > 1){
        risk = length(which(dat$med ==1))
      }
      if(hr <1){
        risk = length(which(dat$med ==0))
      }

      row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)], risk, num_pat)
     	names(row) <- names(results_cox)
    	return(row)
  	}}} #end get_survival function

	genes_survival = llply(genes, get_survival, .progress="text")
	genes_survival_res = ldply(genes_survival, rbind)
	#fdr
	colnames(genes_survival_res) = c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients")
	genes_survival_res$fdr = p.adjust(as.numeric(genes_survival_res$pval), method="fdr")
	genes_survival_res$canc = dato$Cancer[1]
	genes_survival_res = as.data.table(genes_survival_res)
	genes_survival_res = genes_survival_res[order(fdr)]
	return(genes_survival_res)
}

#DO NOT RUN 

#all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
#all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)

#saveRDS(all_cancers_genes_surv_comb, file="lncRNAs_for_plotting_HAzard_Ratios_Pvalues_July9.rds")


##############RUN-----------------------------------------------------------------------------------

all_cancers_genes_surv_comb = readRDS("lncRNAs_for_plotting_HAzard_Ratios_Pvalues_July9.rds")
canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]
colnames(canc_conv)[2] = "canc"
all_cancers_genes_surv_comb = merge(all_cancers_genes_surv_comb, canc_conv, by="canc")
write.csv(all_cancers_genes_surv_comb, file="ALL_lncRNAs_survival_august8.csv", quote=F, row.names=F)

###-------------------------------------------------------------------------------------------------

#plot scatter plot - HR versus p-value draw line for FDR = 0.05
all_cancers_genes_surv_comb$pval = -log10(as.numeric(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb$fdr = -log10(all_cancers_genes_surv_comb$fdr)
all_cancers_genes_surv_comb$HR = as.numeric(all_cancers_genes_surv_comb$HR)

z = which(is.na(all_cancers_genes_surv_comb$pval))
if(!(length(z)==0)){all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]}

z1 = which(all_cancers_genes_surv_comb$fdr == "Inf")
z2 = which(all_cancers_genes_surv_comb$upper95 == "Inf")
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-c(z1,z2),]

z = which(all_cancers_genes_surv_comb$HR > 10)
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

lineval = -log10(0.05)

all_cancers_genes_surv_comb$fdrsig = ""
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr < lineval] = "FDRnotSig"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval > lineval] = "Significant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval < lineval] = "NotSignificant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr >= lineval] = "FDRsig"
#z = which((all_cancers_genes_surv_comb$fdr < lineval) & (all_cancers_genes_surv_comb$fdr > -log10(0.1)))
#all_cancers_genes_surv_comb$fdrsig[z] = "FDR00.1"

#facet by cancer type 

#check overlap of sig prognostic lncRNAs bewteen cancer types
all_cancers_genes_surv_comb$risk_perc = as.numeric(all_cancers_genes_surv_comb$risk_size)/as.numeric(all_cancers_genes_surv_comb$num_patients)
all_cancers_genes_surv_comb$risk_perc_tag[(all_cancers_genes_surv_comb$risk_perc > 0.48) | (all_cancers_genes_surv_comb$risk_perc < 0.52)] = "75%_more_risk_group"
all_cancers_genes_surv_comb$risk_perc_tag[all_cancers_genes_surv_comb$risk_perc > 0.75] = "75%_more_risk_group"


#figure 2B/C?
sig_lncs = as.data.table(all_cancers_genes_surv_comb)
sig_lncs = as.data.table(subset(sig_lncs, fdrsig == "FDRsig"))

#sum freq
sig_lncs = as.data.table(table(sig_lncs$canc, sig_lncs$gene))
sig_lncs = as.data.table(filter(sig_lncs, N > 0))
colnames(sig_lncs)[1] = "canc"
sig_lncs = merge(sig_lncs, canc_conv, by="canc")

#calculate how many lncRNAs overlap between cancer types 
cancs = unique(sig_lncs$type)
#note these are FDR significant 

get_pairs = function(canc){
  canc_lncs = sig_lncs$V2[sig_lncs$type == canc]
  z = which(sig_lncs$V2 %in% canc_lncs)
  canc_pairs = sig_lncs[z,]
  cp = as.data.table(table(canc_pairs$type))
  cp$canc1 = canc
  colnames(cp)[1] = "canc2"
  return(cp)
}

all_canc_pairs = llply(cancs, get_pairs)
all_canc_pairs = ldply(all_canc_pairs)
all_canc_pairs = as.data.table(all_canc_pairs)
all_canc_pairs = all_canc_pairs[order(-N)]
all_canc_pairs$canc1 = factor(all_canc_pairs$canc1, levels=unique(all_canc_pairs$canc1))
all_canc_pairs$canc2 = factor(all_canc_pairs$canc2, levels=unique(all_canc_pairs$canc1))

pdf("overlap_FDR_sig_lncRNA_cands_bw_cancers_aug28.pdf", width=8, height=5)
g = ggplot(all_canc_pairs, aes(canc1, canc2)) +
  geom_tile(aes(fill=N)) +
  geom_text(aes(label = N), size=5) +
  scale_fill_gradient(low = "azure2", high = "orange", na.value = 'transparent') +
    xlab("Cancer 1") + ylab("Cancer 2") + theme_bw()
ggpar(g,
 font.tickslab = c(12,"plain", "black"),
 xtickslab.rt = 45, legend.title="# lncRNAs \noverlap")

dev.off()

#order by most significant to least significant 
order = as.data.table(dplyr::filter(as.data.table(table(all_cancers_genes_surv_comb$canc, all_cancers_genes_surv_comb$fdrsig)), V2 %in% c("FDRsig", "Significant")))
order = order[order(-V2,N)]

z = order[order(order$V1, -order$N),]
# Remove duplicates
z1 = z[!duplicated(z$V1),]
#order again
z1 = z1[order(-V2,N)]
order = z1$V1

all_cancers_genes_surv_comb$canc <- factor(all_cancers_genes_surv_comb$canc, levels = order)
all_cancers_genes_surv_comb$canc  # notice the changed order of factor levels

order_cols = c("NotSignificant", "Significant", "FDRsig")
all_cancers_genes_surv_comb$fdrsig <- factor(all_cancers_genes_surv_comb$fdrsig, levels = order_cols)
all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR > 1] = "Unfavourable"
all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR < 1] = "Favourable"

#Variation 1 of survival overview plot

pdf("HR_vs_pval_survival_all_cancers_scatter_plot_July9.pdf", width=12, height=9)
g = ggscatter(all_cancers_genes_surv_comb, x = "type", y = "HR", color="fdrsig", palett=c("gray34", mypal[1], "lightskyblue3"), size = 0.85) + 
geom_hline(yintercept=1, linetype="dashed", color = "red")
ggpar(g,
 font.xtickslab = c(8,"plain", "black"),
 xtickslab.rt = 90)
dev.off()


#Variation 2 of survival overview plot
head(all_cancers_genes_surv_comb)
all_cancers_genes_surv_comb$HR = log2(all_cancers_genes_surv_comb$HR)

# Change violin plot colors by groups
pdf("HR_vs_pval_survival_all_cancers_scatter_plot_july9_V2.pdf", width=10, height=8)
g = ggplot(all_cancers_genes_surv_comb, aes(type, HR)) +
  geom_violin() + 
  geom_jitter(height = 0.005, width = 0.005, aes(colour = factor(fdrsig)), size=0.15, alpha=0.5) +
  scale_colour_manual(name="colour", values=c("pink", "purple", "orange"))+ geom_hline(yintercept=log2(1), linetype="dashed", color = "red")
ggpar(g,
 font.xtickslab = c(8,"plain", "black"), ylab="log2(HR)",
 xtickslab.rt = 90)
dev.off()

#summarize number favourable and unfabourable lcnRNAs by fdr significance per cancer type

#get order of cancer types by total number of lncRNAs 
order = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig))
order = as.data.table(dplyr::filter(order, N >0))
order = as.data.table(dplyr::filter(order, V2 %in% c("Significant", "FDRsig")))
order = order[order(V2, -N)]
order = unique(order$V1)

summ = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig,
  all_cancers_genes_surv_comb$risk))
colnames(summ) = c("Cancer", "Sig", "Risk", "N")
summ = as.data.table(dplyr::filter(summ, N > 0))

#barplot----summary
#only include significant ones
#how many significant favourable vs unfavourable 
summ = as.data.table(dplyr::filter(summ, Sig %in% c("Significant", "FDRsig")))

#just sig
summ = dplyr::filter(summ, N > 0)
summ$Cancer = factor(summ$Cancer, levels = order)

#pdf("Univariate_summary_28_Cancers_july9.pdf", width=10)
part1 <- ggplot(data=summ, aes(x=Cancer, y=N, fill=Risk)) +
geom_bar(stat="identity")+
  theme_bw() + coord_flip() +
  scale_fill_manual(values=c('darkcyan','orange')) + ggtitle("Number of Univariate Significant lncRNAs, CoxPH p-val < 0.05")

part1 = ggpar(part1, legend="none",
 font.xtickslab = c(8,"plain", "black"), ylab="Number of lncRNAs")
#dev.off()

#just fdrsig
summ = as.data.table(dplyr::filter(summ, Sig == "FDRsig"))
summ = dplyr::filter(summ, N > 0)

#pdf("Univariate_summary_28_Cancers_july9_justfdr.pdf", width=10)
part2 <- ggplot(data=summ, aes(x=Cancer, y=N, fill=Risk)) +
geom_bar(stat="identity")+
  theme_bw() + coord_flip() +
  scale_fill_manual(values=c('darkcyan','orange'))

part2 = ggpar(part2, legend="bottom",
 font.xtickslab = c(8,"plain", "black"), ylab="Number of lncRNAs") + ggtitle("Number of Univariate Significant lncRNAs, CoxPH FDR < 0.05")
#dev.off()

#Figure 2A plot 
pdf("figure2_A_july10.pdf", width=7, height=7)
part1 + part2 + plot_layout(ncol = 1, heights = c(3, 2))
dev.off()

aggregate(summ[, 4], list(summ$Risk), sum)


### Get corrplot of survival candidates within each cancer type for significant lncRNAs ---------------------------------------------------
head(all_cancers_genes_surv_comb)
all_cancers_genes_surv_comb = as.data.table(all_cancers_genes_surv_comb)

cancers = (unique(all_cancers_genes_surv_comb$canc))

#function apply to each cancer type and plot corrplot for significant lncRNAs (sig & fdrsig)

library(corrplot)
library(Hmisc)


get_corplot = function(cancer){
  print(cancer)

  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))
  canc_genes = as.data.table(dplyr::filter(canc_genes, fdrsig %in% c("Significant", "FDRsig")))

  if((dim(canc_genes)[1] >=2)){

    if(dim(canc_genes)[1] >=50){
      canc_genes = canc_genes[order(-fdr)]
      canc_genes = canc_genes[1:50,]
    }

    genes = unique(canc_genes$gene)
    canc_exp = subset(rna, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    res2 <- rcorr(as.matrix(canc_exp), type="spearman")

    # Insignificant correlations are leaved blank
    col<- colorRampPalette(c("blue", "white", "red"))(20)
    c = corrplot(res2$r, order="hclust", col=col, 
      p.mat = res2$P, sig.level = 0.05, insig = "blank", tl.cex=1, method="color", mar=c(0,0,1,0),bg="snow2",
      tl.pos = "n", title=paste(cancer, "top 50 lncRNA correlations"))
    print(c)
  }
}

pdf("correlation_plots.pdf")
llply(cancers, get_corplot, .progress = "text")
dev.off()

#saveRDS(rna, file="lncRNA_expression_for_remote_plotting_july10.rds")
#saveRDS(pcg, file="pcg_expression_for_remote_plotting_july10.rds")

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

get_corplot = function(cancer){
  print(cancer)

  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))
  canc_genes = as.data.table(dplyr::filter(canc_genes, fdrsig %in% c("Significant", "FDRsig")))

  if((dim(canc_genes)[1] >=2)){

    if(dim(canc_genes)[1] >=50){
      canc_genes = canc_genes[order(-fdr)]
      canc_genes = canc_genes[1:50,]
    }

    genes = unique(canc_genes$gene)
    canc_exp = subset(rna, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    
    res2 = rcorr(as.matrix(canc_exp), type="spearman")
    res2 = flattenCorrMatrix(res2$r, res2$P)
    res2$fdr = p.adjust(res2$p, method="fdr")
    res2 = as.data.table(res2)
    res2 = res2[order(fdr)]
    res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05))
  
    #check if lncRNA-lncRNA correlations match HRs 
    check_dir = function(lnc1, lnc2){
      hr_lnc1 = canc_genes$HR[canc_genes$gene == lnc1]
      hr_lnc2 = canc_genes$HR[canc_genes$gene == lnc2]

      check1 = ((hr_lnc1 > 0) & (hr_lnc2 > 0)) 
      check2 = ((hr_lnc1 < 0) & (hr_lnc2 < 0)) 
      
      if(check1){match = "U"
      }else if(check2){
        match = "F"
      }else{match = "D"}

      return(match)
    }

    res2$match = mapply(check_dir, res2$row, res2$column)

    #ordered by strongest correlations to weakest correlations
    res2 = res2[order(match, cor)]
    res2$row = factor(res2$row, levels = unique(res2$row))
    res2$column = factor(res2$column, levels = unique(res2$column))

    g = ggplot(res2, aes(row, column)) +
        geom_tile(aes(fill = cor)) +
        geom_text(aes(label = match), size=2) +
        scale_fill_gradient2(low = "blue", high = "red", mid="white", midpoint = 0, na.value = 'transparent') +
        xlab("lncRNA1") + ylab("lncRNA2") + theme_bw() + coord_fixed() + 
        ggtitle(paste(cancer, "top 50 lncRNA correlations"))       

    g = ggpar(g,
        font.tickslab = c(4,"plain", "black"),
        xtickslab.rt = 45, legend.title="Spearman \nCorrelation")

    print(g)

  }
}

pdf("correlation_plots_with_directions.pdf")
llply(cancers, get_corplot, .progress = "text")
dev.off()

#Calculate how many lncRNA-lncRNA pairs there are per each cancer type 
#and how many of each are significntl correalted 

get_summary = function(cancer){
  print(cancer)

  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))
  canc_genes = as.data.table(dplyr::filter(canc_genes, fdrsig %in% c("Significant", "FDRsig")))

  if((dim(canc_genes)[1] >=2)){

    genes = unique(canc_genes$gene)
    canc_exp = subset(rna, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    
    res2 = rcorr(as.matrix(canc_exp), type="spearman")
    res2 = flattenCorrMatrix(res2$r, res2$P)
    res2$fdr = p.adjust(res2$p, method="fdr")
    res2 = as.data.table(res2)
    res2 = res2[order(fdr)]

    #total pairs 
    tot_pairs = nrow(res2)
    res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05))
    sig_pairs = nrow(res2)    

    #%
    perc = sig_pairs/tot_pairs

    row = c(as.character(cancer), tot_pairs, sig_pairs, perc)
    return(row)
  }
}

canc_results = llply(cancers, get_summary, .progress = "text")
#remove null
canc_results = Filter(Negate(is.null), canc_results)

canc_results = do.call(rbind.data.frame, canc_results)
colnames(canc_results) = c("cancer", "total_pairs", "sig_pairs", "perc")

#Calculate for all significant co-exprseed pairs 
#how many matching correlated unfavourable/favourable
#how many unfavourable or not matching

get_pairs_results = function(cancer){
  print(cancer)

  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))
  canc_genes = as.data.table(dplyr::filter(canc_genes, fdrsig %in% c("Significant", "FDRsig")))

  if((dim(canc_genes)[1] >=2)){

    genes = unique(canc_genes$gene)
    canc_exp = subset(rna, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    
    res2 = rcorr(as.matrix(canc_exp), type="spearman")
    res2 = flattenCorrMatrix(res2$r, res2$P)
    res2$fdr = p.adjust(res2$p, method="fdr")
    res2 = as.data.table(res2)
    res2 = res2[order(fdr)]
    tot_pairs = nrow(res2)
    res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05))
    sig_pairs = nrow(res2)  
    #check if lncRNA-lncRNA correlations match HRs 
    check_dir = function(lnc1, lnc2){
      hr_lnc1 = canc_genes$HR[canc_genes$gene == lnc1]
      hr_lnc2 = canc_genes$HR[canc_genes$gene == lnc2]

      check1 = ((hr_lnc1 > 0) & (hr_lnc2 > 0)) 
      check2 = ((hr_lnc1 < 0) & (hr_lnc2 < 0)) 
      
      if(check1){match = "U"
      }else if(check2){
        match = "F"
      }else{match = "D"}

      return(match)
    }

    res2$match = mapply(check_dir, res2$row, res2$column)

    #ordered by strongest correlations to weakest correlations
    res2 = res2[order(match, cor)]
    res2$cor_sum[res2$cor > 0] = "Pos"
    res2$cor_sum[res2$cor < 0] = "Neg"
    #summarize how many of each kind
    t = table(res2$match, res2$cor_sum)
    t = as.data.table(tidy(t))
    t = t[order(Freq)]
    t$total_sig_pairs = sig_pairs
    t$total_pairs = tot_pairs
    t$perc = t$Freq/sig_pairs
    t$cancer = cancer
    return(t)
  }
}

canc_results_pairs_types = llply(cancers, get_pairs_results, .progress = "text")

#save 
saveRDS(canc_results_pairs_types, file="correlation_lnc_lnc_results_july14.rds")

#remove null
canc_results_pairs_types2 = Filter(Negate(is.null), canc_results_pairs_types)
canc_results_pairs_types2 = ldply(canc_results_pairs_types2)
colnames(canc_results_pairs_types2)[1:3] = c("HR_pair", "Exp_pair", "N")
canc_results_pairs_types2 = as.data.table(canc_results_pairs_types2)
canc_results_pairs_types2 = as.data.table(canc_results_pairs_types2[order(perc)])
canc_results_pairs_types2$perc = round(as.numeric(canc_results_pairs_types2$perc), digits=4)
colnames(canc_conv)[2] = "cancer"
canc_results_pairs_types2 = merge(canc_results_pairs_types2, canc_conv, by="cancer")

#Summarize figure 2 B - summary of correlated pairs across 
#cancer types 
canc_results_pairs_types2$HR_pair = as.character(canc_results_pairs_types2$HR_pair)
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$HR_pair == "F"] = "Both \nFavourable"
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$HR_pair == "U"] = "Both \nUnfavourable"
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$HR_pair == "D"] = "Opposite \nHRs"

#cancer order keep same as first plot
canc_results_pairs_types2$type <- factor(canc_results_pairs_types2$type, levels = rev(order))
canc_results_pairs_types2$column_name = paste(canc_results_pairs_types2$HR_pair, canc_results_pairs_types2$Exp_pair)

#find matched correlation - hazard ratios pairs
canc_results_pairs_types2$match_type = ""

for(i in 1:nrow(canc_results_pairs_types2)){
  hr_pair = canc_results_pairs_types2$HR_pair[i]
  co = as.character(canc_results_pairs_types2$Exp_pair[i])
    
  if((hr_pair == "Both \nFavourable") & (co == "Pos")){
    match = "Yes"
  } else if ((hr_pair == "Both \nUnfavourable") & (co == "Pos")){
    match = "Yes"
  } else if ((hr_pair == "Opposite \nHRs") & (co == "Neg")){
    match = "Yes"
  } else {match = "No"}

  canc_results_pairs_types2$match_type[i] = match
}

saveRDS(canc_results_pairs_types2, file="correlation_lnc_lnc_results_july14.rds")

#canc_results_pairs_types2 = as.data.table(dplyr::filter(canc_results_pairs_types2, match_type == "Yes"))

pdf("Figure2B_summary_types_of_correlations_28_cancers_types.pdf", height=5, width=6)

g = ggplot(canc_results_pairs_types2, aes(type, column_name)) +
  geom_tile(aes(fill = perc), color="grey")  +
    scale_fill_gradient2(low = "turquoise4", high = "tan1") +
    xlab("Cancer") + ylab("lncRNA pairs Hazard Ratios") + theme_bw() +
     labs(fill="% of Significant \nCorrelated Pairs", colour="Type of \ncorrelation")
ggpar(g,
 font.tickslab = c(6,"plain", "black"),
 xtickslab.rt = 45)

dev.off()



















