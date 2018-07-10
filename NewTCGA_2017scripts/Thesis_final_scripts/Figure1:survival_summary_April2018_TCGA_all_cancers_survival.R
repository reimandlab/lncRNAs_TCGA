###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

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
z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]

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

#all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
#all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)

#saveRDS(all_cancers_genes_surv_comb, file="lncRNAs_for_plotting_HAzard_Ratios_Pvalues_July9.rds")


##############RUN-----------------------------------------------------------------------------------

all_cancers_genes_surv_comb = readRDS("lncRNAs_for_plotting_HAzard_Ratios_Pvalues_July9.rds")
canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]
colnames(canc_conv)[2] = "canc"
all_cancers_genes_surv_comb = merge(all_cancers_genes_surv_comb, canc_conv, by="canc")

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

#order by most significant to least significant 
order = as.data.table(filter(as.data.table(table(all_cancers_genes_surv_comb$canc, all_cancers_genes_surv_comb$fdrsig)), V2 %in% c("FDRsig", "Significant")))
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
pdf("HR_vs_pval_survival_all_cancers_scatter_plot_july9.pdf", width=10, height=8)
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
order = as.data.table(filter(order, N >0))
order = as.data.table(filter(order, V2 %in% c("Significant", "FDRsig")))
order = order[order(V2, -N)]
order = unique(order$V1)

summ = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig,
  all_cancers_genes_surv_comb$risk))
colnames(summ) = c("Cancer", "Sig", "Risk", "N")
summ = as.data.table(filter(summ, N > 0))

#barplot----summary
#only include significant ones
#how many significant favourable vs unfavourable 
summ = as.data.table(filter(summ, Sig %in% c("Significant", "FDRsig")))

#just sig
summ = filter(summ, N > 0)
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
summ = as.data.table(filter(summ, Sig == "FDRsig"))
summ = filter(summ, N > 0)

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

  #get sig genes
  canc_genes = as.data.table(filter(all_cancers_genes_surv_comb, canc == cancer))
  canc_genes = as.data.table(filter(canc_genes, fdrsig %in% c("Significant", "FDRsig")))

  if(!(dim(canc_genes)[1] == 0)){
    genes = unique(canc_genes$gene)
    canc_exp = subset(rna, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    res2 <- rcorr(as.matrix(canc_exp))

    # Insignificant correlations are leaved blank
    col<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(res2$r, order="hclust", col=col, 
      p.mat = res2$P, sig.level = 0.01, insig = "blank", tl.cex=1, method="color", 
      cl.pos = "n", tl.pos = "n")


  }

}











