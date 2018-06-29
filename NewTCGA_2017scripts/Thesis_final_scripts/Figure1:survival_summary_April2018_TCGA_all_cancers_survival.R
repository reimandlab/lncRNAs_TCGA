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

#3. function that splits data into cancers 
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. function that calculates survival for each gene 

det_lncs = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
det_lncs =filter(det_lncs, status =="detectable")

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
  	results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
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

    if(dim(table(dat$med)) ==2){
	  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = dat)
  	row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox)
  	return(row)
  	}}

	genes_survival = llply(genes, get_survival, .progress="text")
	genes_survival_res = ldply(genes_survival, rbind)
	#fdr
	colnames(genes_survival_res) = c("gene", "coef", "HR", "pval", "low95", "upper95")
	genes_survival_res$fdr = p.adjust(as.numeric(genes_survival_res$pval), method="fdr")
	genes_survival_res$canc = dato$Cancer[1]
	genes_survival_res = as.data.table(genes_survival_res)
	genes_survival_res = genes_survival_res[order(fdr)]
	return(genes_survival_res)
}

all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)


saveRDS(all_cancers_genes_surv_comb, file="lncRNAs_for_plotting_HAzard_Ratios_Pvalues_June28.rds")


all_cancers_genes_surv_comb = readRDS("lncRNAs_for_plotting_HAzard_Ratios_Pvalues_June28.rds")


###---------------------------------------------------------------

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


#Variation 1 of survival overview plot

pdf("HR_vs_pval_survival_all_cancers_scatter_plot_June28.pdf", width=12, height=9)
g = ggscatter(all_cancers_genes_surv_comb, x = "canc", y = "HR", color="fdrsig", palett=c("gray34", mypal[1], "lightskyblue3"), size = 0.85) + 
geom_hline(yintercept=1, linetype="dashed", color = "red")
ggpar(g,
 font.xtickslab = c(8,"plain", "black"),
 xtickslab.rt = 90)
dev.off()


#Variation 2 of survival overview plot

head(all_cancers_genes_surv_comb)

all_cancers_genes_surv_comb$HR = log2(all_cancers_genes_surv_comb$HR)

# Change violin plot colors by groups
pdf("HR_vs_pval_survival_all_cancers_scatter_plot_May23.pdf", width=10, height=8)

g = ggplot(all_cancers_genes_surv_comb, aes(canc, HR)) +
  geom_violin() + 
  geom_jitter(height = 0.005, width = 0.005, aes(colour = factor(fdrsig)), size=0.15, alpha=0.5) +
  scale_colour_manual(name="colour", values=c("pink", "purple", "orange"))+ geom_hline(yintercept=log2(1), linetype="dashed", color = "red")


ggpar(g,
 font.xtickslab = c(8,"plain", "black"), ylab="log2(HR)",
 xtickslab.rt = 90)

dev.off()



####summarize how many significant lncRNAs per cancer type

all_cancers_genes_surv_comb = as.data.table(all_cancers_genes_surv_comb)

sig = as.data.table(filter(all_cancers_genes_surv_comb, fdrsig %in% c("FDRsig", "Significant")))
sig_counts = as.data.table(table(sig$gene, sig$canc))
sig_counts = as.data.table(filter(sig_counts, N >0))
sig_counts = sig_counts[order(N)]
sig_counts = as.data.table(table(sig_counts$V1))
sig_counts = sig_counts[order(N)]



