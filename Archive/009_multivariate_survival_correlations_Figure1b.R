###---------------------------------------------------------------
###009_multivariate_survival_correlations_Figure1b.R
###---------------------------------------------------------------

#what? 
#Systematic analysis of how prognostic lncRNAs are across cancer types
#After finding prognostic lncRNAs --> look at how co-expressed these pairs of lncRNAs are 
#This script also generates supplementary tables 1 and 2 

#load all required data
source("universal_LASSO_survival_script.R")

#load libraries 
require(caTools)
library(survAUC)
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

date = Sys.Date()

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

#summarize lncRNAs that we studied 
lncs = colnames(rna)[which(str_detect(colnames(rna), "ENSG"))]
z = which(fantom$CAT_geneID %in% lncs)
fantom = fantom[z,]
write.csv(fantom, file="5785_lncRNAs_used_in_study_table2.csv", quote=F, row.names=F)

#summarize patients 
pats = unique(rna[,c("type", "Cancer")])
tt = as.data.table(table(rna$type))
colnames(tt) = c("type", "num_patients")
tt = merge(tt, pats, by="type")
tt = tt[order(num_patients)]
write.csv(tt, file="TCGA_cancer_types_used_in_study_table1.csv", quote=F, row.names=F)

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

#1. list of cancers to apply function to 
cancers = as.list(unique(rna$Cancer))

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = cancers[which(!(cancers %in% canc_rm))]

#2. function that splits data into cancers 
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#3. function that calculates survival for each gene 

canc_survival_genes = function(dato){
	
  #look at all lncRNAs that are expressed in at least some patients 
  z = which(str_detect(colnames(dato), "ENSG"))
  sums = apply(dato[,z], 2, sum)
  rm = names(sums[which(sums == 0)])
  if(!(length(rm)==0)){
  z = which(colnames(dato) %in% rm)
  dato = dato[,-z]
  }

  print(dato$type[1])

  z = which(str_detect(colnames(dato), "ENSG"))
  genes = unique(colnames(dato)[z])	

  #TEST------------------------------------------------------------------------------------------
  #genes = genes[1:100]
	canc_data_genes_analyze = dato 
	
	get_survival = function(gene){
	  print(gene)
  	results_cox <- as.data.frame(matrix(ncol=8)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients")

    #keep gene expresion as well as clinical columns 
    cols_keep = c(gene, "patient", "age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "OS", "OS.time")
    dat = canc_data_genes_analyze[,cols_keep]
    dat$age_at_initial_pathologic_diagnosis = as.numeric(dat$age_at_initial_pathologic_diagnosis)
                
              #remove columns with less than 2 contrasts 
              check_contrasts = function(col){
                check = dim(table(col))
                if(check >1){
                  return("keep")
                }
              }

    keep = unlist(apply(dat, 2, check_contrasts))
    z = which(colnames(dat) %in% names(keep))
    dat = dat[,z]
              
    dat$OS = as.numeric(as.character(dat$OS))
    dat$OS.time = as.numeric(as.character(dat$OS.time))

    z = which(colnames(dat) %in% gene)
    colnames(dat)[z] = "gene"
    dat$gene = as.numeric(dat$gene)
    
    rownames(dat) = dat$patient
    dat$patient = NULL

    #split patients 
    med = median(dat$gene)
    
  	#remove NAs
  	z = which(is.na(dat$OS.time))
  	if(!(length(z) ==0)){
  	dat = dat[-z,]}

	  med_gene = median(dat$gene)	
	  dat$med = ""
	  if(med_gene ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat$gene > 0)
        l2 = which(dat$gene ==0)
        dat$med[l1] = 1
        dat$med[l2] = 0
        }

    if(!(med_gene ==0)){
        l1 = which(dat$gene >= med_gene)
        l2 = which(dat$gene < med_gene)
        dat$med[l1] = 1
        dat$med[l2] = 0
    }

    check1 = table(dat$med)[1] >= 10
    check2 = table(dat$med)[2] >= 10

    if(check1 & check2){
      if(dim(table(dat$med)) ==2){
  	  dat$gene = NULL
      dat = dat[,c("med", colnames(dat)[2:ncol(dat)-1])] #rearrange
      res.cox <- coxph(Surv(OS.time, OS) ~ ., data = dat)
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

#file = paste(date, "009_lncRNAs_prognosis_all_cancers.rds", sep="_")
#saveRDS(all_cancers_genes_surv_comb, file=file) #<---- important file 

##############RUN-----------------------------------------------------------------------------------

all_cancers_genes_surv_comb = readRDS("2019-04-16_009_lncRNAs_prognosis_all_cancers.rds") #load most recent file 
colnames(canc_conv)[2] = "canc"
all_cancers_genes_surv_comb = merge(all_cancers_genes_surv_comb, canc_conv, by="canc")

all_cancers_genes_surv_comb[,c(3:10)] = apply(all_cancers_genes_surv_comb[,c(3:10)], 2, function(x){as.numeric(x)})
all_cancers_genes_surv_comb = as.data.table(all_cancers_genes_surv_comb)

###-------------------------------------------------------------------------------------------------

all_cancers_genes_surv_comb$pval = -log10(as.numeric(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb$fdr = -log10(all_cancers_genes_surv_comb$fdr)
all_cancers_genes_surv_comb$HR = as.numeric(all_cancers_genes_surv_comb$HR)

z = which(is.na(all_cancers_genes_surv_comb$pval))
if(!(length(z)==0)){all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]}

z1 = which(all_cancers_genes_surv_comb$fdr == "Inf")
z2 = which(all_cancers_genes_surv_comb$upper95 == "Inf")
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-c(z1,z2),]

z = which(all_cancers_genes_surv_comb$HR > 20)
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

lineval = -log10(0.1)

all_cancers_genes_surv_comb$fdrsig = ""
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr < lineval] = "FDRnotSig"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval > lineval] = "Significant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval < lineval] = "NotSignificant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr >= lineval] = "FDRSignificant"

#check overlap of sig prognostic lncRNAs bewteen cancer types
all_cancers_genes_surv_comb$risk_perc = as.numeric(all_cancers_genes_surv_comb$risk_size)/as.numeric(all_cancers_genes_surv_comb$num_patients)
all_cancers_genes_surv_comb$risk_perc_tag[(all_cancers_genes_surv_comb$risk_perc > 0.48) | (all_cancers_genes_surv_comb$risk_perc < 0.52)] = "75%_more_risk_group"
all_cancers_genes_surv_comb$risk_perc_tag[all_cancers_genes_surv_comb$risk_perc > 0.75] = "75%_more_risk_group"

sig_lncs = as.data.table(all_cancers_genes_surv_comb)

#--------------------------------------------------------------------------------------------------
#ONLY FDR < 0.1 SIGNIFICNAT# 
#--------------------------------------------------------------------------------------------------

sig_lncs = as.data.table(filter(all_cancers_genes_surv_comb, fdr >= -log10(0.1)))
saveRDS(sig_lncs, file="3671_prognostic_lncRNAs_fdr0.1.rds")

all_cancers_genes_surv_comb = sig_lncs

#summary freq
sig_lncs = as.data.table(table(sig_lncs$canc, sig_lncs$gene))
sig_lncs = as.data.table(filter(sig_lncs, N > 0))
colnames(sig_lncs)[1] = "canc"
sig_lncs = merge(sig_lncs, canc_conv, by="canc")

#calculate how many lncRNAs overlap between cancer types 
cancs = unique(sig_lncs$type) #note these are FDR significant 

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

#SUMMARY OF HOW MANY PROGNOSTIC LNCRNAS IN COMMON BETWEEN CANCER TYPES 
pdf("overlap_ALL_sig_lncRNA_cands_bw_cancers_aug28.pdf", width=8, height=5)
g = ggplot(all_canc_pairs, aes(canc1, canc2)) +
  geom_tile(aes(fill=N)) +
  geom_text(aes(label = N), size=1.5) +
  scale_fill_gradient(low = "grey", high = "orange", na.value = 'transparent') +
    xlab("Cancer 1") + ylab("Cancer 2") + theme_bw()
ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend.title="# lncRNAs \noverlap")

dev.off()

#order by most significant to least significant 
order = as.data.table(table(all_cancers_genes_surv_comb$canc, all_cancers_genes_surv_comb$fdrsig))
order = order[order(-V1,N)]

z = order[order(order$V1, -order$N),]
# Remove duplicates
z1 = z[!duplicated(z$V1),]
#order again
z1 = z1[order(-V2,N)]
order = z1$V1

all_cancers_genes_surv_comb$canc <- factor(all_cancers_genes_surv_comb$canc, levels = order)
all_cancers_genes_surv_comb$canc  # notice the changed order of factor levels

all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR >= 1] = "Unfavourable"
all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR < 1] = "Favourable"

#Variation 2 of survival overview plot
head(all_cancers_genes_surv_comb)
#all_cancers_genes_surv_comb$HR = log2(all_cancers_genes_surv_comb$HR)

#summarize number favourable and unfabourable lcnRNAs by fdr significance per cancer type
all_cancers_genes_surv_comb = as.data.table(filter(all_cancers_genes_surv_comb, fdr >= -log10(0.1)))

#get order of cancer types by total number of lncRNAs 
order = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig))
order = as.data.table(dplyr::filter(order, N >0))
order = order[order(V2, -N)]

order = unique(order$V1)

summ = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig,
  all_cancers_genes_surv_comb$risk))
colnames(summ) = c("Cancer", "Sig", "Risk", "N")
summ = as.data.table(dplyr::filter(summ, N > 0))

#barplot----summary
#only include significant ones
#how many significant favourable vs unfavourable 
summ$Cancer = factor(summ$Cancer, levels = order)
summ$Risk = factor(summ$Risk, levels = c("Unfavourable", "Favourable"))

#how many unique lncrnas per cancer type
unique = as.data.table(table(all_cancers_genes_surv_comb$canc, all_cancers_genes_surv_comb$gene))
unique = unique[order(N)]
unique = as.data.table(filter(unique, N >0))
unique_sum = as.data.table(table(unique$V2))
unique_sum = unique_sum[order(N)]
filter(unique_sum, N ==1)
unique_only_one = as.data.table(filter(unique_sum, N ==1))

colnames(unique_only_one)[1] = "gene"
colnames(unique)[2] = "gene"

unique_sum = merge(unique_only_one, unique, by="gene")
table(unique_sum[,3])


######################################
#FIGURE 1B PART 1---------------------
######################################

pdf("final_figure_1B.pdf", height=6, width=6)
g = ggbarplot(summ, "Cancer", "N",
          fill = "Risk", color = "Risk", 
          palette = "npg")
ggpar(g, 
      font.xtickslab = c(9,"plain", "black"),
      xtickslab.rt = 45) + labs(x="Cancer type", y="Number of prognostic lncRNAs") + 
ggtitle("Number of Univariate Significant lncRNAs, adjusted CoxPH p-val < 0.1")+
scale_y_continuous(breaks=seq(0,2250,250))
dev.off()

#---------------------------------------------------------------------------------
### Figure 1 part 2 - correlations between prognostic lncRNAs in each cancer type 
#---------------------------------------------------------------------------------

### Get corrplot of survival candidates within each cancer type for significant lncRNAs 
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

get_summary = function(cancer){
  print(cancer)
  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))

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

#get correlations 

get_pairs_results = function(cancer){
  print(cancer)

  #get sig genes
  canc_genes = as.data.table(dplyr::filter(all_cancers_genes_surv_comb, canc == cancer))

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
    #res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05))
    sig_pairs = nrow(res2)  
    #check if lncRNA-lncRNA correlations match HRs 
    check_dir = function(lnc1, lnc2){
      
      hr_lnc1 = canc_genes$HR[canc_genes$gene == lnc1]
      hr_lnc2 = canc_genes$HR[canc_genes$gene == lnc2]

      check1 = ((hr_lnc1 > 1) & (hr_lnc2 > 1)) 
      check2 = ((hr_lnc1 < 1) & (hr_lnc2 < 1)) 
      
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
    t = t[order(n)]
    t$total_sig_pairs = sig_pairs
    t$total_pairs = tot_pairs
    t$perc = t$n/sig_pairs
    t$cancer = cancer
    res2$cancer = cancer
    return(res2)
  }
}

canc_results_pairs_types = llply(cancers, get_pairs_results, .progress = "text")

#save 
saveRDS(canc_results_pairs_types, file="correlation_lnc_lnc_results_april10_res2.rds")

#remove null
canc_results_pairs_types2 = Filter(Negate(is.null), canc_results_pairs_types)
canc_results_pairs_types2 = ldply(canc_results_pairs_types2)
canc_results_pairs_types2 = as.data.table(canc_results_pairs_types2)
colnames(canc_conv)[2] = "cancer"
canc_results_pairs_types2 = merge(canc_results_pairs_types2, canc_conv, by="cancer")

canc_results_pairs_types2$HR_pair = ""
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$match == "F"] = "Both \nFavourable"
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$match == "U"] = "Both \nUnfavourable"
canc_results_pairs_types2$HR_pair[canc_results_pairs_types2$match == "D"] = "Opposite \nHRs"

#keep only fdr significant ones
canc_results_pairs_types2 = as.data.table(filter(canc_results_pairs_types2, fdr < 0.05, abs(cor)))

#cancer order keep same as first plot
canc_results_pairs_types2$type <- factor(canc_results_pairs_types2$type, levels = rev(order))
canc_results_pairs_types2$column_name = paste(canc_results_pairs_types2$HR_pair, canc_results_pairs_types2$Exp_pair)

saveRDS(canc_results_pairs_types2, file="correlation_lnc_lnc_results_april10_res2.rds")

canc_results_pairs_types2 = readRDS("correlation_lnc_lnc_results_april10_res2.rds")

######################################
#FIGURE 1B PART 2---------------------
######################################

pdf("final_figure_1B_parttwo.pdf", width=4, height=4)

canc_results_pairs_types2$HR_pair = factor(canc_results_pairs_types2$HR_pair, levels = c("Both \nUnfavourable", "Opposite \nHRs", "Both \nFavourable"))

# Change density plot fill colors by groups
g = ggplot(canc_results_pairs_types2, aes(x=cor, fill=HR_pair), color="black") +
  geom_density(alpha=0.4, aes(x=cor, y=..density..)) + xlab("Spearman Correlation") + scale_fill_brewer(palette="Set1") +
  theme(legend.position="bottom")

ggpar(g, 
      font.tickslab = c(9,"plain", "black"), font.legend=c(8, "plain", "black"))

dev.off()

######################################
#FIGURE 1C PART ---------------------
######################################

gtex_res_risk = readRDS("lncRNAs_risk_groups_correlation_ranks.rds")
gtex_res_risk = ldply(gtex_res_risk)
gtex_res_risk$V2 = as.numeric(gtex_res_risk$V2)
summary(gtex_res_risk$V2)
gtex_res_risk = as.data.table(gtex_res_risk)
colnames(gtex_res_risk)[8:9] = c("Spearman_rho", "Spearman_p")
gtex_res_risk$Spearman_fdr = p.adjust(gtex_res_risk$Spearman_p, method="fdr")

gtex_res_risk$wilcox_p = as.numeric(gtex_res_risk$wilcox_p)
gtex_res_risk$wilcox_p = p.adjust(gtex_res_risk$wilcox_p, method="fdr")
gtex_res_risk = as.data.table(filter(gtex_res_risk, wilcox_p < 0.05))

#gtex_res = readRDS("significant_GTEX_comparisons_april10.rds")
#z = which(gtex_res$canc %in% c("Glioblastoma multiforme", "Brain Lower Grade Glioma"))
#brain  = gtex_res[z,]
#gtex_res = gtex_res[-z,]
#z = which(brain$tis %in% c("Brain - Cerebellum", "Brain - Spinal cord (cervical c-1)"))
#brain = brain[z,]
#gtex_res = rbind(gtex_res, brain)

#get lncRNAs that are prognosic with gtex data 
gtex_res = gtex_res_risk
head(gtex_res)
gtex_res$combo2 = paste(gtex_res$lnc, gtex_res$canc, sep="_")
head(all_cancers_genes_surv_comb)
all_cancers_genes_surv_comb$combo2 = paste(all_cancers_genes_surv_comb$gene, all_cancers_genes_surv_comb$type, sep="_")

gtex_res = merge(gtex_res, all_cancers_genes_surv_comb, by="combo2")
gtex_res$median_diff = as.numeric(gtex_res$median_diff)
gtex_res = as.data.table(filter(gtex_res, abs(median_diff) >= 0.2))

###Data
#gtex = readRDS("allGTEX_lncRNAs_scored_Feb2619.rds")
#tcga = readRDS("TCGA_all_lncRNAs_cancers_scored_byindexMay23.rds")
#gtex = gtex[,c(1:4, 6, 5)]

gtex_res$Hazard = gtex_res$risk.y
gtex_res$Hazard = factor(gtex_res$Hazard, levels = c("Unfavourable", "Favourable"))
gtex_res$Spearman_rho = as.numeric(gtex_res$Spearman_rho)
gtex_res$med[gtex_res$median_diff > 0] = "upreg"
gtex_res$med[gtex_res$median_diff < 0] = "downreg"
table(gtex_res$med, gtex_res$Hazard)

pdf("final_figure_1C.pdf", height=5, width=6)

#unfav 
m <- ggplot(gtex_res, aes(x = median_diff, y = HR)) +
 geom_point(aes(colour = Hazard)) + scale_color_npg() 
m + geom_density_2d(colour="black") + geom_hline(yintercept=0, linetype="dashed", color = "black")+
xlab("Median(High risk rank - GTEx rank)") + ylab("Hazard Ratio")+
geom_vline(xintercept=0, linetype="dashed", color = "black") #+ scale_colour_gradient(low = "blue", high = "red")

dev.off()





