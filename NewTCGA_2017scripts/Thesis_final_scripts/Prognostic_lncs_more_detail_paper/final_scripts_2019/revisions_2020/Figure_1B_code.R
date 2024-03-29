source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#summarize patients
pats = unique(rna[,c("type", "Cancer")])
tt = as.data.table(table(rna$type))
colnames(tt) = c("type", "num_patients")
tt = merge(tt, pats, by="type")
tt = tt[order(num_patients)]
write.csv(tt, file="/u/kisaev/Jan2021/SuppTable1_TCGA_cancer_types_used_in_study.csv", quote=F, row.names=F)

#------FEATURES-----------------------------------------------------

rna = as.data.frame(rna)
dim(rna)
dim(pcg)

table(rna$type)
table(pcg$type)

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

saveRDS(all_genes, file="all_genes_used_in_TCGA.rds")

###---------------------------------------------------------------

#function that tests each lncRNA's survival

cancers = as.list(unique(rna$Cancer))

#3. function that splits data into cancers
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. function that calculates survival for each gene

canc_survival_genes = function(dato){

  lncs_dat = which(str_detect(colnames(dato), "ENSG"))
  medians = apply(dato[,lncs_dat], 2, median)
  dat_meds = dato
  for(k in 1:length(medians)){
        #print(k)
        med = medians[k]
        k = which(colnames(dat_meds) == names(med))
        if(med ==0){
          #if median = 0 then anyone greater than zero is 1
          l1 = which(dat_meds[,k] > 0)
          l2 = which(dat_meds[,k] ==0)
          dat_meds[l1,k] = 1
          dat_meds[l2,k] = 0
        }

        if(!(med ==0)){
          l1 = which(dat_meds[,k] >= med)
          l2 = which(dat_meds[,k] < med)
          dat_meds[l1,k] = 1
          dat_meds[l2,k] = 0
        }
        } #end assigning labels to each patient for each gene

  sums = apply(dat_meds[,lncs_dat], 2, sum)

  #for the full cohort get lncRNAs that have expression in less than 15 patients
  #or less than 10% of patients
  perc_10 = dim(dato)[1]*0.1
  #how many lncRNAs have only 10% of cohort or 15 patients (whichever is greater)
  #with > 0 expression in the cohort
  z1 = which(sums < perc_10)
  print(perc_10)
  z2 = which(sums < 15)
  all_rm = names(sums)[unique(c(z1,z2))]

	print(all_rm)
	z = which(colnames(dato) %in% all_rm)
  dato = dato[,-z]
  print(dato$type[1])

  z = which(str_detect(colnames(dato), "ENSG"))
  genes = unique(colnames(dato)[z])

	canc_data_genes_analyze = dato

	get_survival = function(gene){
	  print(gene)
  	results_cox <- as.data.frame(matrix(ncol=10)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients", "ph_test", "lnc_only_pval")

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
      dat = dat[,c("med", colnames(dat)[2:ncol(dat)-1])]
      res.cox <- coxph(Surv(OS.time, OS) ~ ., data = dat)
      lnc_only_model = coxph(Surv(OS.time, OS) ~ med, data = dat)
      test.ph <- cox.zph(lnc_only_model)
      lnc_test_ph = test.ph$table[1,3]
      lnc_only_pval = summary(lnc_only_model)$coefficients[5]
    	hr = summary(res.cox)$coefficients[1,c(2)]
      num_pat = nrow(dat)
      if(hr > 1){
        risk = length(which(dat$med ==1))
      }
      if(hr <1){
        risk = length(which(dat$med ==0))
      }

      row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)], risk, num_pat, lnc_test_ph, lnc_only_pval)
     	names(row) <- names(results_cox)
    	return(row)
  	}}} #end get_survival function

	genes_survival = llply(genes, get_survival, .progress="text")
	genes_survival_res = ldply(genes_survival, rbind)
	#fdr
	colnames(genes_survival_res) = c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients", "lnc_ph_test", "lnc_only_pval")
	genes_survival_res$fdr = p.adjust(as.numeric(genes_survival_res$pval), method="fdr")
  genes_survival_res$fdr_lnc_only = p.adjust(as.numeric(genes_survival_res$lnc_only_pval), method="fdr")

	genes_survival_res$canc = dato$Cancer[1]
	genes_survival_res = as.data.table(genes_survival_res)
	genes_survival_res = genes_survival_res[order(fdr)]
	return(genes_survival_res)
}

#DO NOT RUN

#all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
#all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)
#saveRDS(all_cancers_genes_surv_comb, file="lncRNAs_all_survival_results_feb27.rds") #<---- important file

##############RUN-----------------------------------------------------------------------------------

all_cancers_genes_surv_comb = readRDS("lncRNAs_all_survival_results_feb27.rds")
all_cancers_genes_surv_comb = as.data.table(all_cancers_genes_surv_comb)
canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]
colnames(canc_conv)[2] = "canc"
all_cancers_genes_surv_comb = merge(all_cancers_genes_surv_comb, canc_conv, by="canc")
write.csv(all_cancers_genes_surv_comb, file="ALL_lncRNAs_survival_Feb262019.csv", quote=F, row.names=F)

all_cancers_genes_surv_comb = as.data.table(all_cancers_genes_surv_comb)

###-------------------------------------------------------------------------------------------------

#plot scatter plot - HR versus p-value draw line for FDR = 0.05
all_cancers_genes_surv_comb$pval_log10 = -log10(as.numeric(all_cancers_genes_surv_comb$pval))
all_cancers_genes_surv_comb$fdr_log10 = -log10(all_cancers_genes_surv_comb$fdr)
all_cancers_genes_surv_comb$HR = as.numeric(all_cancers_genes_surv_comb$HR)

z = which(is.na(all_cancers_genes_surv_comb$pval))
if(!(length(z)==0)){all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]}

z1 = which(all_cancers_genes_surv_comb$fdr == "Inf")
z2 = which(all_cancers_genes_surv_comb$upper95 == "Inf")
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-c(z1,z2),]

z = which(all_cancers_genes_surv_comb$HR > 10)
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

z = which(all_cancers_genes_surv_comb$HR < 0.1)
all_cancers_genes_surv_comb = all_cancers_genes_surv_comb[-z,]

lineval = -log10(0.05)

all_cancers_genes_surv_comb$fdrsig = ""
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr_log10 < lineval] = "FDRnotSig"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval_log10 > lineval] = "Significant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$pval_log10 < lineval] = "NotSignificant"
all_cancers_genes_surv_comb$fdrsig[all_cancers_genes_surv_comb$fdr_log10 >= lineval] = "FDRSignificant"

#check overlap of sig prognostic lncRNAs bewteen cancer types
all_cancers_genes_surv_comb$risk_perc = as.numeric(all_cancers_genes_surv_comb$risk_size)/as.numeric(all_cancers_genes_surv_comb$num_patients)
all_cancers_genes_surv_comb$risk_perc_tag[(all_cancers_genes_surv_comb$risk_perc > 0.48) | (all_cancers_genes_surv_comb$risk_perc < 0.52)] = "75%_more_risk_group"
all_cancers_genes_surv_comb$risk_perc_tag[all_cancers_genes_surv_comb$risk_perc > 0.75] = "75%_more_risk_group"

sig_lncs = as.data.table(all_cancers_genes_surv_comb)
sig_lncs = as.data.table(filter(all_cancers_genes_surv_comb, fdr_log10 >= -log10(0.05)))

#save only those that appear in only one cancer type
print(dim(sig_lncs))

t = filter(as.data.table(table(sig_lncs$gene)), N ==1)
print(length(t$V1)) #lncRNAs signiciant in only one cancer type
print(length(unique(sig_lncs$gene)))
length(t$V1) / length(unique(sig_lncs$gene))
#sig_lncs = as.data.table(filter(sig_lncs, gene %in% t$V1))
sig_lncs$risk_perc_tag[(sig_lncs$risk_perc > 0.48) | (sig_lncs$risk_perc < 0.52)] = "equal"
sig_lncs$risk_perc_tag[sig_lncs$risk_perc > 0.6] = "high_risk"
sig_lncs$risk_perc_tag[sig_lncs$risk_perc < 0.45] = "low_risk"

all_cancers_genes_surv_comb = sig_lncs

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

#SUMMARY OF HOW MANY PROGNOSTIC LNCRNAS IN COMMON BETWEEN CANCER TYPES

#pdf("overlap_ALL_sig_lncRNA_cands_bw_cancers_aug28.pdf", width=8, height=5)
#g = ggplot(all_canc_pairs, aes(canc1, canc2)) +
#  geom_tile(aes(fill=N)) +
#  geom_text(aes(label = N), size=1.5) +
#  scale_fill_gradient(low = "grey", high = "orange", na.value = 'transparent') +
#    xlab("Cancer 1") + ylab("Cancer 2") + theme_bw()
#ggpar(g,
# font.tickslab = c(8,"plain", "black"),
# xtickslab.rt = 45, legend.title="# lncRNAs \noverlap")

#dev.off()

#order by most significant to least significant
order = as.data.table(table(all_cancers_genes_surv_comb$canc, all_cancers_genes_surv_comb$fdrsig))
order = order[order(-V1,N)]

z = order[order(V1, -N)]
# Remove duplicates
z1 = z[!duplicated(z$V1),]
#order again
z1 = z1[order(-V2,N)]
order = z1$V1

all_cancers_genes_surv_comb$canc <- factor(all_cancers_genes_surv_comb$canc, levels = order)
all_cancers_genes_surv_comb$canc  # notice the changed order of factor levels

all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR > 1] = "Unfavourable"
all_cancers_genes_surv_comb$risk[all_cancers_genes_surv_comb$HR < 1] = "Favourable"

#Variation 2 of survival overview plot
head(all_cancers_genes_surv_comb)
all_cancers_genes_surv_comb$HR = log2(all_cancers_genes_surv_comb$HR)

#summarize number favourable and unfabourable lcnRNAs by fdr significance per cancer type
all_cancers_genes_surv_comb = as.data.table(filter(all_cancers_genes_surv_comb, fdr_log10 >= -log10(0.05)))

#get order of cancer types by total number of lncRNAs
order = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig))
order = as.data.table(dplyr::filter(order, N >0)) # <---
order = order[order(V2, -N)]

order = unique(order$V1)

summ = as.data.table(table(all_cancers_genes_surv_comb$type, all_cancers_genes_surv_comb$fdrsig,
  all_cancers_genes_surv_comb$risk))
colnames(summ) = c("Cancer", "Sig", "Risk", "N")
summ = as.data.table(dplyr::filter(summ, N > 0))

#barplot----summary
summ$Cancer = factor(summ$Cancer, levels = order)
summ$Risk = factor(summ$Risk, levels = c("Unfavourable", "Favourable"))

write.csv(all_cancers_genes_surv_comb, file="/u/kisaev/Jan2021/SuppTable3_survival_Associated_lncs.csv", quote=F, row.names=F)

######################################
#FIGURE 1B PART 1---------------------
######################################

write.table(summ, file="figure1B_data_table.txt", quote=F, row.names=F, sep="\t")

pdf("/u/kisaev/Jan2021/final_figure_1B.pdf", height=6, width=6)
g = ggbarplot(summ, "Cancer", "N",
          fill = "Risk", color = "Risk",
          palette = "npg") + ylim(c(0,1850))
g = ggpar(g, #yticks.by = 200,
      font.xtickslab = c(9,"plain", "black"),
      xtickslab.rt = 45) + labs(x="Cancer type", y="Number of prognostic lncRNAs") #+ scale_y_continuous(trans='log10')
print(g)
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

canc_results = llply(cancers, get_summary, .progress = "text")
#remove null
canc_results = Filter(Negate(is.null), canc_results)
canc_results = do.call(rbind.data.frame, canc_results)
colnames(canc_results) = c("cancer", "total_pairs", "sig_pairs", "perc")

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
    print(tot_pairs)
    #res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05))
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
    #t = as.data.table(tidy(t))
		t = as.data.table(t)

		t = t[order(N)] %>% filter(N >0)
    t$total_sig_pairs = sig_pairs
    t$total_pairs = tot_pairs
    t$perc = t$N/sig_pairs
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
canc_results_pairs_types2 = as.data.table(filter(canc_results_pairs_types2, fdr < 0.05)) #, abs(cor) >= 0.3))
table(canc_results_pairs_types2$HR_pair)/dim(canc_results_pairs_types2)[1]

#cancer order keep same as first plot
canc_results_pairs_types2$type <- factor(canc_results_pairs_types2$type, levels = rev(order))
canc_results_pairs_types2$column_name = paste(canc_results_pairs_types2$HR_pair, canc_results_pairs_types2$Exp_pair)

saveRDS(canc_results_pairs_types2, file="correlation_lnc_lnc_results_april10_res2.rds")
canc_results_pairs_types2 = readRDS("correlation_lnc_lnc_results_april10_res2.rds")

######################################
#FIGURE 1B PART 2---------------------
######################################

canc_results_pairs_types2$HR_pair = factor(canc_results_pairs_types2$HR_pair, levels = c("Both \nUnfavourable", "Opposite \nHRs", "Both \nFavourable"))
canc_results_pairs_types2$fdr_sig = ""
z = which(((canc_results_pairs_types2$fdr < 0.05) & (abs(canc_results_pairs_types2$cor) >= 0.3)))
canc_results_pairs_types2$fdr_sig[z]= "sig"

cols = RColorBrewer::brewer.pal(8, "Set1")

pdf("/u/kisaev/Jan2021/final_figure_1B_parttwoa_2019.pdf", width=5, height=5)
# Change density plot fill colors by groups
g1 = ggplot(canc_results_pairs_types2[canc_results_pairs_types2$HR_pair == "Both \nUnfavourable"], aes(x=cor, fill=HR_pair), color="black") +
  geom_density(alpha=0.4, aes(x=cor, y=..scaled..)) + xlab("") + scale_fill_manual(values=cols[1]) +
  theme(legend.position="bottom")+theme_classic()
g1 = ggpar(g1,
      font.tickslab = c(15,"plain", "black"), font.legend=c(4, "plain", "black"), xlim=c(-1,1))+
theme(legend.position="none")
print(g1)
dev.off()

pdf("/u/kisaev/Jan2021/final_figure_1B_parttwob_2019.pdf", width=5, height=5)

# Change density plot fill colors by groups
g2 = ggplot(canc_results_pairs_types2[canc_results_pairs_types2$HR_pair == "Opposite \nHRs"], aes(x=cor, fill=HR_pair), color="black") +
  geom_density(alpha=0.4, aes(x=cor, y=..scaled..)) + xlab("") + scale_fill_manual(values=cols[2]) +
  theme(legend.position="bottom")+theme_classic()

g2 = ggpar(g2,
      font.tickslab = c(15,"plain", "black"), font.legend=c(4, "plain", "black"), xlim=c(-1,1))+
theme(legend.position="none")
print(g2)
dev.off()

pdf("/u/kisaev/Jan2021/final_figure_1B_parttwoc_2019.pdf", width=5, height=5)
# Change density plot fill colors by groups
g3 = ggplot(canc_results_pairs_types2[canc_results_pairs_types2$HR_pair == "Both \nFavourable"], aes(x=cor, fill=HR_pair), color="black") +
  geom_density(alpha=0.4, aes(x=cor, y=..scaled..)) + xlab("Spearman Correlation") + scale_fill_manual(values=cols[3]) +
  theme(legend.position="bottom")+theme_classic()

g3 = ggpar(g3,
      font.tickslab = c(15,"plain", "black"), font.legend=c(4, "plain", "black"), xlim=c(-1,1))+
theme(legend.position="none")
print(g3)
dev.off()

pdf("/u/kisaev/Jan2021/final_figure2b_2019.pdf")
plot_grid(g1, g2, g3, ncol=1, nrow=3,  align = "v")
dev.off()
