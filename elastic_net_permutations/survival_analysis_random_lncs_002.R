source("check_lnc_exp_cancers.R")

source("permutation_universal_LASSO_survival_script.R")
print("done source script")

#need to combine randomly generated datasets back into an "all" file
#matrix with all gene expression 

#canc_datas #<-- into one matrix 

#1. bind all randomly shuffled cancer datasets into one file 

dff <- ldply(canc_datas, data.frame)
dim(dff)

all = dff
all_permuted = all
all = all_permuted

#load elastic net cross-validation results 

#cindices
res = readRDS("permutations_elastic_net_all_cancers_Jan22.rds")
res = as.data.table(res)
res$type_test = "permutation"
realres = readRDS("noFDR_all_cindices_june22.rds")
realres = as.data.table(realres)
realres$type_test = "real"

#how do cindices from random shuffled samples compare to real survival information 
#for each cancer type, compare combo cindex real to shuffled, lncRNA only, clinical only 
#expect real data to perform better 

cancers = as.list(as.character(unique(res$canc)))

###########
#function 1 - compare cindex via wilcoxon and plot boxplot for each cancer type
###########

#pdf("results_permutations_cindices_vs_real.pdf")
get_comp = function(cancer){
     canc_dat_real = as.data.table(filter(realres, canc == cancer))
    canc_dat_perm = as.data.table(filter(res, canc == cancer))
    all_dat = rbind(canc_dat_real, canc_dat_perm)
    
    #look at the three types of tests seperatley (lncs only, clinical only, combined)
    tests = unique(all_dat$type)
    
    my_comparisons <- list( c("cinds_justlncs", "cinds_clin"), c("cinds_clin", "cinds_combined"), c("cinds_justlncs", "cinds_combined"))
    
    #boxplot 
    g = ggboxplot(all_dat, x="type_test", y="cindex", add = "jitter", facet.by = "type", color="type_test", title = cancer) +
      geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
      stat_compare_means(label = "p.format") +
      stat_n_text()
      
    print(g)
    
    print(paste(cancer, "done"))
    
    #wilcoxon test 
}

#llply(cancers, get_comp)
#dev.off()

#how many lncRNAs were selected by elastic net in the random permutations and what is their performance? 

real_genes = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")
perm_genes = readRDS("permutations_genes_keep_100CV_No_FDR_Jan22nd2018.rds")

real_genes$type = "real"
perm_genes$type = "perm"
all_genes = rbind(real_genes, perm_genes)
all_genes$combo = paste(all_genes$GeneName, all_genes$Cancer)

###########
#function 2 - compare num of genes selected in real dataset and in permutations 
###########

get_genes = function(cancer){
 
  all_dat = as.data.table(filter(all_genes, Cancer == cancer))
  summ = as.data.table(table(all_dat$type))
  
  #barplot
  g = ggbarplot(summ, x = "V1", y = "N", title = cancer, 
                label = TRUE, label.pos = "out") + labs(x="Elastic Net", y = "Num of genes")
  
  print(g)
  
  print(paste(cancer, "done"))
  
  #wilcoxon test 
}

#RUN 

cancers = as.list(as.character(unique(all_genes$Cancer)))

#pdf("results_permutations_num_genes_selected_vs_real.pdf")
#llply(cancers, get_genes)
#dev.off()


###########
#function 3 - get survival performance of eahc of the cands selected in permutatins 
###########

all_genes = as.data.table(all_genes)
perms = as.data.table(filter(all_genes, type == "perm"))

#list of gene-cancer combos 
perms$comb = paste(perms$Geneid, perms$Cancer, sep="_")
combos = as.list(unique(perms$comb))

get_km_plot = function(comb){
  print(comb)
  gene = unlist(strsplit(comb, "_"))[1]
  canc = unlist(strsplit(comb, "_"))[2]
  cancer = canc_conv$type[canc_conv$cancer == canc]

  all_g = all
  all_g = as.data.frame(all_g)
  dat = all[,c(which(colnames(all) %in% c("type", gene, "OS", "OS.time")))]

  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  dat[,1] = as.numeric(dat[,1])
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
  wald_p = round(glance(cox_mod)$p.value.wald, digits=4)
  hr = round(summary(cox_mod)$coefficients[2], digits=4)

  results = c(cancer, gene, gene_name, canc, conc, wald_p, hr)
  names(results) = c("type", "gene", "gene_name", "cancer", "concordance", "wald_p", "hr")

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
    return(results)
    print(paste("done", comb))
}
} 

pdf("km_plots_permutations_results_lncRNAs.pdf", width=9)
perms_surv_results = llply(combos, get_km_plot, .progress="text")
dev.off()

#how many are significant prior to fdr? 

#228 unique genes, 31/228 sig before fdr 14.5%
perms_surv_results1 = as.data.table(ldply(perms_surv_results)) #35/241 
perms_surv_results1$wald_p = as.numeric(perms_surv_results1$wald_p)
sig = filter(perms_surv_results1, wald_p < 0.05)
perms_surv_results1$fdr = p.adjust(perms_surv_results1$wald_p, method="fdr")
fdr_sig = filter(perms_surv_results1, fdr < 0.05)

