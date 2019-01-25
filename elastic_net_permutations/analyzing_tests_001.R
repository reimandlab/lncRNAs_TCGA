library(data.table)
source("check_lnc_exp_cancers.R")

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

pdf("results_permutations_cindices_vs_real.pdf")
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

llply(cancers, get_comp)
dev.off()

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
llply(cancers, get_genes)
#dev.off()







