setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/march4_EN_perms")

library(data.table)
library(dplyr)
library(plyr)


cands = read.csv("168_lncRNA_cancers_combos_22_cancer_types_feb19.csv")
cands = as.data.table(cands)

#read in all files 
results = list.files()
print(length(results))

#break into cancer types 
get_canc = function(file){
  dat = readRDS(file)
  return(dat)
}

all_res = llply(results, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)
all_res$comb = paste(all_res$type, all_res$round, sep="_")


u = unique(all_res$comb)
u = sapply(u, function(x){unlist(strsplit(x, "_"))[1]})
u = as.data.table(table(u))
u = u[order(N)]
colnames(u) = c("cancer", "num_runs_successful_outof100")

#do fdr wtihin cancer type 
#calculate how many fdr sig in each one 

cancers = unique(all_res$cancer)

get_fdr = function(canc){

  canc_dat = as.data.table(filter(all_res, cancer == canc))
  canc_dat$fdr = p.adjust(canc_dat$wald_p, method="bonferroni")
  rounds = length(unique(canc_dat$round))
  genes = length(unique(canc_dat$gene))
  fdrs = length(which(canc_dat$fdr < 0.05))
  total = nrow(canc_dat)
  row = c(canc, rounds, genes, fdrs, total)
  names(row) = c("canc", "rounds", "genes", "fdrs_sig", "total")
  print(row)
  return(row)
}

l = llply(cancers, get_fdr)
l = ldply(l)
l = as.data.table(l)
l=l[order(rounds)]
l$fdr_sig_por = as.numeric(l$fdrs_sig)/as.numeric(l$total)
l=l[order(fdr_sig_por, rounds)]

#write.csv(l, file="results_permutations_march7.csv", quote=F, row.names=F)
#write.csv(all_res,file="all_res_results_permutations.csv", quote=F, row.names=F)

#save bonferoni version 

write.csv(l, file="results_bonferoni_permutations_march7.csv", quote=F, row.names=F)
write.csv(all_res,file="all_bonferoni_res_results_permutations.csv", quote=F, row.names=F)


#summarize in how many rounds each gene appeared 
#expectation that it would appear many times? 

get_gene_counts = function(canc){

  canc_dat = as.data.table(filter(all_res, cancer == canc))
  canc_dat$fdr = p.adjust(canc_dat$wald_p, method="bonferroni")
  rounds = length(unique(canc_dat$round))
  genes = length(unique(canc_dat$gene))
  fdrs = length(which(canc_dat$fdr < 0.05))
  total = nrow(canc_dat)

  #sum of genes
  t = as.data.table(table(canc_dat$round))
  t = t[order(N)]
  #t$total_rounds = rounds
  t$canc = canc
  colnames(t) =c("Round_id", "num_lncRNAs_selected_in_round", "cancer")

  return(t)
}

res = llply(cancers, get_gene_counts)
res = ldply(res)
res = as.data.table(res)
res = res[order(-num_lncRNAs_selected_in_round)]







