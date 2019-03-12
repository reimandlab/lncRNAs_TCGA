setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/full_gene_exp_perm")

library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)

real_file = readRDS("TCGA_results_multivariate_results_Oct3.rds")

cands = read.csv("168_lncRNA_cancers_combos_22_cancer_types_feb19.csv")
cands = as.data.table(cands)
cands = as.data.table(filter(cands, data == "TCGA"))
t = as.data.table(table(real_file$cancer))
t = t[order(N)]

#read in all files 
results = list.files(pattern="_.rds")
print(length(results))

#break into cancer types 
get_canc = function(file){
  dat = readRDS(file)
  #add round id 
  round = unlist(strsplit(file, "_"))[2]
  dat$round_id = round
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
rounds = unique(all_res$round_id)
all_res$round_id = as.numeric(all_res$round_id)
rounds = as.numeric(rounds)

get_fdr_v2 = function(round){

  round = as.character(round)
  z = which(all_res$round_id == round)

  canc_dat = as.data.table(all_res[z,])
  canc_dat$fdr = p.adjust(canc_dat$wald_p, method="fdr")
  rounds = length(unique(canc_dat$round))
  genes = length(unique(canc_dat$gene))
  fdrs = length(which(canc_dat$fdr < 0.05))
  total = nrow(canc_dat)
  row = c(rounds, genes, fdrs, total)
  names(row) = c("rounds", "genes", "fdrs_sig", "total")
  print(row)

}

rounds_fdrs = llply(rounds, get_fdr_v2)
rounds_fdrs = ldply(rounds_fdrs)

get_fdr = function(canc){

  canc_dat = as.data.table(filter(all_res, cancer == canc))
  canc_dat$fdr = p.adjust(canc_dat$wald_p, method="fdr")
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

#write.csv(l, file="results_bonferoni_permutations_march7.csv", quote=F, row.names=F)
#write.csv(all_res,file="all_bonferoni_res_results_permutations.csv", quote=F, row.names=F)

#write.csv(l, file="results_clin_shuffle_fdrs_match10th.csv", quote=F, row.names=F)

#summarize in how many rounds each gene appeared 
#expectation that it would appear many times? 

get_gene_counts = function(canc){

  canc_dat = as.data.table(filter(all_res, cancer == canc))
  canc_dat$fdr = p.adjust(canc_dat$wald_p, method="fdr")
  rounds = length(unique(canc_dat$round))
  genes = length(unique(canc_dat$gene))
  fdrs = length(which(canc_dat$fdr < 0.05))
  total = nrow(canc_dat)

  #sum of genes
  new = as.data.table(table(canc_dat$round))
  new = new[order(N)]
  #t$total_rounds = rounds
  new$canc = canc
  new$type = "permutation" 
  colnames(new) =c("Round_id", "num_lncRNAs_selected_in_round", "cancer", "type")

  return(new)
}

res = llply(cancers, get_gene_counts)
res = ldply(res)
res = as.data.table(res)
res = res[order(-num_lncRNAs_selected_in_round)]

t = as.data.table(table(cands$cancer))
t = t[order(N)]

plot_real_vs_permutation = function(canc){

  canc_dat = as.data.table(filter(all_res, cancer == canc))

  real_val = as.data.table(filter(t, V1 == canc))
  real_val$Round_id = "roundid"
  colnames(real_val)[2] = "num_lncRNAs_selected_in_round"
  colnames(real_val)[1] = "cancer"
  real_val$type = "real"
  cols = colnames(canc_dat)
  real_val = real_val[,..cols]

  #combine real and permutation
  canc_dat = rbind(canc_dat, real_val)

  #plot
  #x-axis: rounds 1-100 for example
  #y-axis: num of lncRNAs selected at each round 
  canc_dat = canc_dat[order(num_lncRNAs_selected_in_round)]
  canc_dat$rounds = 1:nrow(canc_dat)

  ggdotplot(canc_dat, x = "rounds", y = "num_lncRNAs_selected_in_round", size=0.4,
    color = "type", fill = "type", title=canc) 

  #Tobs <- 1.837
  #hist(Tstar, labels=T)
  #abline(v=Tobs, lwd=2, col="red")
  #plot(density(Tstar), main="Density curve of Tstar")
  #abline(v=Tobs, lwd=2, col="red")

}


pdf("summary_elastic_net_selected_real_vs_permutations.pdf", width=10, height=9)
llply(cancers, plot_real_vs_permutation)
dev.off()














