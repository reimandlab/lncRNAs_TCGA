setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/feb22_EN_perms")

library(data.table)
library(dplyr)
library(plyr)

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

all_res$fdr = p.adjust(all_res$wald_p, method="fdr")

#do fdr wtihin cancer type 
#calculate how many fdr sig in each one 

get_split = function(canc){

  z = which(all_res$type == canc)
  split = c()

  for(i in 1:(length(z)-1)){
    new = z[i+1]
    print(z[i])
    print(new)
    if(!(new-z[i] ==1)){
      print("stop")
      split = c(split, z[i+1])
    }
  }

  for(i in 1:length(split)){
    if(i ==1){
      start = split[i]-(split[i]-1)
      end = split[i]
      all_res[start:end,] = paste("round", i, sep="_")
    }
  }

}








