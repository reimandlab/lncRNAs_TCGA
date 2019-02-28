setwd("/Users/kisaev/remote10")

library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(EnvStats)

#summmary results 

res = readRDS("results_analysis_Feb26.rds")
res = as.data.table(ldply(res))
res$combo = paste(res$canc, res$tis)

#save pairs for sup table 3

ss3 = unique(res[,c("canc", "tis")])
write.csv(ss3, file="sup_table3_gtex_tis_cancers_studied.csv", quote=F, row.names=F)

#keep only sig ones
res$median_difference = as.numeric(res$median_difference)
res = as.data.table(filter(res, fdr < 0.05, abs(median_difference) >= 0.25))

#figure 1b 
res$sign = ""
res$sign[res$median_difference > 0] = "upreg"
res$sign[res$median_difference < 0] = "downreg"


