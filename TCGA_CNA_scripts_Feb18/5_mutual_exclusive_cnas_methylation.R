#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(parallel)
library(limma)
library(EnvStats)

mypal = pal_npg("nrc", alpha = 0.7)(10)


###DATA--------------------------------------------------------

cna_df = readRDS("all_dfs_CNAs_status_exp_patients_sept12.rds")
meth_df = readRDS("all_dfs_methylaion_status_exp_patients_sept12.rds")
meth_df$combo = paste(meth_df$gene, meth_df$cancer, sep="_")

#these have both methylation and cnas?
meth_df = merge(meth_df, cna_df, by=colnames(meth_df)[which(colnames(meth_df) %in% colnames(cna_df))])

#12 combos have patients with both CNAs and methylation 

#x-axis --> combo
#y-axis --> methlation vs cna 
#fill --> type of meth or 

meth_df$probe = NULL

meth = meth_df
meth$cna_status = NULL
meth$median = NULL
meth$test = "Methylation"


cna = meth_df
cna$mut_status = NULL
cna$median = NULL
cna = cna[,c(1:4, 7, 5,6, 8)]
colnames(cna)[5] = "mut_status"
cna$test = "CopyNumber"

exp = meth_df
exp$mut_status = NULL
exp$cna_status = NULL
exp = exp[,c(1, 3:5, 2, 6,7,8)]
colnames(exp)[5] = "mut_status"
exp$test = "Expression"

all_ppl = rbind(cna, meth, exp)

#split by gene 
all_ppl_list = split(all_ppl, by="gene")

get_overlap = function(gene_dat){

z = which(gene_dat$mut_status %in% c("noCNA"))
if(!(length(z)==0)){gene_dat = gene_dat[-z,]}

gene_dat$test = factor(gene_dat$test, levels=c("Expression", "CopyNumber", "Methylation"))
gene_dat$mut_status = factor(gene_dat$mut_status, levels=c("High", "Low", "DUP", "noCNA", "DEL", "Methylated", "Unmethylated"))

g = ggplot(gene_dat, aes(patient, test)) +
  
  geom_tile(aes(fill = mut_status)) + ggtitle(paste(gene_dat$name[1], gene_dat$cancer[1])) + 
    
    scale_fill_manual(values=c("#D6CAC3", "#E0A05B" ,"#B75BDB", "#9ADEB1" ,"#D485B5", "#AAE361", "#88AFD8")) +
    
    xlab("Patients") + ylab("") +    
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#get chisq value
gene_ch = as.data.table(filter(meth_df, gene == gene_dat$gene[1]))
t = table(gene_ch$mut_status, gene_ch$cna_status)
print(g)
return(chisq.test(t))

}

pdf("plots_lncRNAs_wSig_CNAs_and_Mehtylation.pdf", width=9)
llply(all_ppl_list, get_overlap)
dev.off()


