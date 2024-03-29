###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)
library(mixtools)

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

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
rna = rna[-which(rna$Cancer %in% canc_rm),]

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer", "type", "OS", "OS.time"))
#all = all[,1:25170]

#lncRNA candidates, n = 166, n=173 combos 
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#look at only those with at least 10 risk patients 
#remove those with -Inf hazard ratios 
#then do spearman correlation of Hazard Ratios 

t = as.data.table(filter(allCands, data == "TCGA"))
p = as.data.table(filter(allCands, data == "PCAWG"))
p$num_risk = as.numeric(p$num_risk)
p = as.data.table(filter(p, num_risk >= 5)) #--> 44 (LUSC excluded bc was actually all patients in pcawg)

allCands = rbind(t,p)

z = which(p$upper95 == "Inf")
if(!(length(z)==0)){
p = p[-z,]}

p$HR = as.numeric(p$HR)
p$perc_risk = as.numeric(p$perc_risk)
p$perc_risk[(p$perc_risk >0.45) & (p$perc_risk <0.55)] = 0.5
p$perc_type[!(p$perc_risk==0.5)] = "unbal"
p$perc_type[(p$perc_risk==0.5)] = "bal"
colnames(p)[23] = "pcawg_perc_type"
colnames(p)[3] = "pcawg_hr"

#p = as.data.table(filter(p, pval <= 0.25))

t = as.data.table(filter(t, combo %in% p$combo))
t$perc_risk = as.numeric(t$perc_risk)
t$perc_risk[(p$perc_risk >0.45) & (t$perc_risk <0.55)] = 0.5
t$perc_type[!(p$perc_risk==0.5)] = "unbal"
t$perc_type[(p$perc_risk==0.5)] = "bal"

tp = merge(t, p , by=c("gene", "cancer", "combo", "CAT_geneName"))
tp$HR = as.numeric(tp$HR)
tp$HR = log2(tp$HR)
tp$pcawg_hr = log2(tp$pcawg_hr)
tp$fav[tp$HR > 0] = "unfav"
tp$fav[tp$HR < 0] = "fav"

#change cancer type
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2]= "cancer"
tp = merge(canc_conv, tp, by = "cancer")


mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")


pdf("summary_PCAWG_TCGA_HR_correlation_all_p_all_44_lncs.pdf", height=5, width=6)
g = ggscatter(tp, x = "HR", y = "pcawg_hr", 
   color = "black", fill="type", shape = 21, size = 2, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   #conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.sep = "\n")
   ) + geom_label_repel(data=subset(tp, pval.y <= 0.05), aes(label=CAT_geneName), size=2)+
   scale_fill_manual(values = sample(mypal, 10)) 
ggpar(g, legend="bottom") + xlab("TCGA Hazard Ratio") + ylab("PCAWG Hazard Ratio")+
geom_hline(yintercept=0, linetype="dashed", color = "black")+
geom_vline(xintercept=0, linetype="dashed", color = "black")


dev.off()






