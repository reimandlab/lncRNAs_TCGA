#------Load libraries and scripts-----------------------------------

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#this script below prepares the RNA and clinical files for analysis 
source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

source("check_lnc_exp_cancers.R")

###Data

#[1] outlier based results 
res = readRDS("TCGA_ION_CHANNEL_results_Jan3119.rds")
res$gene_name = sapply(res$gene, get_name_pcg)

#filter to keep at least 10 people in risk group 
res = as.data.table(res, num_risk >= 10)

res$HR = as.numeric(res$HR)
res$HR = log2(res$HR)

#keep only sig ones
res = as.data.table(filter(res, pval < 0.05))
res$type = ""
res$type[res$HR > 0] = "Unfavourable"
res$type[res$HR < 0] = "Favourable"

#how many ion channles are always unfavourable or always favourable?
t = as.data.table(table(res$gene, res$type))
t = as.data.table(filter(t, N > 0)) #250 unique genes
tt = as.data.table(table(t$V1))
tt = as.data.table(filter(tt, N ==1)) 

ttt = as.data.table(filter(t, V1 %in% tt$V1))
ttt = as.data.table(filter(ttt, V2 == "Unfavourable"))
table(ttt$N)

cands = c("GJB2", "CATSPER1", "AQP9", "SCN9A")
cands = sapply(cands, get_ensg_pcg)

ttt$gene = sapply(ttt$V1, get_name_pcg)
tttt = as.data.table(table(ttt$N))

#summarize uniqueness of genes 
file = paste("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/KI_plots", "Figure1e.png", sep="_")
png(file)
ggbarplot(tttt, x="V1", "N") + labs(x="Number of cancer types", y="Number of ion channels") + theme_bw()
dev.off()

#82 ion channels were unfavourable in some cancers and favourable in others  
ics_unfav = as.data.table(filter(res, gene %in% ttt$V1))

summ = as.data.table(table(ics_unfav$cancer))
colnames(summ) = c("Cancer", "Unfavourable_ion_channels")
summ = summ[order(-(Unfavourable_ion_channels))]

file = paste("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/KI_plots", "Figure1d.png", sep="_")
png(file, width = 2500, height = 2000, res = 500)
g = ggbarplot(summ, x="Cancer", "Unfavourable_ion_channels") +
theme_bw() +
labs(y="Number of unfavourable \n ion channels")
ggpar(g, xtickslab.rt = 45, font.tickslab = c(6,"plain", "black")) 
dev.off()

#include HR in figure 
orderr = summ$Cancer

ics_unfav$cancer = factor(ics_unfav$cancer, levels= orderr)

ttt = ttt[order(-N)]
orderr = ttt$gene
ics_unfav$gene_name = factor(ics_unfav$gene_name, levels= orderr)
ics_unfav$fdr = ""
ics_unfav$fdr[ics_unfav$fdr_pval < 0.05] = "*"

#file = paste("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/KI_plots", "Figure1Sup.png", sep="_")

#png(file, width = 3750, height = 2000, res = 450)

file = paste("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/KI_plots", "Figure1Sup.pdf", sep="_")

pdf(file, width = 14, height = 7)

g = ggplot(ics_unfav, aes(gene_name, cancer)) + theme_bw()+
  geom_tile(aes(fill = HR)) + 
	geom_text(aes(label = fdr), size=2) +
  scale_fill_gradient(low="burlywood1", high="red", na.value = 'transparent') + labs(x = "Cancer", y="Ion channel")

g = ggpar(g, xtickslab.rt = 90, font.xtickslab = c(6,"plain", "black"), font.ytickslab = c(9,"plain", "black"))+
theme(legend.position="none")

xplot = ggplot(ics_unfav, aes(gene_name)) + geom_bar(fill = "#0073C2FF") + theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

yplot = ggplot(ics_unfav, aes(cancer)) + geom_bar(fill = "#0073C2FF") + theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + coord_flip()

ggarrange(xplot, NULL, g, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(3, 1), heights = c(1, 3),
          common.legend = FALSE)

dev.off()

