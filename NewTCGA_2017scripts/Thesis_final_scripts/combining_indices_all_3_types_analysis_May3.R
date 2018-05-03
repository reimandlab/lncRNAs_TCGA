source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)


#files of genes from 3 types of analysis 

#1. No prefiltering 
np = readRDS("genes_keep_100CV_Noprefiltering_May2nd2018.rds")
np$analysis = "noprefilter"

#2. Prefiltering 
pre = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")
pre$analysis = "noFDR"

#3. Prefiltering + FDR 0.2 
prefdr = readRDS("genes_keep_100CV_FDR0.2_May2nd2018.rds")
prefdr$analysis = "FDR0.2"

allgenes = rbind(np, pre, prefdr)

t = as.data.table(table(allgenes$Geneid, allgenes$Cancer, allgenes$analysis))
t = filter(t, N >0)
t$combo = paste(t$V1, t$V2, sep="_")
across = as.data.table(table(t$combo))
across = across[order(N)]

justonce = filter(across, N==1)
z = which(t$combo %in% justonce$V1)
once = t[z,]
once$freq = 1

twotypes = filter(across, N ==2)
z = which(t$combo %in% twotypes$V1)
twotwo = t[z,]
twotwo$freq = 2

threetypes = filter(across, N==3)
z = which(t$combo %in% threetypes$V1)
threethree = t[z,]
threethree$freq = 3

all_combined = rbind(once, twotwo, threethree)
#need to generate new palette 
n <- length(unique(all_combined$V2))
mypal2 <- viridis_pal(option = "D")(n)  

all_combined = as.data.table(all_combined)
all_combined = all_combined[order(-freq)]
myorder = unique(all_combined$V1)
all_combined$V1 <- factor(all_combined$V1, levels = myorder)
myorder = unique(all_combined$V2)
all_combined$V2 <- factor(all_combined$V2, levels = myorder)
	
pdf("summary_howmanytimes_across_types_filtering_candidates.pdf", width=20)
g = ggdotchart(all_combined, x = "V1", y = "freq",
          color = "V2", size = 3,
          #dd = "segment",
          sorting = "descending",  
          #group = "V2",  
          #add.params = list(color = "lightgray", size = 1.5),
          #position = position_dodge(0.3),
          palette = "mypal2",
          ggtheme = theme_pubr(),
          xlab = "Candidate lncRNAs",
          ylab = "Number of times selected across analysis" 
)
gg = g + theme(legend.position = "bottom", legend.text = element_text(colour="black", size=7))
gg = gg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1, size=6)) + ggtitle("Top Candidates - Selection across types of analysis and filtering")
print(gg)
dev.off()

all_combined$N = NULL
colnames(all_combined) = c("Name", "Cancer", "AnalysisType", "Combo", "freq")

saveRDS(all_combined, file="all_candidates_combined_cancers_typesAnalysis_May3rd.rds")

