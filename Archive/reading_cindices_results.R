
source("universal_LASSO_survival_script.R")


files = list.files(pattern = "cindices_1000CV_1000_no_fdr_ELASTICNET.rds")

res = as.data.frame(matrix(ncol=3)) ; colnames(res) = c("cindex", "canc", "type")

for(i in 1:length(files)){
	f = readRDS(files[[i]])
	res = rbind(res, f)
}

res = res[-1,]

library(ggpubr)

pdf("cindices_results_1000CVs_FDR0.2_cancers.pdf", width=10, height=10)
g = ggboxplot(res, x="type", y="cindex", facet.by="canc", color="type", order=c("cinds_clin", "cinds_justlncs", "cinds_combined"))
g + geom_hline(yintercept=0.5, linetype="dashed", color = "red") + 
stat_summary(fun.y=median, geom="line", aes(group=1))  + 
stat_summary(fun.y=median, geom="point")+
theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
theme(strip.text.x = element_text(size = 8, colour = "black"))
dev.off()


#reading genes results 

###change to be universal 

genes = ""


kirc_features = as.data.table(table(unlist(kirc_genes_results)))
kirc_features = kirc_features[order(N)]
kirc_features = dplyr::filter(kirc_features, N >=500)
kirc_features$canc = "kirc"
kirc_features$name = ""
for(i in 1:nrow(kirc_features)){
	z = which(fantom$gene == kirc_features$V1[i])
	kirc_features$name[i] = fantom$CAT_geneName[z]
}
