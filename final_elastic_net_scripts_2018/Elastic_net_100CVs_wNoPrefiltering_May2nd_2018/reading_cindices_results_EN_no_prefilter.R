
source("universal_LASSO_survival_script.R")

files = list.files(pattern = "cindices_100CV_May2_ElasticNet_no_prefilter.rds")

res = as.data.frame(matrix(ncol=3)) ; colnames(res) = c("cindex", "canc", "type")

for(i in 1:length(files)){
	f = readRDS(files[[i]])
	res = rbind(res, f)
}

res = res[-1,]

res = as.data.table(res)
t = as.data.table(table(res$canc))
t = filter(t, N >100)

res = filter(res, canc %in% t$V1)

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)

#to compare to get pvalues 
my_comparisons <- list( c("cinds_clin", "cinds_justlncs"), c("cinds_justlncs", "cinds_combined"), c("cinds_clin", "cinds_combined") )

pdf("cindices_results_elastic_net_100CVs_noprefilter_cancers_May2nd.pdf", width=12, height=12)
g = ggboxplot(res, x="type", y="cindex", facet.by="canc", color="type", order=c("cinds_clin", "cinds_justlncs", "cinds_combined"))

g =  g + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "cinds_clin") 

g + geom_hline(yintercept=0.5, linetype="dashed", color = "red") + 

stat_summary(fun.y=median, geom="line", aes(group=1))  + 
stat_summary(fun.y=median, geom="point")+

theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
theme(strip.text.x = element_text(size = 6.5, colour = "black"))

dev.off()


###reading genes results 

files = list.files(pattern = "lncRNAs_selected_100CV_May2_ElasticNet_no_prefiler.rds")

###change to be universal 

genes_keep = as.data.frame(matrix(ncol = 4)) ; colnames(genes_keep) = c("Geneid", "NumtimesChosen" , "Cancer", "GeneName")

pdf("summary_of_genes_results_elastic_net_100CVs_noprefilter_cancers_May2nd.pdf")

for(i in 1:length(files)){

f = files[[i]]
genes = readRDS(f)

genes_list = as.data.table(table(unlist(genes)))
genes_list = genes_list[order(N)]
#genes_list = dplyr::filter(genes_list, N >=50)
genes_list$canc = unlist(strsplit(f, "_"))[1]
genes_list$name = ""
for(i in 1:nrow(genes_list)){
	z = which(fantom$gene == genes_list$V1[i])
	genes_list$name[i] = fantom$CAT_geneName[z]
}

#draw histogram of how many times lncRNA got selected in 100 runs 
g = ggdotchart(genes_list, x = "name", y = "N",
           color = "canc",          
           palette = mypal, # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           #rotate = TRUE,                                # Rotate vertically
           #group = "cyl",                                # Order by groups
           dot.size = 1.3,                                 # Large dot size
           #label = round(genes_list$N),                        # Add mpg values as dot labels
           #font.label = list(color = "white", size = 1, 
           #                  vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(), 
           xlab = "lncRNAs selected by Elastic Net", 
           ylab = "% of runs selected (100 runs total)"                        # ggplot2 theme
           )

g = g + geom_text_repel(data = subset(genes_list, N > 50), aes(label = name), size=2) 

gg = g + scale_y_continuous(breaks = seq(0, 100, by = 10))
	
print(gg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1, size=4)) +
          geom_hline(yintercept = 50, linetype = 2, color = "blue"))

genes_list = dplyr::filter(genes_list, N >=50)
if(!(dim(genes_list)[1]==0)){
colnames(genes_list) = colnames(genes_keep) 
genes_keep = rbind(genes_keep, genes_list)
}
}
dev.off()

genes_keep = genes_keep[-1,]
genes_keep = as.data.table(genes_keep)

dups = genes_keep$GeneName[which(duplicated(genes_keep$GeneName))]
if(!(length(dups)==0)){
	z = which(genes_keep$GeneName %in% dups)
	for(y in 1:length(z)){
		genes_keep$GeneName[z][y] = paste(genes_keep$GeneName[z][y], y, sep="_")
	}	
}


#need to generate new palette 
n <- length(unique(genes_keep$Cancer))
mypal2 <- viridis_pal(option = "D")(n)  

#get actual counts and put beside dotchart 
t = as.data.table(table(genes_keep$Cancer))
t = t[order(N)]

canc_order = t$V1
genes_keep$Cancer = factor(genes_keep$Cancer, levels = canc_order)

pdf("summary_of_final_lncRNA_candidates_elastic_net_100CVs_noprefilter_cancers_May2nd.pdf", width=14, height=11)

overg = ggbarplot(t, x = "V1", y = "N",
          fill = "grey",               # change fill color by cyl
          palette = "mypal2", 
          color = "white",            # Set bar border colors to white
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90 ,
          ggtheme = theme_pubr(), 
          xlab = "Cancers",
          legend = "none",  
          ylab = "# of final lncRNAs selected more than 50%"           # Rotate vertically x axis texts
          )
overg = overg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=55, hjust=1, size=7)) +
          ggtitle("Summary of genes selected 100 CVs with No Pre-filtering")

#sumamrize final list of lncrnas 
g = ggdotchart(genes_keep, x = "GeneName", y = "NumtimesChosen", 
           color = "Cancer",         
           palette = "mypal2", # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           #rotate = TRUE,                                # Rotate vertically
           group = "Cancer",                                # Order by groups
           dot.size = 1.5,                                 # Large dot size
           #label = round(genes_keep$NumtimesChosen),                        # Add mpg values as dot labels
           #font.label = list(color = "white", size = 1, 
           #                  vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(), 
           xlab = "lncRNAs selected by Elastic Net", 
           ylab = "% of runs selected (100 runs total)"                        # ggplot2 theme
           )
gg = g + theme(legend.position = "bottom", legend.text = element_text(colour="black", size=7))
gg = gg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1, size=6)) +
          geom_hline(yintercept = 50, linetype = 2, color = "blue") + ggtitle("Summary of genes selected 100 CVs with No Pre-filtering")

#combine the two side by side
gg +overg + plot_layout(ncol = 1, heights = c(3, 1))

dev.off()

saveRDS(genes_keep, file="genes_keep_100CV_Noprefiltering_May2nd2018.rds")



