
source("universal_LASSO_survival_script.R")

files = list.files(pattern = "cindices_1000CV_1000_no_fdr_ELASTICNET.rds")

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
library(wesanderson)
library(tidyr)


#to compare to get pvalues 
my_comparisons <- list( c("cinds_clin", "cinds_justlncs"), c("cinds_justlncs", "cinds_combined"), c("cinds_clin", "cinds_combined") )

#horizontal boxplot of all c-indices distributions 
res$cindex = as.numeric(res$cindex)

#remove NAs
z = which(is.na(res$cindex))
res= res[-z,]
#order cancer types by worse lncRNA performance to best lncRNA performance 
res = as.data.table(res)
lncs = filter(res, type =="cinds_justlncs")
lncs = aggregate(lncs[,1], list(lncs$canc), median)
lncs = as.data.table(lncs)
lncs = lncs[order(x)]

canc_order = unique(lncs$Group.1)
res$canc = factor(res$canc, levels = canc_order)

res$type = as.factor(res$type)
order = c("cinds_justlncs", "cinds_clin", "cinds_combined")
res$type = factor(res$type, levels = order)

saveRDS(res, file="noFDR_all_cindices_june22.rds")

#for each cancer type plot the distribution of cindices 
plots = res %>% dplyr::group_by(canc, type) %>% do(plots=ggplot(data=.) +
         aes(x=cindex) + geom_histogram())

#to view plots -----------------------------------------------------------
#DO NOT RUN
#plots$plots
#end 

#but basically distirbution is normal 
#so can calculate standard error 

#for horizontal plot
#just plot the cancer types that have median just_lncRNA greater than 0.5
#later can decide if only want to inlucdee those significanly better than clinical variables

#need to calcualte confidence interval 
cisum = as.data.table(res %>% dplyr::group_by(canc, type) %>% 
  dplyr::summarize(n = n(), sd = sd(cindex), mean = mean(cindex), median=median(cindex), low=quantile(cindex, 0.025), high=quantile(cindex, 0.975)))
cisum = cisum[order(-median, sd)]
cisum = as.data.table(dplyr::filter(cisum, type == "cinds_justlncs", low >=0.5, high>=0.5))
as.data.table(dplyr::filter(cisum, type == "cinds_justlncs", median >= 0.5))

summ = aggregate(res[, 1], list(res$canc, res$type), median)
summ = as.data.table(summ)
summ = filter(summ, Group.2 == "cinds_justlncs", cindex >= 0.5)

res = filter(res, canc %in% summ$Group.1)
canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]
colnames(canc_conv)[2] = "canc"
res = merge(res, canc_conv, by="canc")
colnames(res)[3:4] = c("type", "TYPE")
res = as.data.table(res)
colnames(summ)[1] = "canc"
summ = merge(summ, canc_conv, by="canc")
summ = as.data.table(summ)
summ = summ[order(cindex)]
res$canc <- factor(res$canc, levels =summ$canc)
res$TYPE <- factor(res$TYPE, levels = summ$type)

mypal = wes_palette("FantasticFox")

pdf("cindices_results_elastic_net_100CVs_No_FDR_cancers_July11_nofacet_horizontal_boxplots.pdf", width=11, height=10)

g = ggboxplot(res, "TYPE", "cindex", orientation = "horizontal", fill="type", color="black", palette=mypal, notch = TRUE)
g =  g + stat_compare_means(aes(group = type), label = "p.signif") + theme_minimal()
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red"))

dev.off()

#V2 using ggplot

ggplot(res, aes(TYPE, cindex)) +
geom_boxplot(aes(fill = type), color="black", outlier.size = 0.01) + theme_bw() + xlab("Cancer") + ylab("C-index") +
stat_compare_means(aes(group = type), label = "p.signif") + 
coord_flip() + scale_fill_manual(values=mypal) +
geom_hline(yintercept=0.5, linetype="dashed", color = "red")
dev.off()

#for each cancer type do wilcoxon test 
#keep those cancers in which, lncRNA median is greater than 0.5 (done)
#and those with higher c-index than clinical or just as good

cancers = unique(res$TYPE)
check_perform = function(cancer){
  dat = subset(res, TYPE==cancer)
  w = wilcox.test(dat$cindex[dat$type=="cinds_justlncs"], dat$cindex[dat$type=="cinds_clin"], alternative="greater")
  med_lnc = median(dat$cindex[dat$type=="cinds_justlncs"])
  med_clin =  median(dat$cindex[dat$type=="cinds_clin"])
  t = tidy(w)
  t$cancer = cancer
  t$med_lnc = med_lnc
  t$med_clin = med_clin
  imp = ((med_lnc/med_clin)-1) * 100
  t$imp = imp
  return(t)
}

wil = llply(cancers, check_perform)
wil = ldply(wil)
wil = as.data.table(wil)
colnames(wil)[2] = "pval"
wil = wil[order(pval, -med_lnc)]
wil_sig = filter(wil, pval <= 0.05)
wil_sig$imp = round(wil_sig$imp, digits=2)
wil_sig = as.data.table(wil_sig)
wil_sig = wil_sig[order(med_lnc)]

#remove those that have no significant difference between 3 groups
z = which(res$TYPE %in% wil_sig$cancer)
res = res[z,]

res$type = as.character(res$type)
res$TYPE = factor(res$TYPE, levels = wil_sig$cancer)

res$type[res$type == "cinds_justlncs"] = "lncRNA \nExpression"
res$type[res$type == "cinds_clin"] = "Clinical"
res$type[res$type == "cinds_combined"] = "lncRNA \nExpression & Clinical"

#label = % increase
labels = wil_sig$imp
colnames(wil_sig)[5] = "TYPE"
res = merge(res, wil_sig, by="TYPE")

wil_sig = as.data.table(wil_sig)
wil_sig = wil_sig[order(med_lnc)]
res$TYPE = factor(res$TYPE, levels = wil_sig$TYPE)
z = which(duplicated(res$imp))
res$imp[z] = ""


#-----------final plot figure 2D-------------------------------------------------------------------------------------- 

pdf("cindices_results_elastic_net_100CVs_No_FDR_cancers_July11_nofacet_horizontal_boxplots.pdf", width=7, height=7)
g = ggplot(res, aes(TYPE, cindex)) + 
geom_boxplot(aes(fill = type), color="black", outlier.size = 0.01) + theme_bw() + xlab("Cancer") + ylab("C-index") +
stat_compare_means(aes(group = type), label = "p.signif") + 
coord_flip() + scale_fill_manual(values=mypal) +
geom_hline(yintercept=0.5, linetype="dashed", color = "red")
ggpar(g, legend.title= "Predictor \nVariables") +
 theme(text = element_text(size=10),
        axis.text.x = element_text(hjust=1, size=10))
dev.off()

#----------------------------------------------------------------------------------------------------------------------



#pdf("cindices_results_elastic_net_100CVs_No_FDR_cancers_May2nd_nofacet.pdf", width=5, height=5)

  for(i in 1:length(unique(res$canc))){
      newres = subset(res, canc %in% unique(res$canc)[i])

    g = ggboxplot(newres, x="type", y="cindex", color="type", order=c("cinds_clin", "cinds_justlncs", "cinds_combined"), title=newres$canc[1]) + 
    coord_cartesian(ylim = c(0, 1))

    g =  g + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "cinds_clin") 

    print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red") + 

    stat_summary(fun.y=median, geom="line", aes(group=1))  + 
    stat_summary(fun.y=median, geom="point")+

    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
    theme(strip.text.x = element_text(size = 6.5, colour = "black")))
    
    }
dev.off()

#----------------------------------------------------------------------------------------------------------------------
### Reading in files with genes that were selected via feature selection
#----------------------------------------------------------------------------------------------------------------------

files = list.files(pattern = "lncRNAs_selected_1000CV_1000_no_fdr_ELASTICNET.rds")

###change to be universal 

genes_keep = as.data.frame(matrix(ncol = 4)) ; colnames(genes_keep) = c("Geneid", "NumtimesChosen" , "Cancer", "GeneName")
all_genes = as.data.frame(matrix(ncol = 4)) ; colnames(genes_keep) = c("Geneid", "NumtimesChosen" , "Cancer", "GeneName")


#DO NOT RUN
#pdf("summary_of_genes_results_elastic_net_100CVs_No_FDR_cancers_May2nd.pdf", width= 7, height=5)

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

#DO NOT RUN
	
#print(gg +  theme(text = element_text(size=10),
#        axis.text.x = element_text(angle=90, hjust=1, size=4)) +
#          geom_hline(yintercept = 50, linetype = 2, color = "blue"))

#genes_list = dplyr::filter(genes_list, N >=50)
if(!(dim(genes_list)[1]==0)){
colnames(genes_list) = colnames(genes_keep) 
genes_keep = rbind(genes_keep, genes_list)
}
}

#dev.off()

genes_keep = genes_keep[-1,]
genes_keep = as.data.table(genes_keep)
all_genes = genes_keep
write.csv(all_genes, file="all_genes_selected_by_elastic_net_final_run_aug7.csv", quote=F, row.names=F)
genes_keep = filter(genes_keep, NumtimesChosen >=50)

dups = unique(genes_keep$GeneName[which(duplicated(genes_keep$GeneName))])
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

pdf("summary_of_final_lncRNA_candidates_elastic_net_100CVs_No_FDR_cancers_May2nd.pdf", width=12, height=8)

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
          ylab = "# of lncRNAs selected more than 50%",  label = TRUE, label.pos = "out")           # Rotate vertically x axis texts

overg = overg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=55, hjust=1, size=8))
         

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
#gg = g + theme(legend.position = "bottom", legend.text = element_text(colour="black", size=7))
gg = ggpar(g, legend="none")
gg = gg +  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1, size=4)) +
          geom_hline(yintercept = 50, linetype = 2, color = "blue") + ggtitle("Summary of genes selected 100 CVs with No FDR")

#combine the two side by side
gg +overg + plot_layout(ncol = 1, heights = c(3, 1))

dev.off()

saveRDS(genes_keep, file="genes_keep_100CV_No_FDR_May2nd2018.rds")

colnames(t) = c("TYPE", "Nun_lncRNAs")
cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "TYPE"
summary = merge(t, cancers_conv, by="TYPE")
summary = summary[order(Nun_lncRNAs)]

pdf("barplot_summary_lncs_selected_100CVs_elastic_net_july13.pdf", width=8, height=5)
overg = ggbarplot(summary, x = "type", y = "Nun_lncRNAs",
          fill = "grey",               # change fill color by cyl
          palette = "mypal2", 
          color = "white",            # Set bar border colors to white
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          ggtheme = theme_bw(), 
          xlab = "Cancers",
          legend = "none",  
          x.text.angle = 45 , 
          ylab = "# of lncRNAs selected more than 50%",  label = TRUE, label.pos = "out")           # Rotate vertically x axis texts

#overg = overg +  theme(text = element_text(size=10),
#        axis.text.x = element_text(angle=55, hjust=1, size=8))
overg
dev.off()




