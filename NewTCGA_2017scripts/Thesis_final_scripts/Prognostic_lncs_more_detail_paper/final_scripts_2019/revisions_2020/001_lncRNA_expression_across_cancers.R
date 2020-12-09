source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
library(ggridges)

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

dim(rna)
rownames(rna) = rna$patient

mypal = c("#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% c("type", "patient"))
rna = as.data.frame(rna)
rna = rna[,c(z1, z2)]
z1 = which(str_detect(colnames(rna), "ENSG"))
rna[,z1] = log1p(rna[,z1])

#get summary of gene expression across each cancer type
rna = as.data.table(melt(rna))

sums = as.data.table(rna %>% group_by(type, variable) %>% dplyr::summarize(sum = sum(value)))
sums = filter(sums, sum == 0)
sums$combo = paste(sums$type, sums$variable)
rna$combo =  paste(rna$type, rna$variable)
rna = as.data.table(filter(rna, !(combo %in% sums$combo)))

medians= as.data.table(rna %>% group_by(type) %>% dplyr::summarize(median = median(value)))
medians=medians[order(-median)]
rna$type = factor(rna$type, levels=medians$type)

theme_set(theme_minimal())

pdf("/u/kisaev/lncRNA_expression_across_cancers.pdf", height=10)
g = ggviolin(rna, x="type", y="value", fill="type",
legend="none", xlab="Cancer", ylab="log1p(FPKM-UQ)",
 palette=mypal)
g=ggpar(g, xtickslab.rt = 90)+theme_bw()

ggplot(rna, aes(x = value, y = type)) +
  geom_density_ridges(aes(fill = type)) +
  scale_fill_manual(values = mypal) +
  theme(legend.position = "none")+
  ylab("Cancer") + xlab("log1p(FPKM-UQ)")

dev.off()

#just KIRC
kirc = as.data.table(filter(kirc, type == "KIRC"))
