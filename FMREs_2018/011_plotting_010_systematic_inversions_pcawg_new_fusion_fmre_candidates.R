library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(Hmisc)
library(broom)
library(EnvStats)
library(plyr)
library(dplyr)
library(patchwork)


#Data---------------------------------------------------

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion<- fread("pcawgConversion.tsv", data.table=F)

#RNA data 
pcg_rna = readRDS(file="all_rna_may8th.rds")

colnames(patient_table) = c("patient", "cancer")
pcg_rna = merge(pcg_rna, patient_table, by = "patient")

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg_rna$cancer))

#ucsc gene ids
ucsc = fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt")

#------read in results-------------------------------------------------------------------------------

load("merged_translocation_frame.rsav")
res = as.data.table(k.frame)
res = as.data.table(filter(res, p.value < 0.1))

#function 1 - get fold change between mean of non-mutation group and translocation person

get_fc = function(row){

	pat = row[[2]]
	canc = row[[13]]
	#subset to PCG
	pcg = row[[1]]

	#get ensembl id
	ens = ucsc$hg19.ensGene.name2[ucsc$hg19.ensemblToGeneName.value == pcg][1]
	canc_exp = subset(pcg_rna, cancer == canc)
	z = which(colnames(canc_exp) == ens)
	#canc_exp[,z] = log1p(canc_exp[,z])

	#expression of patient
	pat_exp =canc_exp[which(canc_exp$patient == pat),z]

	#expression of everyone else
	else_ppl = mean(canc_exp[which(!(canc_exp$patient ==pat)),z])

	#fold change 
	fc = log2(pat_exp/else_ppl)

	#pvalue
	pval = row[[10]]
	#convert to star
	if(pval < 0.1){
		newpval = "."
	}
	if(pval < 0.05){
		newpval = "*"
	}

	fmre = row[[6]]

	dat = as.data.frame(matrix(ncol = 10))
	dat[1,] = c(pat, canc, pcg, ens, pat_exp, else_ppl, fc, pval, newpval, fmre)
	return(dat)

}

plotting_dat = apply(res, 1, get_fc)
plotting_dat = ldply(plotting_dat)
colnames(plotting_dat) = c("patient", "cancer", "name", "ens", "pat_exp", "mean_others", "fc", "pval", "newpval", "fmre")
plotting_dat$fc = as.numeric(plotting_dat$fc)
z = which(duplicated(plotting_dat$name))
plotting_dat$name[z] = paste(plotting_dat$name[z], ".")
plotting_dat = as.data.table(plotting_dat)
plotting_dat = plotting_dat[order(fc)]
plotting_dat$name = factor(plotting_dat$name, levels = plotting_dat$name)

#------barplot------------------------------------------------------------------------------------------

labels = plotting_dat$newpval

barplot = ggbarplot(plotting_dat, x="name", y="fc", lab.size = 6, color="fmre", label = labels, lab.vjust=-0.1) +
#scale_fill_gradient2(low='darkcyan', mid='snow3', high= "orange", space='Lab') +
 xlab("Translocation target") + ylab("Fold Change") + theme_bw() +
  scale_colour_brewer(palette="Dark2", breaks=plotting_dat$fmre) 

barplot = ggpar(barplot,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45, legend.title="FMRE")

#------fmre id covariate---------------------------------------------------------------------------------

mypal = readRDS("palette_23_cancer_types.rds")
mypal = mypal[c(1, 2, 6, 9, 21, 20, 8, 14)]

fmres = ggplot(plotting_dat, aes(name, 0.2)) +
    geom_tile(aes(fill = fmre)) + geom_text(aes(label = fmre), size=2) +
    theme_classic() + scale_fill_manual(values=mypal)

fmres = ggpar(fmres, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#------cancer covariate---------------------------------------------------------------------------------

mypal = fread("colors_cancer_types.txt", header=F)
colnames(mypal) = c("cancer", "color_name", "color_code")
plotting_dat = merge(plotting_dat, mypal, by="cancer")
plotting_dat = as.data.table(plotting_dat)
plotting_dat = plotting_dat[order(fc)]
plotting_dat$name = factor(plotting_dat$name, levels = plotting_dat$name)
mypal = filter(mypal, cancer %in% plotting_dat$cancer)
mypal$cancer = factor(mypal$cancer, levels = unique(plotting_dat$cancer))
plotting_dat$cancer = factor(plotting_dat$cancer, levels = unique(plotting_dat$cancer))

cancers = ggplot(plotting_dat, aes(name, 0.2)) +
    geom_tile(aes(fill = cancer)) + geom_text(aes(label = cancer), size=1.7) +
    theme_classic() + scale_fill_manual(values=unique(plotting_dat$color_code)) 

cancers = ggpar(cancers, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

pdf("summary_barplot_translocation_exp_fc.pdf", width=6, height=6)
barplot + cancers + plot_layout(ncol = 1, heights = c(10, 1))
dev.off()


#------barplot with x-axis ticks pcg + fmre--------------------------------------------------------------

plotting_dat$translocation = paste(plotting_dat$fmre, plotting_dat$name, sep=" - ")
labels = plotting_dat$newpval
plotting_dat$expression[plotting_dat$fc >0] = "Upregulated"
plotting_dat$expression[plotting_dat$fc <0] = "Downregulated"
plotting_dat$expression = factor(plotting_dat$expression, levels = c("Upregulated", "Downregulated"))

barplot = ggbarplot(plotting_dat, x="name", y="fc", lab.size = 6, color="fmre", label = labels, lab.vjust=-0.1, fill="fc") +
scale_fill_gradient2(low='darkcyan', mid='snow3', high= "orange", space='Lab')+
 xlab("Translocation target") + ylab("Fold Change") + theme_bw() +
 scale_colour_brewer(palette="Dark2", breaks=plotting_dat$fmre) 
 
barplot = ggpar(barplot,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

pdf("summary_barplot_translocation_exp_fc_version_2_.pdf", width=6, height=6)
barplot + cancers + plot_layout(ncol = 1, heights = c(10, 1))
dev.off()


#####3
barplot = ggbarplot(plotting_dat, x="translocation", y="fc", lab.size = 6, label = labels, lab.vjust=-0.1, fill="expression") +
 	xlab("Translocation") + ylab("log2(Fold Change)") + theme_bw() +
 	scale_fill_manual(values=c("#E69F00", "#56B4E9"))
 
barplot = ggpar(barplot,
 font.tickslab = c(7,"plain", "black"), legend.title = "Expression", font.legend = c(8, "plain", "black"), 
 xtickslab.rt = 45)

pdf("summary_barplot_translocation_exp_fc_version_3_.pdf", width=7, height=7)
barplot + cancers + plot_layout(ncol = 1, heights = c(10, 1))
dev.off()



####4

plotting_dat$labels = paste(plotting_dat$fmre, plotting_dat$newpval, sep=" ")
labels = plotting_dat$labels

barplot = ggbarplot(plotting_dat, x="name", y="fc", lab.size = 2, label = labels, lab.vjust=-0.3, fill="fc") +
scale_fill_gradient2(low='darkcyan', mid='snow3', high= "orange", space='Lab')+
 xlab("Translocation target") + ylab("log2(Fold Change)") + theme_bw() 
 barplot = ggpar(barplot, legend.title="log2(FC)",
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

pdf("summary_barplot_translocation_exp_fc_version_4_.pdf", width=6, height=6)
barplot + cancers + plot_layout(ncol = 1, heights = c(10, 1))
dev.off()


