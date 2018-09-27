library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(gridExtra)
library(grid)
library(Hmisc)
library(broom)
library(EnvStats)
library(plyr)
library(dplyr)

#Data---------------------------------------------------

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion<- fread("pcawgConversion.tsv", data.table=F)

#RNA data 
pcg_rna = readRDS(file="all_mrna_may8th.rds")

colnames(patient_table) = c("patient", "cancer")
pcg_rna = merge(pcg_rna, patient_table, by = "patient")

lncrna = readRDS(file="all_lncRNA_may8th.rds")
pcg_rna = merge(pcg_rna, lncrna, by = colnames(pcg_rna)[which(colnames(pcg_rna) %in% colnames(lncrna))])

#for each cancer type look at correlation between the two genes 
cancers = as.list(unique(pcg_rna$cancer))

#ucsc gene ids
ucsc = fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt")

#**** translocations
load("k.frame_complete_25000.rsav")
head(k.frame_complete_25000)
k.frame_complete = k.frame_complete_25000
k.frame_complete$gene_id = unlist(lapply(k.frame_complete$gene_id, function(x){unlist(strsplit(x, "::"))[1]}))
k.frame_complete$FMRE_id = unlist(lapply(k.frame_complete$FMRE_id, function(x){unlist(strsplit(x, "::"))[2]}))
k.frame_complete$fmre_start = unlist(lapply(k.frame_complete$FMRE_id, function(x){unlist(strsplit(x, ":"))[2]}))
k.frame_complete$fmre_start = unlist(lapply(k.frame_complete$fmre_start, function(x){unlist(strsplit(x, "-"))[1]}))
k.frame_complete$fmre_flank = as.numeric(k.frame_complete$fmre_start) - as.numeric(k.frame_complete$FMRE_vicinity_start)

#------FUNCTION---------------------------------------------------------------------------------

pcg_rna = as.data.frame(pcg_rna)

#subset to genes 
ens = unique(ucsc$hg19.ensGene.name2[ucsc$hg19.ensemblToGeneName.value %in% k.frame_complete$gene_id])
library(stringr)
z = which(str_detect(colnames(pcg_rna), "ENSG"))
clin = colnames(pcg_rna)[-z]

pcg_rna = pcg_rna[,which(colnames(pcg_rna) %in% c(clin, ens))]
pcg_rna = as.data.table(pcg_rna)

check_exp = function(row) {

	#get cancer type of patient 
	pat = row[[11]] #patient id 
	canc_pat = pcg_rna$cancer[which(pcg_rna$patient == pat)]

	if(!(length(canc_pat) ==0)){
	#compare to just that cancer
	canc_exp = (filter(pcg_rna, cancer == canc_pat))

	#subset to PCG
	pcg = row[[5]]

	#get ensembl id
	ens = ucsc$hg19.ensGene.name2[ucsc$hg19.ensemblToGeneName.value == pcg][1]
	z = which(colnames(pcg_rna) == ens)
	if(!(length(z)==0)){

	#expression of patient
	p = which(canc_exp$patient == pat)
	
	if(!(length(p) ==0)){
		pat_exp = canc_exp[which(canc_exp$patient == pat),z]

	#expression of everyone else
	else_ppl = mean(canc_exp[which(!(canc_exp$patient ==pat)),z])

	#fold change 
	fc = (pat_exp/else_ppl)
	#fc = log2(pat_exp/else_ppl)	

	#do wilcoxon test
	w = wilcox.test(canc_exp[which(canc_exp$patient == pat),z], 
	canc_exp[which(!(canc_exp$patient ==pat)),z])
	print(tidy(w))
	w = tidy(w)
	w$canc = canc_pat
	w$pat = pat
	w$pcg = pcg
	w$flank = row[[14]]
	w$fc = fc
	print(pcg)
    return(w)
}
}
}
}

exp_results = apply(k.frame_complete, 1, check_exp)
exp_results = ldply(exp_results) #50/83 genes were able to be evaluated 
exp_results = as.data.table(exp_results)
exp_results = exp_results[order(p.value)]


pdf("translocations_fmres_expression_plots.pdf")

for(i in 1:nrow(exp_results)){

	#get cancer type of patient 
	pat = exp_results$pat[i]
	canc_pat = exp_results$canc[i]
	pval = round(exp_results$p.value[i], digits=4)

	#compare to just that cancer
	canc_exp = subset(pcg_rna, cancer == canc_pat)

	#subset to PCG
	pcg = exp_results$pcg[i]
	print(pcg)
	#get ensembl id
	ens = ucsc$hg19.ensGene.name2[ucsc$hg19.ensemblToGeneName.value == pcg][1]
	z = which(colnames(pcg_rna) == ens)

	#expression of patient
	p = which(canc_exp$patient == pat)
	canc_exp = as.data.frame(canc_exp)
	#plot boxplot
	plot = canc_exp[, which(colnames(canc_exp) %in% c(ens, "cancer", "patient"))]
	plot$tran = ""
	plot$tran[which(plot$patient == pat)] = "Yes"
	plot$tran[!(plot$patient == pat)] = "No"
	colnames(plot)[2] = "GeneExpression"

	g = ggboxplot(plot, x="tran", y="GeneExpression", fill="tran", palette=c("lightcyan4", "darksalmon")) + 
	ggtitle(paste(canc_pat, pcg, "U-test p-val = ", pval)) + 
	stat_n_text(size=5) + theme_bw() 
    print(g)

}

dev.off()

write.table(exp_results, file="233_translocations_wexp.txt", quote=F, row.names=F, sep="\t")
saveRDS(exp_results, file="233_translocations_wexp.rds")





