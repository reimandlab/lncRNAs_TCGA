###---------------------------------------------------------------
###TCGA_cancers_data4.R
###---------------------------------------------------------------

###Compare relative expression differences between liver
###cancer and gtex liver samples

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Data
###---------------------------------------------------------------

liver_tum = readRDS("TCGA_liver_ranked_genes_Dec14.rds")
liver_norm = readRDS("GTEX_liver_ranked_genes_Dec14.rds")

#subset liver_tum to only genes in liver_norm 
liver_tum = as.data.table(liver_tum)
liver_tum = filter(liver_tum, gene %in% liver_norm$Name)
colnames(liver_tum)[1] = "Name"
liver_tum = as.data.frame(liver_tum)

#Fantom data 
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

###---------------------------------------------------------------
###Wilcoxon test compare difference 
###---------------------------------------------------------------

#1. Only care about lncRNAs for this 
all_data = rbind(liver_tum, liver_norm[,c(1,4:9)])
all_data = as.data.table(all_data)
all_data = filter(all_data, type == "lncRNA")

#2. For each gene, want to compare the difference in expression 
#between normal and tumours 

genes = as.list(unique(all_data$Name))

#3. Seperate TCGA into low and high, 4 different ways and compare 
#to normal tissues 

add_exp_tag = function(gene){
	data = subset(all_data, all_data$Name == gene)
	#subset to tcga 
	data_t = subset(data, canc == "liver_tcga")
	mean = mean(data_t$order)
	median = median(data_t$order)
	quantile1 = quantile(data_t$order)[2]  
	quantile3 = quantile(data_t$order)[4]  
	
	data_t$mean[data_t$order > mean] = "High"
	data_t$mean[data_t$order <= mean] = "Low"

	data_t$median[data_t$order > median] = "High"
	data_t$median[data_t$order <= median] = "Low"

	data_t$quantile1[data_t$order > quantile1] = "High"
	data_t$quantile1[data_t$order <= quantile1] = "Low"

	data_t$quantile3[data_t$order > quantile3] = "High"
	data_t$quantile3[data_t$order <= quantile3] = "Low"

	#if not possible to divide, ie everyone has super low expression but one person
	#has expression slightly greater than 1

	quant = 0.25*nrow(data_t)
	quant = floor(quant)
	meantest = (table(data_t$mean)[1] >= quant) & (table(data_t$mean)[2] >= quant) 
	mediantest = (table(data_t$median)[1] >= quant) & (table(data_t$median)[2] >= quant) 
	quantile1test = (table(data_t$quantile1)[1] >= quant) & (table(data_t$quantile1)[2] >= quant) 
	quantile3test = (table(data_t$quantile3)[1] >= quant) & (table(data_t$quantile3)[2] >= quant) 

	if(meantest & mediantest & quantile1test & quantile3test){
		return(gene)
	}
}

filters = llply(genes, add_exp_tag, .progress = "text")
filters = Filter(Negate(is.null), filters) #2930 remain out of 5376, need to go back to them and check why they were removed
#most likely because they had very low read TPMs 

#4. For the genes that are seperable into low and high groups 
#add tag 

add_exp_tag_to_all_data = function(gene){
	data = subset(all_data, all_data$Name == gene)
	#subset to tcga 
	data_t = subset(data, canc == "liver_tcga")
	mean = mean(data_t$order)
	median = median(data_t$order)
	quantile1 = quantile(data_t$order)[2]  
	quantile3 = quantile(data_t$order)[4]  
	
	data_t$mean[data_t$order > mean] = "High"
	data_t$mean[data_t$order <= mean] = "Low"

	data_t$median[data_t$order > median] = "High"
	data_t$median[data_t$order <= median] = "Low"

	data_t$quantile1[data_t$order > quantile1] = "High"
	data_t$quantile1[data_t$order <= quantile1] = "Low"

	data_t$quantile3[data_t$order > quantile3] = "High"
	data_t$quantile3[data_t$order <= quantile3] = "Low"

	#gtex 
	data_g = subset(data, canc == "liver_gtex")
	data_g$mean = "GTEx"
	data_g$median  = "GTEx"
	data_g$quantile1 = "GTEx"
	data_g$quantile3 = "GTEx"

	#combine
	data = rbind(data_t, data_g)
	return(data)
}

filtered_data = llply(filters, add_exp_tag_to_all_data, .progress = "text")

#5. Wilcoxon test

get_wilcoxon = function(data){
	
	#1. Compare wilcoxon high group vesus gtex
	high_mean = data$score[which(data$mean == "High")]
	pval_mean = wilcox.test(high_mean, data$score[data$canc == "liver_gtex"])$p.value

	high_median = data$score[which(data$median == "High")]
	pval_median = wilcox.test(high_median, data$score[data$canc == "liver_gtex"])$p.value

	high_quantile1 = data$score[which(data$quantile1 == "High")]
	pval_quantile1 = wilcox.test(high_quantile1, data$score[data$canc == "liver_gtex"])$p.value

	high_quantile3 = data$score[which(data$quantile3 == "High")]
	pval_quantile3 = wilcox.test(high_quantile3, data$score[data$canc == "liver_gtex"])$p.value

	test_high = (pval_mean & pval_median & pval_quantile1 & pval_quantile3)
	test_high_pval = mean(pval_mean, pval_median, pval_quantile1, pval_quantile3) <= 0.01
	high_mean_pval = mean(pval_mean, pval_median, pval_quantile1, pval_quantile3)

	#2. Compare wilcoxon low group versus gtex
	low_mean = data$score[which(data$mean == "Low")]
	pval_mean = wilcox.test(low_mean, data$score[data$canc == "liver_gtex"])$p.value

	low_median = data$score[which(data$median == "Low")]
	pval_median = wilcox.test(low_median, data$score[data$canc == "liver_gtex"])$p.value

	low_quantile1 = data$score[which(data$quantile1 == "Low")]
	pval_quantile1 = wilcox.test(low_quantile1, data$score[data$canc == "liver_gtex"])$p.value

	low_quantile3 = data$score[which(data$quantile3 == "Low")]
	pval_quantile3 = wilcox.test(low_quantile3, data$score[data$canc == "liver_gtex"])$p.value

	test_low = (pval_mean & pval_median & pval_quantile1 & pval_quantile3)
	test_low_pval = mean(pval_mean, pval_median, pval_quantile1, pval_quantile3) <= 0.01
	low_mean_pval = mean(pval_mean, pval_median, pval_quantile1, pval_quantile3)

	#If at least one of them different plot boxplot and save data
	final_test = test_high_pval | test_low_pval

	if(final_test){
		return(data)
	}
}

sig_data = llply(filtered_data, get_wilcoxon,  .progress = "text")
# 2930/2930 are different than normal tissues 
filtered_data = ldply(filtered_data, data.frame)
filtered_data = as.data.table(filtered_data)

#6. Get pvalue and change in means between each type of test 

get_wilcoxon_data = function(data){
	#1. Compare wilcoxon high group vesus gtex
	high_mean = data$score[which(data$mean == "High")]
	pval_mean = wilcox.test(high_mean, data$score[data$canc == "liver_gtex"])$p.value
	change_mean = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(high_mean))

	high_median = data$score[which(data$median == "High")]
	pval_median = wilcox.test(high_median, data$score[data$canc == "liver_gtex"])$p.value
	change_median = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(high_median))

	high_quantile1 = data$score[which(data$quantile1 == "High")]
	pval_quantile1 = wilcox.test(high_quantile1, data$score[data$canc == "liver_gtex"])$p.value
	change_quantile1 = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(high_quantile1))

	high_quantile3 = data$score[which(data$quantile3 == "High")]
	pval_quantile3 = wilcox.test(high_quantile3, data$score[data$canc == "liver_gtex"])$p.value
	change_quantile3 = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(high_quantile3))

	high_data = c("high", change_mean, pval_mean, change_median, pval_median, change_quantile1, pval_quantile1, change_quantile3, pval_quantile3)
	names(high_data) = c("group", "change_mean", "pval_mean", "change_median", "pval_median", "change_quantile1", "pval_quantile1", "change_quantile3", "pval_quantile3")

	#2. Compare wilcoxon low group versus gtex
	low_mean = data$score[which(data$mean == "Low")]
	pval_mean = wilcox.test(low_mean, data$score[data$canc == "liver_gtex"])$p.value
	change_mean = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(low_mean))

	low_median = data$score[which(data$median == "Low")]
	pval_median = wilcox.test(low_median, data$score[data$canc == "liver_gtex"])$p.value
	change_median = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(low_median))

	low_quantile1 = data$score[which(data$quantile1 == "Low")]
	pval_quantile1 = wilcox.test(low_quantile1, data$score[data$canc == "liver_gtex"])$p.value
	change_quantile1= log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(low_quantile1))

	low_quantile3 = data$score[which(data$quantile3 == "Low")]
	pval_quantile3 = wilcox.test(low_quantile3, data$score[data$canc == "liver_gtex"])$p.value
	change_quantile3  = log2(mean(data$score[data$canc == "liver_gtex"])) - log2(mean(low_quantile3))

	low_data = c("low", change_mean, pval_mean, change_median, pval_median, change_quantile1, pval_quantile1, change_quantile3, pval_quantile3)
	names(low_data) = c("group", "change_mean", "pval_mean", "change_median", "pval_median", "change_quantile1", "pval_quantile1", "change_quantile3", "pval_quantile3")

	results = rbind(low_data, high_data)
	results = as.data.frame(results)
	results$gene = data$Name[1]

	return(results)
	
}

sig_data = llply(sig_data, get_wilcoxon_data,  .progress = "text")
sig_data = Filter(Negate(is.null), sig_data)
sig_data = as.data.frame(do.call(rbind, sig_data))

#7. FDR adjustment and save the genes that have all fdr sig pvalues 
sig_data$fdr_mean_pval = p.adjust(sig_data$pval_mean, method="fdr")
sig_data$fdr_median_pval = p.adjust(sig_data$pval_median, method="fdr")
sig_data$fdr_quantile1_pval = p.adjust(sig_data$pval_quantile1, method="fdr")
sig_data$fdr_quantile3_pval = p.adjust(sig_data$pval_quantile3, method="fdr")

sig_data = as.data.table(sig_data)
sig_data = filter(sig_data, fdr_mean_pval <=0.01, fdr_median_pval <=0.01, fdr_quantile1_pval <=0.01, fdr_quantile3_pval <=0.01)
#all still significant 
sig_data = as.data.table(sig_data)

#8. Use median for now and plot boxplots showing differences between expression of low, high and normal 
genes = unique(sig_data$gene)

make_boxplots = function(g){
	data = filter(sig_data, gene == g)
	data = filter(filtered_data, Name == g)
	my_comparisons <- list( c("GTEx", "Low"), c("Low", "High"), c("GTEx", "High") )
	name = fantom$CAT_geneName[which(fantom$CAT_geneID %in% g)]
	f = ggboxplot(data, title=paste(name, "Expression"), x="quantile3", y="score", color="quantile3", fill="quantile3", palette=mypal[c(3,4,1)], ggtheme=theme_bw(), order=c("GTEx", "Low", "High"), add="jitter")
	f = ggpar(f, ylab="Score", x.text.angle=65, font.tickslab=c(14, "plain", "black"), legend="right", 
		font.x = c(18, "plain", "black"),
   		font.y = c(18, "plain", "black"))
	#f = f + rremove("y.grid") 
	f = f + stat_compare_means(comparisons = my_comparisons, label.y = c(0.9, 0.85, 0.8))
	return(f)
}

pdf("new_QUANTIL3_Diff_exp_genes_normal_liver_tumour.pdf")
llply(genes, make_boxplots, .progress = "text")
dev.off()

#Order by greatest differences in means 
sig_data = sig_data[order(change_median)]

#1. Divide into high tumour greater than gtex
high_tum_great = as.data.table(filter(sig_data, group == "high", change_median < 0))
high_tum_great$type = "high_tum_great"

#low tumour greater than gtex
low_tum_great = as.data.table(filter(sig_data, group == "low", change_median < 0))
low_tum_great$type = "low_tum_great"

#high tumour lower than gtex
high_tum_lower = as.data.table(filter(sig_data, group == "high", change_median > 0))
high_tum_lower$type = "high_tum_lower"

#low tumour lower tha gtex 
low_tum_lower = as.data.table(filter(sig_data, group == "low", change_median > 0))
low_tum_lower$type = "low_tum_lower"

#save as a list of dataframes 
#within each one conduct survival analysis 

compiled_data = list(high_tum_great, low_tum_great, high_tum_lower, low_tum_lower)

saveRDS(compiled_data, file="compiled_data_lncRNAs_highlow_gtex.rds")

saveRDS(filtered_data, file="filtered_data_with_TAGS.rds")











