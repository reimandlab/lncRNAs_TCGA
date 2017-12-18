###---------------------------------------------------------------
###TCGA_cancers_data5.R
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

compiled_data = readRDS("compiled_data_lncRNAs_highlow_gtex.rds")
clin_exp = readRDS("lnc_rna_ovary_liver_plus_clinical.rds")
clin_exp$patient = rownames(clin_exp)
clin_exp = clin_exp[,c(5920:5924)]
clin_exp = subset(clin_exp, canc == "Liver hepatocellular carcinoma")
clin_exp$time = as.numeric(clin_exp$time)
z = which(is.na(clin_exp$time))
clin_exp = clin_exp[-z,]

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
fantom = as.data.table(fantom)

#filtered data
filtered_data = readRDS("filtered_data_with_TAGS.rds")
filtered_data = filter(filtered_data, canc == "liver_tcga")
filtered_data = merge(filtered_data, clin_exp, by = "patient")
colnames(filtered_data)[2] = "gene"

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

#1. Get list of genes 

genes_to_test = function(dat){
	genes = unique(dat$gene)
	genes = filter(fantom, CAT_geneID %in% genes)
	colnames(genes)[1] = "gene"
	genes = merge(genes, filtered_data, by = "gene")
	return(genes)
}

genes_test = llply(compiled_data, genes_to_test)

#2. Calculate survival 

calc = function(genee){
		
		results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "type_predictor")
		dat = as.data.table(dat)
		dat = filter(dat, gene %in% genee)
		dat$status[dat$status=="Alive"] <- 0
        dat$status[dat$status=="Dead"] <- 1
        dat$status <- as.numeric(dat$status)
        dat$time <- as.numeric(dat$time)
        dat$order = log1p(dat$order)
		#model 1 - using logged expression values  
		model1 = coxph(Surv(time, status) ~ order, data = dat)
		row <- c(genee, summary(model1)$coefficients[1,c(1,2,5)], "logged_tpm")
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)

		#model 2 - using mean cutoff
		model2 = coxph(Surv(time, status) ~ mean, data = dat)
		row <- c(genee, summary(model2)$coefficients[1,c(1,2,5)], "mean_cutoff")
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)

		#model 3 - using median cutoff 
		model3 = coxph(Surv(time, status) ~ median, data = dat)
		row <- c(genee, summary(model3)$coefficients[1,c(1,2,5)], "median_cutoff")
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)

		#model 4 - using quantile 1 cutoff 
		model4 = coxph(Surv(time, status) ~ quantile1, data = dat)
		row <- c(genee, summary(model4)$coefficients[1,c(1,2,5)], "quantile1_cutoff")
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)

		#model 5 - using quantile 3 cutoff
		model5 = coxph(Surv(time, status) ~ quantile3, data = dat)
		row <- c(genee, summary(model5)$coefficients[1,c(1,2,5)], "quantile3_cutoff")
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)
        results_cox = results_cox[-1,]

        return(results_cox)
	}


calc_survival = function(dat){
	genes = as.list(unique(dat$gene))
	results = llply(genes, calc, .progress = "text")	
	results = ldply(results, data.frame)
	results$fdr = p.adjust(results$pval, method="fdr")
	results = as.data.table(results)
	results = filter(results, fdr <=0.1)
	return(results)
}

survival_results = as.data.frame(matrix(ncol=6)) ; colnames(survival_results) <- c("gene", "coef", "HR", "pval", "type_predictor", "fdr")

colnames(survival_results) = colnames(genes_test[[1]])
for(i in 1:length(genes_test)){
	dat = genes_test[[i]]
	new = calc_survival(dat)
	survival_results = rbind(survival_results, new)
}

survival_results = survival_results[-1,]
compiled_data = ldply(compiled_data, data.frame)
compiled_data = as.data.table(compiled_data)

genes1 = filter(compiled_data, type_test == "high_tum_great")$gene
genes2 = filter(compiled_data, type_test == "high_tum_lower")$gene
genes3 = filter(compiled_data, type_test == "low_tum_great")$gene
genes4 = filter(compiled_data, type_test == "low_tum_lower")$gene

survival_results$group1 = ""
survival_results$group2 = ""
survival_results$group3 = ""
survival_results$group4 = ""

survival_results$group1[which(survival_results$gene %in% genes1)] = "high_tum_great"
survival_results$group2[which(survival_results$gene %in% genes2)] = "high_tum_lower"
survival_results$group3[which(survival_results$gene %in% genes3)] = "low_tum_great"
survival_results$group4[which(survival_results$gene %in% genes4)] = "low_tum_lower"

#4 lists of cox results, adjsut pvalues within each one

survival_results = as.data.table(survival_results)
survival_results = survival_results[order(fdr)]

z <- which(duplicated(survival_results$gene))
dups = survival_results$gene[z]
dups = filter(survival_results, gene %in% dups)

unique_sig_genes = unique(survival_results$gene)

##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("ALL_including_justtpm_survival_results_Dec18_tcga_liver.pdf", pointsize=6, width=9, height=8)
require(gridExtra)

filtered_data = as.data.table(filtered_data)

for(i in 1:nrow(survival_results)){
 	
	genex = survival_results$gene[i]
	type = survival_results$type_predictor[i]
	#if(!(type == "logged_tpm")){
		df = filter(filtered_data, gene == genex)

 	 	#cox
        df$status[df$status=="Alive"] <- 0
        df$status[df$status=="Dead"] <- 1
        df$status <- as.numeric(df$status)
        df$time <- as.numeric(df$time)
      	
      	genexname = fantom$CAT_geneName[which(fantom$CAT_geneID %in% genex)]

          #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = df)
          s <- ggsurvplot(
          title = paste(genexname, "Survival Curve Split by Median"),
          fit, 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          main = paste(gene, df$canc[1]),       
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = df,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,2000),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 500,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)  
      #}
}


dev.off()

saveRDS(survival_results, file="survival_results_DEC15_liver_TCGA.rds")





















