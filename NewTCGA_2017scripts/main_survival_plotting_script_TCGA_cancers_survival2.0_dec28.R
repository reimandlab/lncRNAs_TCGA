###---------------------------------------------------------------
###TCGA_cancers_survival1.1_dec20.R
###---------------------------------------------------------------

###December 20th, 2017
###establish a list of lncRNAs to evaluate using survival analysis
###by looking at mean versus variance density for each lncRNA 
###within each cancer type 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. lncRNA expression in different cancers 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

lncrnas = colnames(rna)[1:(ncol(rna)-5)]

rna$patient = rownames(rna) ; 

pcg = readRDS("19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")
pcgs = colnames(pcg)[1:(ncol(pcg)-5)]
pcg$patient = rownames(pcg)
rna = merge(rna, pcg, by = c("canc", "time", "status", "sex", "patient"))

#matched normal 
norm = readRDS("5919_lncs4matched_normal_tissues_TCGAnew.rds")

#3. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

#4. Fantom data 
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
colnames(fantom)[1] = "gene"

#5. List of lncRNA survival associated candidates (from PCAWG)
#cands = read.csv("42_original_candidates.csv", header=F)
#colnames(cands) = c("lncRNA", "cancer")

#6. Cands from monte carlo CV 
cands = data.frame(gene=c("ENSG00000232018" ,"ENSG00000237149", "ENSG00000153363", "ENSG00000234967",
 "ENSG00000236393", "ENSG00000258082", "ENSG00000248210", "ENSG00000249662",
 "ENSG00000272243", "ENSG00000261189", "ENSG00000236543", "ENSG00000234634",
 "ENSG00000215483", "ENSG00000231265", "ENSG00000259410", "ENSG00000237686",
 "ENSG00000235572", "ENSG00000249036", "ENSG00000165655"), cancer=rep("Ovarian serous cystadenocarcinoma",19))

###---------------------------------------------------------------
###Get lncRNAs for each cancer 
###---------------------------------------------------------------

#1. Remove patients with unlabelled cancer type 
#z <- which(rna$canc == "")
#rna = rna[-z,]
#2. Remove patients without survival status 
#z <- which(rna$status == "")
#rna = rna[-z,]

#for(i in 1:nrow(rna)){
	#canc = unlist(strsplit(rna$canc[i], " "))[[1]]
	#if(canc == "Ovarian"){
		#canc = "Ovary"
	#}
	#rna$canc[i] = canc
#}

rna = subset(rna, canc %in% cands$cancer)

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")
cands = subset(cands, cancer %in% rna$canc)

for(i in 1:nrow(cands)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% cands$cancer[i])
  z <- which(colnames(df) %in% cands$gene[i])
  
  if(!(length(z)==0)){
  df <- df[,c(z,1:5)]  

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
  if(median2 ==0){
  	median2 = mean(as.numeric(df[,1]))
  }

  #median2 <- median(df[,1])
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 

  gene <- colnames(df)[1]
  #cox
        df$status[df$status=="Alive"] <- 0
        df$status[df$status=="Dead"] <- 1
        df$status <- as.numeric(df$status)
        df$time <- as.numeric(df$time)
      
        #cox regression 
        res.cox <- coxph(Surv(time, status) ~ median, data = df)
        row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)], df$canc[1])
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)  
}
}

results_cox <- results_cox[-1,]
results_cox$pval <- as.numeric(results_cox$pval)

results_cox <- as.data.table(results_cox)
results_cox <- results_cox[order(pval)]
#results_cox = subset(results_cox, pval <=0.05)
results_cox$name = ""
#for(i in 1:nrow(results_cox)){
for(i in 1:18){
	results_cox$name[i] = fantom$CAT_geneName[which(fantom$gene %in% results_cox$gene[i])]
}
results_cox$name[19] = "ZNF503"

##+++++++++++++++++++++++++++++
##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("ZNF503inc_results_lncRNAs_from_MCCV_10rounds.pdf", pointsize=6, width=10, height=8)
require(gridExtra)

results_cox = as.data.frame(results_cox)

for(i in 1:nrow(results_cox)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% results_cox$canc[i])
  #gene = fantom$gene[which(fantom$CAT_geneName %in% results_cox$name[i])]
  z <- which(colnames(df) %in% results_cox$gene[i])
  if(!(length(z)==0)){
  df <- df[,c(z,1:5)]  

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
  if(median2 ==0){
  	median2 = mean(as.numeric(df[,1]))
  }
  #median2 <- median(df[,1])
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 

  gene <- results_cox$name[i]
  #plot boxplot showing difference between the two groups and sex
  title <- paste(gene, df$canc[1], "Expression")
  colnames(df)[1] <- "Gene"
  g <- ggboxplot(df, x= "median", y="Gene", palette=mypal[c(4,1)], order=c("0", "1"), fill = "median",  add = "jitter")
  g <- g + stat_compare_means()
  g <- ggpar(g, font.legend = c(8, "plain", "black")) 
  g <- g + labs(title = title, y="log1p(FPKM)", x="Median") + 
      theme(plot.title = element_text(hjust = 0.5))
      print(g)

  #cox
        df$status[df$status=="Alive"] <- 0
        df$status[df$status=="Dead"] <- 1
        df$status <- as.numeric(df$status)
        df$time <- as.numeric(df$time)
        df$time = df$time/365
      
          #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = df)
          s <- ggsurvplot(
          title = paste(gene, df$canc[1]),
          fit, 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = df,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)  

}
}


dev.off()

###ZNF503-AS2

rna = readRDS("5919_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")

lncrnas = colnames(rna)[1:(ncol(rna)-5)]

rna$patient = rownames(rna) ; 

z = which(colnames(rna) == "ENSG00000237149")
rna = rna[,c(z, 5920:ncol(rna))]

pcg = readRDS("19438_lncRNAs_tcga_all_cancers_Dec20_wclinical_data.rds")
pcgs = colnames(pcg)[1:(ncol(pcg)-5)]
pcg$patient = rownames(pcg)

#ZNF503 
z = which(colnames(pcg) == "ENSG00000165655")
pcg = pcg[,c(z, 19439:ncol(pcg))]

znf = merge(rna, pcg, by = c("canc", "time", "status", "sex", "patient"))

#additional clinical data (stage, grade)
clin2 = fread("OV_clinical_core.txt", data.table=F)

#Check how they are correlated 

pdf("ZNF503_ZNF503AS2_correlation.pdf")

ggscatter(znf, x = "ENSG00000237149", y = "ENSG00000165655",
   title= "Pearson Correlation",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.sep = "\n")
   )

ggscatter(znf, x = "ENSG00000237149", y = "ENSG00000165655",
   title= "Spearman Correlation",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )


znf[,6:7] = log1p(znf[,6:7])

ggscatter(znf, x = "ENSG00000237149", y = "ENSG00000165655",
   title= "Pearson Correlation",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.sep = "\n")
   )

ggscatter(znf, x = "ENSG00000237149", y = "ENSG00000165655",
   title= "Spearman Correlation",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )

dev.off()











