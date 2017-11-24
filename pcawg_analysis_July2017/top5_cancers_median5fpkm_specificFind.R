#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 11th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#last updated August 28th

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library("colorout")
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#fantom 
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

#Clinical file 
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 

#---------------------------------------------------------
#Find cancer sepcific lncRNAs with median E >= 5 FPKM 
#---------------------------------------------------------

#Write function that takes a dataframe, calculates medians and 
#output list of genes with median greater than that 

check_medians <- function(column){
  med <- median(column)
  if(med >=5){
    return(med)
  } 
}

#save results
high_lncs <- as.data.frame(matrix(ncol=3))
colnames(high_lncs) <- c("median", "gene", "canc")

#apply to dataframe 

for(i in 1:length(unique(lnc_rna$canc))){

#subset RNA-dataset to one cancer type
df <- subset(lnc_rna, lnc_rna$canc %in% unique(lnc_rna$canc)[i])

#apply function
res <- apply(df[,1:5607], 2, check_medians)
res <- Filter(Negate(is.null), res)  
res <- data.frame(median=(matrix(unlist(res), nrow=length(res), byrow=T)), gene = names(res), canc=unique(lnc_rna$canc)[i])
high_lncs <- rbind(high_lncs, res)

}#end loop

high_lncs <- high_lncs[-1,]

#---------------------------------------------------------
#Visualize high lncRNA overlap between cancers 
#---------------------------------------------------------

g <- as.data.table(table(high_lncs$gene, high_lncs$canc))
gg <- matrix(nrow=7, ncol=215)
colnames(gg) <- unique(g$V1)
rownames(gg) <- unique(g$V2)

for(i in 1:length(unique(g$V2))){  
  canc <- unique(g$V2)[i]
  z <- which(g$V2 == canc)
  gg[which(rownames(gg)==canc),] <- as.numeric(g$N[z])
}

# Plot 
pdf("215lncRNAs_overlap_Matrix.pdf", pointsize=8, width=9, height=10)
heatmap.2(as.matrix(t(gg)), scale="none", col=mypal[c(2,1)], trace="none", hclustfun = function(x) hclust(x,method = 'ward.D2'),
  distfun = function(x) dist(x,method = 'euclidean'), srtCol=21, cexCol=1.4, key.title=NA, keysize=0.6, 
          margins = c(10, 7), main = list("215 lncRNAs Expression Status in Cancers", cex = 1.2))
dev.off()

#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------

lnc_rna <- lnc_rna[,c((which(colnames(lnc_rna) %in% high_lncs$gene)), 5608,5609)] #215 lncRNAs remain 

#For each patient add survival status and days since last seen 
lnc_rna$status <- ""
lnc_rna$time <- ""
lnc_rna$sex <- ""

#lncs
for(i in 1:nrow(lnc_rna)){
  pat <- rownames(lnc_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  lnc_rna$status[i] <- clin$donor_vital_status[z]
  lnc_rna$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

#---------------------------------------------------------
#Run survival analysis on each gene-cancer combo
#from high_lncs 
#---------------------------------------------------------
results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")

for(i in 1:nrow(high_lncs)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(lnc_rna, lnc_rna$canc %in% high_lncs$canc[i])
  z <- which(colnames(df) %in% high_lncs$gene[i])
  df <- df[,c(z,216:220)]  

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  #median2 <- quantile(as.numeric(df[,1]), 0.5)
  #median2 <- median(df[,1])
  median2 = mean(df[,1])
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
        df$status[df$status=="alive"] <- 0
        df$status[df$status=="deceased"] <- 1
        df$status <- as.numeric(df$status)
        df$time <- as.numeric(df$time)
      
        #cox regression 
        res.cox <- coxph(Surv(time, status) ~ median, data = df)
        row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)], df$canc[1])
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)  
}

results_cox <- results_cox[-1,]
#results_cox$fdr <- p.adjust(results_cox$pval, method="fdr")
#results_cox$fdr <- as.numeric(results_cox$fdr)
#results_cox$pval <- as.numeric(results_cox$pval)

results_cox$fdr = ""
results_cox <- as.data.table(results_cox)
results_cox$pval = as.numeric(results_cox$pval)

cancers = unique(results_cox$canc)

adjustment = function(cancer){
  z <- which(results_cox$canc %in% cancer)  
  #data = filter(data, fdr <= 0.05)
  data = results_cox[z,]
  data$fdr = p.adjust(data$pval, method="fdr")
  z <- which(data$fdr <= 0.1)
  data = data[z,]

  if(!(dim(data)[1] == 0)){
    return(data)
  }
}

adjusted = llply(cancers, adjustment)
adjusted = ldply(adjusted, data.frame)
adjusted = as.data.table(adjusted)
adjusted = adjusted[order(fdr)]

#results_cox <- results_cox[order(fdr)]
#write.table(results_cox, file="results_coxAug28_median5fpkmMin.txt",sep=";", quote=F, row.names=F)

#save as image, dataframe of lncRNAs with pvalue < 0.05 
#sig <- results_cox[pval < 0.05]
#pdf("42_sig_lncRNA_predictors_usingMedian5.pdf", pointsize=8, width=12, height=14)
#p<-tableGrob(sig)
#grid.arrange(p)
#dev.off()

#tier 1 lncRNAs - fdr sig < 0.1
#tier 2 lncRNAs - pvalue < 0.05 

#tier1 <- filter(results_cox, fdr < 0.1) #7 
#tier2 <- filter(results_cox, pval < 0.05 & fdr > 0.1) #35 

#all candidates 
#all <- rbind(tier1, tier2)
#write.table(all, file="7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";", quote=F, row.names=F)
write.table(adjusted, file="lncRNAs_sig_FDR_0.1_Nov23.txt", sep=";", quote=F, row.names=F)


##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("survival_results_usingMean_medGreatThan5_lncRNAs_Nov23.pdf", pointsize=6, width=9, height=8)
require(gridExtra)

for(i in 1:nrow(adjusted)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(lnc_rna, lnc_rna$canc %in% adjusted$canc[i])
  z <- which(colnames(df) %in% adjusted$gene[i])
  df <- df[,c(z,216:220)]  

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  #median2 <- quantile(as.numeric(df[,1]), 0.5)
  #median2 <- median(df[,1])
  median2 = mean(df[,1])
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
  #plot boxplot showing difference between the two groups and sex
  title <- paste(gene, df$canc[1], "Expression")
  colnames(df)[1] <- "Gene"
  g <- ggboxplot(df, x= "median", y="Gene", palette=mypal[c(4,1)], order=c("0", "1"), fill = "median",  add = "jitter")
  g <- g + stat_compare_means()
  g <- ggpar(g, font.legend = c(10, "plain", "black")) 
  g <- g + labs(title = title, y="log1p(FPKM)", x="Mean") + 
      theme(plot.title = element_text(hjust = 0.5))
      print(g)

  #cox
        df$status[df$status=="alive"] <- 0
        df$status[df$status=="deceased"] <- 1
        df$status <- as.numeric(df$status)
        df$time <- as.numeric(df$time)
      
          #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = df)
          s <- ggsurvplot(
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
}


dev.off()






