#!/usr/bin/Rscript

###---------------------------------------------------------------
###TCGA_cancers_survival1.R
###---------------------------------------------------------------

###October 11th, 2017
###Goal: Using all high expressing lncRNAs in dataset for top 5 cancers 
###analyze their assocation to survival

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. RNA data 

rna = readRDS("5919_lncs4cancers_TCGAnew.rds")
rownames(rna) = rna$gene
rna$gene = NULL
rna = t(rna)

#2. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]

#3. Fantom data 
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

#4. List of lncRNA survival associated candidates 
cands = read.csv("42_original_candidates.csv", header=F)
colnames(cands) = c("lncRNA", "cancer")

#5. TCGA ID cancer type conversion 
canc_conversion = readRDS("tcga_id_cancer_type_conversion.txt")

#6. List of TCGA IDs used in PCAWG - to remove
ids_remove = fread("819_unique_TCGAids_usedbyPCAWG.txt")


###---------------------------------------------------------------
###Process Data 
###---------------------------------------------------------------

#Change patient ids to shorted id
for(i in 1:nrow(rna)){
	rowname = rownames(rna)[i]
	new = canc_conversion$id[which(canc_conversion$TCGA_id %in% rowname)]
	rownames(rna)[i] = new
}

#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_uuid %in% ids_remove$V2)]) #600 to remove 
z <- which(rownames(rna) %in% ids_remove) #69 PCAWG samples in this TCGA RNA file
rna = rna[-z,]

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(rna) %in% clin$bcr_patient_barcode)
rna = rna[z,] #all have clinical data

#Add survival info to rna file
rna = as.data.frame(rna)
rna$canc = "" 
rna$time = ""
rna$status = ""
rna$sex = ""

for(i in 1:dim(rna)[1]){
	pat = rownames(rna)[i]
	z <- which(clin$bcr_patient_barcode == pat)
	if(length(z) >1){z = z[1]}
	status = clin$vital_status[z]
	time = clin$days_to_death[z]
	if(is.na(time)){
		time = clin$days_to_last_followup[z]
	}
	rna$status[i] = status
	rna$time[i] = time
	rna$sex[i] = clin$gender[z]
}

for(i in 1:dim(rna)[1]){
	pat = rownames(rna)[i]
	z = which(canc_conversion$id == pat)
	canc = canc_conversion$Cancer[z]
	rna$canc[i] = canc
}


###---------------------------------------------------------------
###Looking at candidate lncRNAs
###---------------------------------------------------------------

#Find cancer sepcific lncRNAs with median E >= 5 FPKM 
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
for(i in 1:length(unique(rna$canc))){

#subset RNA-dataset to one cancer type
df <- subset(rna, rna$canc %in% unique(rna$canc)[i])

#apply function
res <- apply(df[,1:5919], 2, check_medians)
res <- Filter(Negate(is.null), res)  
res <- data.frame(median=(matrix(unlist(res), nrow=length(res), byrow=T)), gene = names(res), canc=unique(rna$canc)[i])
high_lncs <- rbind(high_lncs, res)

}#end loop

high_lncs <- high_lncs[-1,]

#---------------------------------------------------------
#Visualize high lncRNA overlap between cancers 
#---------------------------------------------------------

g <- as.data.table(table(high_lncs$gene, high_lncs$canc))
gg <- matrix(nrow=4, ncol=4617)
colnames(gg) <- unique(g$V1)
rownames(gg) <- unique(g$V2)

for(i in 1:length(unique(g$V2))){  
  canc <- unique(g$V2)[i]
  z <- which(g$V2 == canc)
  gg[which(rownames(gg)==canc),] <- as.numeric(g$N[z])
}

# Plot 
pdf("4617lncRNAs_overlap_Matrix.pdf", pointsize=8, width=9, height=10)
heatmap.2(as.matrix(t(gg)), scale="none", col=mypal[c(2,1)], trace="none", hclustfun = function(x) hclust(x,method = 'ward.D2'),
  distfun = function(x) dist(x,method = 'euclidean'), srtCol=21, cexCol=1.4, key.title=NA, keysize=0.6, 
          margins = c(10, 7), main = list("4,617 lncRNAs Expression Status in Cancers", cex = 1.2))
dev.off()

#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------

rna <- rna[,c((which(colnames(rna) %in% high_lncs$gene)), 5920,5923)] #4617 lncRNAs remain 

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")

for(i in 1:nrow(high_lncs)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% high_lncs$canc[i])
  z <- which(colnames(df) %in% high_lncs$gene[i])
  df <- df[,c(z,4618:4621)]  

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
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

results_cox <- results_cox[-1,]
results_cox$fdr <- p.adjust(results_cox$pval, method="fdr")
results_cox$fdr <- as.numeric(results_cox$fdr)
results_cox$pval <- as.numeric(results_cox$pval)

results_cox <- as.data.table(results_cox)
results_cox <- results_cox[order(fdr)]
write.table(results_cox, file="4617lncRNAs_results_Octo16_median5FPKM_TCGAanalysis.txt",sep=";", quote=F, row.names=F)


###Change gene names to Hugo IDs
#for(i in 1:nrow(results_cox)){
	#g = results_cox$gene[i]
	#hugo = fantom$CAT_geneName[which(fantom$CAT_geneID %in% g)]
	#results_cox$gene[i] = hugo
#}


##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("TCGA_survival_highExplncRNAs_results_4617lncRNAs.pdf", pointsize=6, width=15, height=14)
require(gridExtra)

for(i in 1:nrow(results_cox)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% high_lncs$canc[i])
  z <- which(colnames(df) %in% high_lncs$gene[i])
  df <- df[,c(z,4618:4621)] 

  df[,1] <- log1p(df[,1])

  #3. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
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
  gene = fantom$CAT_geneName[which(fantom$CAT_geneID %in% gene)]

  #plot boxplot showing difference between the two groups and sex
  title <- paste(gene, df$canc[1], "Expression")
  colnames(df)[1] <- "Gene"
  g <- ggboxplot(df, x= "median", y="Gene", palette=mypal, order=c("0", "1"), fill = "median",  add = "jitter")
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
      
          #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = df)
          s <- ggsurvplot(
          fit, 
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
          palette = colorRampPalette(mypal)(14), 
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)  
}


dev.off()





















































































