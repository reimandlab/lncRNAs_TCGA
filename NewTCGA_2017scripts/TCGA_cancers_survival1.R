###---------------------------------------------------------------
###TCGA_cancers_survival1.R
###---------------------------------------------------------------

###October 11th, 2017
###Goal: Using list of candidate lncRNAs obtained in PCAWG analysis
#validate their survival association in this new TCGA data

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(stringr)

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
#cands = fread("7tier1_35tier2_lncRNA_candidates_August28th.txt")

#5. TCGA ID cancer type conversion 
canc_conversion = readRDS("tcga_id_cancer_type_conversion.txt")

#6. List of TCGA IDs used in PCAWG - to remove
ids_remove = fread("TCGA_IDs_usedinPCAWG.txt")

###---------------------------------------------------------------
###Process Data 
###---------------------------------------------------------------

#Change patient ids to shorted id

change = function(rowname){
  new = canc_conversion$id[which(canc_conversion$TCGA_id %in% rowname)]
  return(new)  
}

rownames(rna) = sapply(rownames(rna), change)

#remove thos patients already used in PCAWG
ids_remove = unique(clin$bcr_patient_barcode[which(clin$bcr_patient_barcode %in% ids_remove$bcr_patient_barcode)]) #600 to remove 
z <- which(rownames(rna) %in% ids_remove) #666 PCAWG samples in this TCGA RNA file
rna = rna[-z,]

#Keep only those patients with both RNA-Seq AND clinical data
z <- which(rownames(rna) %in% clin$bcr_patient_barcode)
rna = rna[z,] #all have clinical data - 7387 patients 

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

#change all the cancer type names so they match up
#rna$canc = substr(rna$canc , 1, 4)
#cands$canc = substr(cands$canc , 1, 4)

#subset RNA to candidate genes 
#ensg = fantom$CAT_geneID[which(fantom$CAT_geneName %in% cands$gene)]
#z <- which(colnames(rna) %in% ensg)
#rna = rna[,c(z, 5920:5923)]


###Get list of high expressing lncRNAs in each cancer 

#1. First remove lncRNAs with 0 expression in all patients 
sums = apply(rna[,1:5919], 2, sum)
sums = as.numeric(sums)
z <- which(sums == 0) #134
#remove - MAYBE INSTEAD OF MEDIAN = 0 , SHOULD REMOVE THE ONES 
#THAT HAVE SUM OF 0 MEANING IT HAS 0 EXPRESSION IN EVERY SINGLE PATIENT 
rna = rna[,-z]
#meds = apply(rna, 2, median)
#meds = as.numeric(meds)
#z <- which(meds > 20000)
#remove - potential outliers 
#rna = rna[,-z]

#---------------------------------------------------------
#Find cancer sepcific lncRNAs with median E >= 5 FPKM 
#---------------------------------------------------------

#Write function that takes a dataframe, calculates medians and 
#output list of genes with median greater than that 

check_medians <- function(column){
  med <- median(column)
  if(med >=2000){
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
res <- apply(df[,1:5785], 2, check_medians)
res <- Filter(Negate(is.null), res)  
res <- data.frame(median=(matrix(unlist(res), nrow=length(res), byrow=T)), gene = names(res), canc=unique(rna$canc)[i])
high_lncs <- rbind(high_lncs, res)

}#end loop

high_lncs <- high_lncs[-1,]
#subset rna file to these lncRNAs 
z = which(colnames(rna) %in% high_lncs$gene)
z = c(z, 5786:5789)
rna = rna[,z]
rna[,1:4116] = log1p(rna[,1:4116])

#why some don't have a cancer type?
z <- which(rna$canc == "")
rna = rna[-z,]

z <- which(high_lncs$canc == "")
high_lncs = high_lncs[-z,]

#PCA
library(factoextra)
#res.pca <- prcomp(rna[1:3232], scale = TRUE)
#res.pca <- prcomp(rna[1:3232])

#p = fviz_pca_ind(res.pca,
  #label="none", habillage=rna$canc,
  #           addEllipses=TRUE, ellipse.level=0.95)

#p + theme_minimal()
#dev.off()

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

#For now just to reduce number of tests - I will keep only those genes 
#with median expression greater than the median of medain individual gene expression
#ie = summary(high_lncs$median)$Median == 8252
high_lncs = high_lncs[high_lncs$median >= 8252,]

results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")

survival_analysis = function(row){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% row[[3]])
  ens = fantom$CAT_geneID[which(fantom$CAT_geneID %in% row[[2]])]
  z <- which(colnames(df) %in% ens)
  
  if(!(length(z) ==0)){

  df <- df[,c(z,4117:4120)]  

  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
  if(!(median2 == 0)){
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

  #3. check if there is significant difference in expression between high and low 
  #ie - see that it's not flat 
  p = wilcox.test(df[which(df$median==0),1], df[which(df$median==1),1], paired=F)$p.value
  if(p <= 0.05){

  gene <- colnames(df)[1]
  #cox
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
      
  #cox regression 
  res.cox <- coxph(Surv(time, status) ~ median, data = df)
  row <- data.frame(gene = gene, coef = summary(res.cox)$coefficients[1,1], 
    HR = summary(res.cox)$coefficients[1,2],
    pval = summary(res.cox)$coefficients[1,5],
    canc = df$canc[1])
  if(length(row) > 1){
  names(row) <- names(results_cox)
  print(gene)
  return(row) 
}
}
}
}
}

results = apply(high_lncs, 1, survival_analysis)

#turn list of results into a dataframe 
df <- ldply (results, data.frame)

#need to adjust fdr for each cancer seperatley 
#results_cox = as.data.table(results$cox)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++RESULTS

results = readRDS("results_cox_high_lncs_TCGA.rds")
results$fdr = ""
results = as.data.table(results)
cancers = unique(results$canc)

adjustment = function(cancer){
  data = filter(results, canc %in% cancer)
  data$fdr = p.adjust(data$pval, method="fdr")
  data = filter(data, fdr <= 0.05)
  if(!(dim(data)[1] == 0)){
    return(data)
  }
}

adjusted = llply(cancers, adjustment)
adjusted = ldply (adjusted, data.frame)
adjusted = as.data.table(adjusted)

###Change gene names to Hugo IDs
#for(i in 1:nrow(adjusted)){
	#g = adjusted$gene[i]
	#hugo = fantom$CAT_geneName[which(fantom$CAT_geneID %in% g)]
	#adjusted$gene[i] = hugo
#}

cancs_wlncs = as.data.table(table(adjusted$canc))
cancs_wlncs = cancs_wlncs[order(N)]

##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("TCGA_toplncRNAS_survival_associations.pdf", pointsize=6, width=15, height=14)
require(gridExtra)

adjusted = adjusted[order(fdr)]

for(i in 1:nrow(adjusted)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% adjusted$canc[i])
  z <- which(colnames(df) %in% adjusted$gene[i])
  df <- df[,c(z,4117:4120)]  

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





















































































