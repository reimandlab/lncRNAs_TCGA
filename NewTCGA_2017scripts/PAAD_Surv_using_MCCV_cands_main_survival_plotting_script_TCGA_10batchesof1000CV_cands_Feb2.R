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

source("source_code_Cox_MonteCarlo_CV_PAAD_Feb5.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5749:5753)]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095

library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]] #PAAD
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,2336:ncol(canc_data)]

#New clinical file from Firehose 
#New clinical file from Firehose 
newclin = readRDS("185_PAAD_pats_clinical_data_Feb5_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = c(which(str_detect(colss, "stage")), which(str_detect(colss, "grade")), which(str_detect(colss, "days")), which(str_detect(colss, "followup")))
newclin = newclin[,c(1:25, colsskeep)]

#see if time and status mathch old information 
clin$newtime = ""
clin$newstatus =""
for(i in 1:nrow(newclin)){
  pat = newclin$patient.bcr_patient_barcode[i]
  z = which(clin$patient %in% pat)
  t = newclin$patient.days_to_death[i]
  if(!(is.na(t))){
    clin$newtime[z] = t
    clin$newstatus[z] = "Dead"
  }
  if(is.na(t)){
    clin$newtime[z] = newclin$patient.days_to_last_followup[i]
    clin$newstatus[z] = "Alive"
  }
}

z = which(is.na(clin$time))
#clin = clin[-z,]
ind <- clin$time == clin$newtime
clin[ind, "newtime"] <- "EQUAL"
ind <- clin$status == clin$newstatus
clin[ind, "newstatus"] <- "EQUAL"

clin = subset(clin, time >0)

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% clin$patient)

clin$grade = ""
clin$stage = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient.bcr_patient_barcode ==pat)
  clin$stage[i] = newclin$patient.stage_event.pathologic_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
}

#292 LIHC patients with clinical data 

#change stage and grade to numeric values
clin$grade[clin$grade == "g1"] = 1
clin$grade[clin$grade == "g2"] = 2
clin$grade[clin$grade == "g3"] = 3
clin$grade[clin$grade == "g4"] = 4
z = which(clin$grade == "gx")
clin = clin[-z,]

clin$stage[clin$stage == "stage i"] = 1
clin$stage[clin$stage == "stage ia"] = 1
clin$stage[clin$stage == "stage ib"] = 1

clin$stage[clin$stage == "stage iia"] = 2
clin$stage[clin$stage == "stage iib"] = 2

clin$stage[clin$stage == "stage iii"] = 3
clin$stage[clin$stage == "stage iv"] = 4

z = which(is.na(clin$stage))
clin = clin[-z,]
z = which(is.na(clin$grade))
#clin = clin[-z,]

clin$grade = as.numeric(clin$grade)
clin$stage = as.numeric(clin$stage)

###---------------------------------------------------------------
###Cands - 10 batches of 1000 CV 
###---------------------------------------------------------------

genes = readRDS("final_candidates_10batches_of1000CV_PAAD_cancer_patientsFeb6.RDS")
cands = as.data.frame(table(genes$gene))
colnames(cands)[1] = "gene"
cands$name = ""
for(i in 1:nrow(cands)){
  n = fantom$CAT_geneName[which(fantom$gene %in% cands$gene[i])]
  cands$name[i] = n
}

#only keep those that appreared in 9-10 of 10 batches 
cands = subset(cands, Freq >=9)

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

#1. Data set-up
genes = as.list(as.character(cands$gene))
canc_data[,1:2335] = log1p(canc_data[,1:2335])
sums = apply(canc_data[,1:2335], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
z <- which(colnames(canc_data) %in% zeroes)
if(!(length(z)==0)){
  canc_data = canc_data[,-z]
}

canc_data = canc_data[,-c(2335:2338)]
canc_data = merge(canc_data, clin, by=c("patient"))
rownames(canc_data) = canc_data$patient
canc_data = canc_data[,-1]

#-----------------------------------------------------------------

surv_test = function(gene){
  print(gene)
  results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
  #1. Subset lnc_rna to those patients in cancer
  z <- which(colnames(canc_data) %in% gene)
  if(!(length(z)==0)){
  df = as.data.frame(canc_data)
  df <- df[,c(z,(ncol(canc_data)-7):ncol(canc_data))]  

  #2. Add Median cutoff tag High or Low to each patient per each gene 
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
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
      
  #cox regression 
  res.cox <- coxph(Surv(time, status) ~ median + grade + stage, data = df)
  row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)])
  names(row) <- names(results_cox)
}
 return(row)
}

results = llply(genes, surv_test)
results = ldply(results, rbind)
results$pval = as.numeric(results$pval)
results$pval = round(results$pval, digits=6)
results$low95 = as.numeric(results$low95)
results$low95 = round(results$low95, digits=4)
results$upper95 = as.numeric(results$upper95)
results$upper95 = round(results$upper95, digits=4)
results$canc = "Pancreatic adenocarcinoma"
results$name = ""
for(i in 1:nrow(results)){
  n = fantom$CAT_geneName[which(fantom$gene %in% results$gene[i])]
  results$name[i] = n
}

##+++++++++++++++++++++++++++++
##Full - order plots by increasing absolute HR 
##+++++++++++++++++++++++++++++
results$HR = as.numeric(results$HR)
results = as.data.table(results)
sort.field = "HR"
sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
results = sortme(results, sort.field)
results = as.data.frame(results)

pdf("2_cands_10batches_1000CV_PAAD_170ishpatients_Feb6.pdf", pointsize=6, width=10, height=8)
require(gridExtra)

for(i in 1:nrow(results)){
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% results$canc[i])
  #gene = fantom$gene[which(fantom$CAT_geneName %in% results_cox$name[i])]
  z <- which(colnames(df) %in% results$gene[i])
  if(!(length(z)==0)){
  df <- df[,c(z,2336:ncol(df))]  

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

      gene <- results$name[i]
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
          title = paste(gene, "HR =", round(results$HR[i], digits=4), df$canc[1]),
          fit, 
          xlab = "Time (Years)", 
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


###---------------------------------------------------------------
###Are these 7 genes correlated with each other? 
###---------------------------------------------------------------

lncrnas = unlist(genes)
lncrnas = canc_data[,which(colnames(canc_data) %in% lncrnas)]
for(i in 1:ncol(lncrnas)){
  colnames(lncrnas)[i] = fantom$CAT_geneName[which(fantom$gene %in% colnames(lncrnas)[i])]
}

library("Hmisc")
res2 <- rcorr(as.matrix(lncrnas), type = c("spearman"))
res2$r
# Extract p-values
res2$P
# Insignificant correlation are crossed
library(corrplot)
pdf("2cands_from10batches_PAAD_1000CV_corrplotFeb6.pdf")
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.05, insig = "blank")
dev.off()





































