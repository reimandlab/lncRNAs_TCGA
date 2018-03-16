###KIRC survival plots

source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

rna = subset(rna, canc %in% c("Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma", 
  "Ovarian serous cystadenocarcinoma", "Pancreatic adenocarcinoma"))

#remove 0 sums
sums = apply(rna[,1:(ncol(rna)-5)], 2, sum)
z = which(sums==0)
#rna = rna[,-(which(colnames(rna) == names(z)))]

#Going to work on each cancer seperatley 
#for now start with one cancer
cancers = cancers[which(cancers %in% rna$canc)] 
cancer = cancers[[3]] #KIRC 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,(ncol(canc_data)-4):ncol(canc_data)]

#New clinical file from Firehose 
newclin = readRDS("537_KIRC_pats_clinical_data_Jan23_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = c(which(str_detect(colss, "stage")), which(str_detect(colss, "age")), which(str_detect(colss, "grade")), which(str_detect(colss, "days")), which(str_detect(colss, "followup")))
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
clin$age = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient.bcr_patient_barcode ==pat)
  clin$stage[i] = newclin$patient.stage_event.pathologic_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
  clin$age[i] = newclin$patient.age_at_initial_pathologic_diagnosis[z]
}

#change stage and grade to numeric values
clin$grade[clin$grade == "g1"] = 1
clin$grade[clin$grade == "g2"] = 2
clin$grade[clin$grade == "g3"] = 3
clin$grade[clin$grade == "g4"] = 4
z = which(clin$grade == "gx")
clin = clin[-z,]

clin$stage[clin$stage == "stage i"] = 1
clin$stage[clin$stage == "stage ii"] = 2
clin$stage[clin$stage == "stage iii"] = 3
clin$stage[clin$stage == "stage iv"] = 4

z = which(is.na(clin$stage))
clin = clin[-z,]
z = which(is.na(clin$grade))
clin = clin[-z,]

clin$grade = as.numeric(clin$grade)
clin$stage = as.numeric(clin$stage)
clin$age = as.numeric(clin$age)
clin$time = as.numeric(clin$time)
clin$status[clin$status=="Alive"] <- 0
clin$status[clin$status=="Dead"] <- 1

############
###START####
############

exp_canc_data = canc_data

for(k in 1:(ncol(canc_data)-5)){
    median2 <- quantile(as.numeric(canc_data[,k]), 0.5)
    if(median2 ==0){
    #if median = 0 then anyone greater than zero is 1 
    for(m in 1:nrow(canc_data)){
    genexp <- canc_data[m,k]
    if(genexp > 0){
      canc_data[m,k] <- 1
      }
    if(genexp == 0 ){
      canc_data[m,k] <- 0
      }
    } 
    }
    if(!(median2 ==0)){
    for(m in 1:nrow(canc_data)){
    genexp <- canc_data[m,k]
    if(genexp >= median2){
      canc_data[m,k] <- 1
      }
    if(genexp < median2){
      canc_data[m,k] <- 0
      }
    } 
    }
}    

###---------------------------------------------------------------
###Cands - 10 batches of 1000 CV 
###---------------------------------------------------------------

genes = readRDS("chosen_features_100CVs_prebinary_predictors_March14.rds")
genes = subset(genes, canc == "kirc")
genes$name = ""
genes = genes[which(str_detect(genes$V1, "ENSG")),]
for(i in 1:nrow(genes)){
  n = fantom$CAT_geneName[which(fantom$gene %in% genes$V1[i])]
  genes$name[i] = n
}

###---------------------------------------------------------------
###Survival analysis 
###---------------------------------------------------------------

#1. Data set-up
cands = genes
genes = as.list(as.character(cands$V1))

#-----------------------------------------------------------------

surv_test = function(gene){
  print(gene)
  results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
  #1. Subset lnc_rna to those patients in cancer
  z <- which(colnames(canc_data) %in% gene)
  if(!(length(z)==0)){
  df = as.data.frame(canc_data)
  df <- df[,c(z,(ncol(canc_data)-4):ncol(canc_data))]  

  gene <- colnames(df)[1]
  colnames(df)[1] = "median"
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
  
  df = merge(df, clin, by=c("patient", "canc", "time", "status", "sex"))
  df$grade <- as.numeric(df$grade)
  df$stage <- as.numeric(df$stage)
  df$age <- as.numeric(df$age)
  #cox regression 
  res.cox <- coxph(Surv(time, status) ~ median + grade + stage + age, data = df)
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
results$canc = "Kidney renal clear cell carcinoma"
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

pdf("11_cands_100CV_KIRC_463patients_March14_preLabelled_binarypred.pdf", pointsize=6, width=10, height=8)
require(gridExtra)

for(i in 1:nrow(results)){
  gene <- results$name[i]
  
  z <- which(colnames(canc_data) %in% results$gene[i])
  if(!(length(z)==0)){
  df = as.data.frame(canc_data)
  df <- df[,c(z,(ncol(canc_data)-4):ncol(canc_data))]  

  colnames(df)[1] = "median"
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
  
  df = merge(df, clin, by=c("patient", "canc", "time", "status", "sex"))
  df$grade <- as.numeric(df$grade)
  df$stage <- as.numeric(df$stage)
  df$age <- as.numeric(df$age)

      #plot boxplot showing difference between the two groups and sex
      title <- paste(gene, df$canc[1], "Expression")
      
      exp = exp_canc_data[,c(which(colnames(exp_canc_data) == results$gene[i]), ncol(exp_canc_data))]
      df = merge(df, exp, by="patient")
      colnames(df)[ncol(df)] = "Gene"
      df$Gene = log1p(df$Gene)
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
pdf("14cands_from10batches_KIRC_1000CV_corrplotJan29.pdf")
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.05, insig = "blank")
dev.off()





































