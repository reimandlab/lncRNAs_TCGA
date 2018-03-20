library(survAUC)

source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#start with only lncRNA_intergenic
#lincs = subset(fantom, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))

lnc_info = read.csv("fantom_genebased_evidence_supp_table_17.csv")
lnc_info = lnc_info[which(str_detect(lnc_info$CAT_geneID, "ENSG")),]
lnc_info = subset(lnc_info, lnc_info$CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_divergent", "p_lncRNA_intergenic"))

#shorten the gene names 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
lnc_info[,1] <- apply(lnc_info[,1:2], 1, extract3) #5049 lncRNAs 
lnc_info = subset(lnc_info, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))
#2091 

#how many conserved 
#1207 have conserved exons 
z1 = which(!(lnc_info$exon_RS_score == "__na"))

#924 have conserved TIRs 
z2 = which(!(lnc_info$TIR_RS_score == "__na"))

conserved = unique(c(z1, z2))
lnc_info = lnc_info[conserved, ]
colnames(lnc_info)[1] = "gene"

#z = which(colnames(rna) %in% lnc_info$gene) <-------- 
rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))] <---------

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

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

#how many lncRNAs have median < 1? 
canc_data[,1:(ncol(canc_data)-5)] = log1p(canc_data[,1:(ncol(canc_data)-5)])

detectable = readRDS("PCAWG_detectable_genes_4cancers_March20.rds")
detectable = subset(detectable, canc == "Kidney-RCC")
detectable$canc = "Kidney renal clear cell carcinoma"

z = which(colnames(canc_data) %in% c(detectable$gene, "canc", "time", "status", "sex", "patient"))
canc_data = canc_data[,z]

#------FEATURES-----------------------------------------------------
kirc_genes_results = readRDS(file="ALL_KIRC_pats_prebin_predictors_list_of_sig_genes_CV_100March19nokeep.rds")
kirc_features = as.data.table(table(unlist(kirc_genes_results)))
kirc_features = kirc_features[order(N)]
kirc_features = dplyr::filter(kirc_features, N >=40)
kirc_features$canc = "kirc"
kirc_features$name = ""
for(i in 1:nrow(kirc_features)){
	z = which(fantom$gene == kirc_features$V1[i])
	kirc_features$name[i] = fantom$CAT_geneName[z]
}

z = which(colnames(canc_data) %in% c(kirc_features$V1, "canc", "time", "status", "sex", "patient"))
canc_data = canc_data[,z]

#add high low tag
medians = apply(canc_data[,1:(ncol(canc_data)-5)], 2, median)
for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(canc_data[,k] > 0)
    l2 = which(canc_data[,k] ==0)
    canc_data[l1,k] = 1
    canc_data[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(canc_data[,k] >= med)
    l2 = which(canc_data[,k] < med)
    canc_data[l1,k] = 1
    canc_data[l2, k] = 0
    }
}  

canc_data$time = as.numeric(canc_data$time)
canc_data$status[canc_data$status=="Alive"] <- 0
canc_data$status[canc_data$status=="Dead"] <- 1
canc_data$status = as.numeric(canc_data$status)

#####Train model using all TCGA data and the chosen predictor lncRNAs 
canc_data = canc_data[,which(colnames(canc_data) %in% c(kirc_features$V1, "time", "status"))]
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)
#keep = names(which(summary(justlncs)$coefficients[,5] <=0.05))
#if(!(length(keep)==0)){
#canc_data = canc_data[,c(which(colnames(canc_data) %in% c(keep,"time", "status")))]
#}
#updated model
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)

#------PCAWG DATA---------------------------------------------------
pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March20.rds")
pcawg_data = subset(pcawg_data, canc == "Kidney Adenocarcinoma, clear cell type")
z = which(colnames(pcawg_data) %in% c(colnames(canc_data), "time", "status"))
pcawg_data = pcawg_data[,z]
#add high low tag
medians = apply(pcawg_data[,1:(ncol(pcawg_data)-2)], 2, median)
for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(pcawg_data[,k] > 0)
    l2 = which(pcawg_data[,k] ==0)
    pcawg_data[l1,k] = 1
    pcawg_data[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(pcawg_data[,k] >= med)
    l2 = which(pcawg_data[,k] < med)
    pcawg_data[l1,k] = 1
    pcawg_data[l2, k] = 0
    }
}  

pcawg_data$time = as.numeric(pcawg_data$time)
pcawg_data$status = as.numeric(pcawg_data$status)

####Test on the PCAWG data

#using all 8 candidate lncRNAs 
pdf("timedependentAUC_externalPCAWG_KIRC.pdf")
lpnew <- predict(justlncs, newdata=pcawg_data)
Surv.rsp <- Surv(canc_data$time, canc_data$status)
Surv.rsp.new <- Surv(pcawg_data$time, pcawg_data$status)
times <- seq(10, 1000, 10)
AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
names(AUC_Uno)
AUC_Uno$iauc
plot(AUC_Uno)
abline(h = 0.5)
BrierScore <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times,
type = "brier", int.type = "weighted")
plot(BrierScore)
abline(h = 0.25)
dev.off()

#------individual lncs-----------------------------------------

pdf("topCands_KIRC_individual_survivaplots_TCGA.pdf")
results_cox1 <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
for(i in 1:(ncol(canc_data)-2)){
	dat = canc_data[,c(i, ncol(canc_data), (ncol(canc_data)-1))]
	justlncs_pcawg = coxph(Surv(time, status)  ~ ., data = dat)
	row <- c(colnames(canc_data)[i], summary(justlncs_pcawg)$coefficients[1,c(1,2,5)],  summary(justlncs_pcawg)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox1)	
  	results_cox1 = rbind(results_cox1, row)
  	gene = colnames(canc_data)[i]
  	colnames(dat)[1] = "gene"
  	dat$time = dat$time/365
  	fit <- survfit(Surv(time, status) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene, "HR =", round(as.numeric(row[3]), digits=4)),
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
          data = dat,      # data used to fit survival curves. 
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
dev.off()
results_cox1 = results_cox1[-1,]

#------individual lncs-----------------------------------------

pdf("topCands_KIRC_individual_survivaplots_PCAWG.pdf")
results_cox2 <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox2) <- c("gene", "coef", "HR", "pval", "low95", "upper95")
for(i in 1:(ncol(pcawg_data)-2)){
	dat = pcawg_data[,c(i, ncol(pcawg_data), (ncol(pcawg_data)-1))]
	justlncs_pcawg = coxph(Surv(time, status)  ~ ., data = dat)
	row <- c(colnames(pcawg_data)[i], summary(justlncs_pcawg)$coefficients[1,c(1,2,5)],  summary(justlncs_pcawg)$conf.int[1,c(3,4)])
  	names(row) <- names(results_cox2)	
  	results_cox2 = rbind(results_cox2, row)
  	gene = colnames(pcawg_data)[i]
  	colnames(dat)[1] = "gene"
  	dat$time = dat$time/365
  	fit <- survfit(Surv(time, status) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene, "HR =", round(as.numeric(row[3]), digits=4)),
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
          data = dat,      # data used to fit survival curves. 
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
dev.off()
results_cox2 = results_cox2[-1,]




