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

###canc data
source("ov_source_canc_dataMar21.R")
canc_data = readRDS("OV_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
corlncs = readRDS("OV_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21_mostcorrelated_lncs.rds")

#------FEATURES-----------------------------------------------------
ov_genes_results = readRDS(file="OV_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
ov_features = as.data.table(table(unlist(ov_genes_results)))
ov_features = ov_features[order(N)]
ov_features = dplyr::filter(ov_features, N >=40)
ov_features$canc = "ov"
ov_features$name = ""
for(i in 1:nrow(ov_features)){
  z = which(fantom$gene == ov_features$V1[i])
  ov_features$name[i] = fantom$CAT_geneName[z]
}

z = which(colnames(canc_data) %in% c(ov_features$V1, "canc", "time", "status", "sex", "patient"))
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
canc_data = canc_data[,which(colnames(canc_data) %in% c(ov_features$V1, "time", "status"))]
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)
#keep = names(which(summary(justlncs)$coefficients[,5] <=0.05))
#if(!(length(keep)==0)){
#canc_data = canc_data[,c(which(colnames(canc_data) %in% c(keep,"time", "status")))]
#}
#updated model
justlncs = coxph(Surv(time, status)  ~ ., data = canc_data)

#------PCAWG DATA---------------------------------------------------
pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March20.rds")
pcawg_data = subset(pcawg_data, canc == "Ovary Serous cystadenocarcinoma")
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
pdf("timedependentAUC_externalPCAWG_OV.pdf")
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

pdf("topCands_OV_individual_survivaplots_TCGA.pdf")
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

pdf("topCands_OV_individual_survivaplots_PCAWG.pdf")
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




