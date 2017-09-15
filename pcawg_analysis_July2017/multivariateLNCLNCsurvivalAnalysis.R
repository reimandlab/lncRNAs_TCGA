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

#multivaraite analysis
lncsTostudy <- c("LINC00665", "LINC00657")
lncsTostudy <- lnc_rna[,c(which(colnames(lnc_rna) %in% lncsTostudy),216:220)]

#only ovarian
lncsTostudy <- as.data.table(lncsTostudy)
lncsTostudy <- filter(lncsTostudy, canc == "Ovary Serous cystadenocarcinoma")

##-------------------------------------------------------------
##TP53 and survival 
pcg_rna <- pcg_rna[,c((which(colnames(pcg_rna) %in% "TP53")), 19958,19959)] #215 lncRNAs remain 
##-------------------------------------------------------------

#For each patient add survival status and days since last seen 
pcg_rna$status <- ""
pcg_rna$time <- ""
pcg_rna$sex <- ""
for(i in 1:nrow(pcg_rna)){
  pat <- rownames(pcg_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  pcg_rna$status[i] <- clin$donor_vital_status[z]
  pcg_rna$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        pcg_rna$time[i] <- t
}

#only ovarian
lncsTostudy <- as.data.table(pcg_rna)
lncsTostudy <- filter(lncsTostudy, canc == "Ovary Serous cystadenocarcinoma")


#---------------------------------------------------------
#Multivariate survival analysis - Kaplain Mier 
#---------------------------------------------------------

lncsTostudy$status[lncsTostudy$status=="deceased"] = 1
lncsTostudy$status[lncsTostudy$status=="alive"] = 0
lncsTostudy$time <- as.numeric(lncsTostudy$time)
lncsTostudy$status <- as.numeric(lncsTostudy$status)
model4 <- lncsTostudy
lncsTostudy[,1:2] <- log1p(lncsTostudy[,1:2])

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(lncsTostudy, time = "time", event = "status",
   variables = c("LINC00657", "LINC00665"))

# 2. Plot cutpoint for DEPDC1
palette = "npg"
pdf("optimalCutoffs.pdf")
plot(res.cut, "LINC00657", palette = "npg")
plot(res.cut, "LINC00665", palette = "npg")
dev.off()

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize
pdf("survival2LNCS.pdf", pointsize=8, width=10, height=9)

#plot survival plot
          fit <- survfit(Surv(time, status) ~ LINC00665 + LINC00657, data = res.cat)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, 2 top lncRNA Candidates",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = res.cat,      # data used to fit survival curves. 
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

dev.off()


#---------------------------------------------------------
#Multivariate survival analysis - Cox HR 
#---------------------------------------------------------

#ggforest(): Draws forest plot for CoxPH model.

pdf("ggforests_multivaraiteANDunivariateCOX.pdf", pointsize=10, width=10, height=9)

#univariate - LINC00665
model <- coxph( Surv(time, status) ~ LINC00665, 
                data = res.cat )
ggforest(model)
          fit <- survfit(Surv(time, status) ~ LINC00665 , data = res.cat)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, LINC00665",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = res.cat,      # data used to fit survival curves. 
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


#univariate - LINC00657
model <- coxph( Surv(time, status) ~ LINC00657, 
                data = res.cat )
ggforest(model)
          fit <- survfit(Surv(time, status) ~ LINC00657 , data = res.cat)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, LINC00665",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = res.cat,      # data used to fit survival curves. 
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


#multivariate 
model <- coxph( Surv(time, status) ~ LINC00665 + LINC00657, 
                data = res.cat )
ggforest(model)
 fit <- survfit(Surv(time, status) ~ LINC00665 + LINC00657, data = res.cat)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, 2 top lncRNA Candidates",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = res.cat,      # data used to fit survival curves. 
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


dev.off()

##Using continous values instead of tags 
#univariate - LINC00665

pdf("CONTINOUSvalues_ggforests_multivaraiteANDunivariateCOX.pdf", pointsize=10, width=10, height=9)

#univariate - LINC00665
model <- coxph( Surv(time, status) ~ LINC00665, data = lncsTostudy )
ggforest(model)

#univariate - LINC00657
model <- coxph( Surv(time, status) ~ LINC00657, data = lncsTostudy )
ggforest(model)

##multivariate 
model <- coxph( Surv(time, status) ~ LINC00657 + LINC00665, data = lncsTostudy )
ggforest(model)

dev.off()

#---------------------------------------------------------
#Multivariate survival analysis - Using Median 
#---------------------------------------------------------

medianScored <- lncsTostudy
med1 <- median(medianScored$LINC00657)
med2 <- median(medianScored$LINC00665)
medianScored$LINC00657[medianScored$LINC00657 >=med1] = "high"
medianScored$LINC00657[medianScored$LINC00657 <med1] = "low"
medianScored$LINC00665[medianScored$LINC00665 >=med2] = "high"
medianScored$LINC00665[medianScored$LINC00665 <med2] = "low"
test <- medianScored[,c(6,5,1,2)]

pdf("logMinusRatopSept14th_USINGMEDIANSggforests_multivaraiteANDunivariateCOX.pdf", pointsize=10, width=10, height=9)

#using median logRatio between lnc1 and lnc2 
model4$logratio <- log(model4[,1]/model4[,2])
hist(model4$logratio)
med2 <- median(model4$logratio)
model4$logRatio[model4$logratio >= med2] = "high"
model4$logRatio[model4$logratio <  med2] = "low"
test2 <- model4[,c(6,5,9)]

splots <- list()

#univariate - logRATIO between two lncs 
model <- coxph( Surv(time, status) ~ logRatio, 
                data = test2)
ggforest(model)
          fit <- survfit(Surv(time, status) ~ logRatio , data = test2)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, log(lnc1/lnc2)",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = test2,      # data used to fit survival curves. 
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
          splots[[1]] <- s

#univariate - LINC00665
model <- coxph( Surv(time, status) ~ LINC00665, 
                data = test)
ggforest(model)
          fit <- survfit(Surv(time, status) ~ LINC00665 , data = test)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, LINC00665",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = test,      # data used to fit survival curves. 
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
          splots[[2]] <- s


#univariate - LINC00657
model <- coxph( Surv(time, status) ~ LINC00657, 
                data = test )
ggforest(model)
          fit <- survfit(Surv(time, status) ~ LINC00657 , data = test)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, LINC00665",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = test,      # data used to fit survival curves. 
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
          splots[[3]] <- s 


#multivariate 
model <- coxph( Surv(time, status) ~ LINC00665 + LINC00657, 
                data = test )
ggforest(model)
 fit <- survfit(Surv(time, status) ~ LINC00665 + LINC00657, data = test)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, 2 top lncRNA Candidates",       
          #legend.labs = c("H LINC00665, L LINC00657", "High Expression"),             # survfit object with calculated statistics.
          data = test,      # data used to fit survival curves. 
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
          splots[[4]] <- s
 
dev.off()

pdf("4survival_plots.pdf", width=18, height=15, pointsize=8)
 arrange_ggsurvplots(splots, print = TRUE,
  ncol = 2, nrow = 2, risk.table.height = 0.4)
dev.off()



#---------------------------------------------------------
#Comparing models ANOVA 
#---------------------------------------------------------


m1 =  coxph( Surv(time, status) ~ LINC00657, 
                data = test )
m2 = coxph( Surv(time, status) ~ LINC00665, 
                data = test)

#ln1 vs null
anova(m1)

#lnc2 vs null
anova(m2)

#lnc1 vs lnc1+lnc2
m3 =  coxph( Surv(time, status) ~ LINC00657 + LINC00665, 
                data = test )
anova(m3)
anova(m1, m3)

#lnc2 vs lnc1+lnc2
anova(m2, m3)

#lnc1 vs lnc1*lnc2 
m4 =  coxph( Surv(time, status) ~ LINC00657 * LINC00665, 
                data = test )
anova(m4)
anova(m3,m4)
anova(m1,m4)
anova(m2, m4)


##PLOT lnc1 against lnc2 and draw lines on the medians to divide the plot
#into 4 qudrants 

medianScored <- lncsTostudy
med1 <- median(medianScored$LINC00657)
med2 <- median(medianScored$LINC00665)

pdf("lnc1lnc2correlations_quadrants.pdf")

ggplot(medianScored, aes_string(x = "LINC00657", y = "LINC00665")) + 
      geom_point() +
      geom_vline(xintercept = med1) + 
      geom_hline(yintercept = med2)
dev.off()































