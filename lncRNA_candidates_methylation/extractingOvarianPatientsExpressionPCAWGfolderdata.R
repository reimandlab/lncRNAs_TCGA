#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 11th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#last updated August 28th


##----PCAWG folder------------––––––––––––––––––––––––––––

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
lncsTostudy <- c("LINC00665", "ZNF503-AS2", "GS1-251I9.4")
lncsTostudy <- lnc_rna[,c(which(colnames(lnc_rna) %in% lncsTostudy),216:220)]

#only ovarian
lncsTostudy <- as.data.table(lncsTostudy)
lncsTostudy <- filter(lncsTostudy, canc == "Ovary Serous cystadenocarcinoma")

#---------------------------------------------------------
#Multivariate survival analysis - Kaplain Mier 
#---------------------------------------------------------

lncsTostudy$status[lncsTostudy$status=="deceased"] = 1
lncsTostudy$status[lncsTostudy$status=="alive"] = 0
lncsTostudy$time <- as.numeric(lncsTostudy$time)
lncsTostudy$status <- as.numeric(lncsTostudy$status)
model4 <- lncsTostudy
lncsTostudy[,1:3] <- log1p(lncsTostudy[,1:3])

# 1. Determine the optimal cutpoint of variables
colnames(lncsTostudy)[3] = "ZNF503"
colnames(lncsTostudy)[1] = "GS1"

res.cut <- surv_cutpoint(lncsTostudy, time = "time", event = "status", variables = c("ZNF503", "LINC00665", "GS1"))

# 2. Plot cutpoint for DEPDC1
palette = "npg"
pdf("optimalCutoffs.pdf")
plot(res.cut, "ZNF503", palette = "npg")
plot(res.cut, "LINC00665", palette = "npg")
plot(res.cut, "GS1", palette = "npg")
dev.off()

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize
pdf("survival3LNCSZNF503_lnc00665.pdf", pointsize=8, width=10, height=9)

#plot survival plot
          fit <- survfit(Surv(time, status) ~ LINC00665 + ZNF503 + GS1, data = res.cat)
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

pdf("ZNF503ggforests_multivaraiteANDunivariateCOX.pdf", pointsize=10, width=10, height=9)

#multivariate 
model <- coxph( Surv(time, status) ~ LINC00665 + ZNF503 + GS1, 
                data = res.cat )
ggforest(model)
 fit <- survfit(Surv(time, status) ~ LINC00665 + ZNF503 + GS1, data = res.cat)
          s <- ggsurvplot(
          fit, 
          main = "Survival Curves, 3 top lncRNA Candidates",       
          #legend.labs = c("H LINC00665, L ZNF503", "High Expression"),             # survfit object with calculated statistics.
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
#Multivariate survival analysis - Using Median 
#---------------------------------------------------------

medianScored <- lncsTostudy
med1 <- median(medianScored$ZNF503)
med2 <- median(medianScored$LINC00665)
med3 <- median(medianScored$GS1)
medianScored$ZNF503[medianScored$ZNF503 >=med1] = "high"
medianScored$ZNF503[medianScored$ZNF503 <med1] = "low"
medianScored$LINC00665[medianScored$LINC00665 >=med2] = "high"
medianScored$LINC00665[medianScored$LINC00665 <med2] = "low"
medianScored$GS1[medianScored$GS1 >=med3] = "high"
medianScored$GS1[medianScored$GS1 <med3] = "low"

#Save patient status for each lncRNA 
write.table(medianScored, file="medianScoresOvarianCancerTop3_lncRNAs.txt", quote=F, sep=";")

test <- medianScored[,c(7,6,1,2,3)]

#---------------------------------------------------------
#Comparing models ANOVA 
#---------------------------------------------------------


m1 =  coxph( Surv(time, status) ~ ZNF503, 
                data = test )
m2 = coxph( Surv(time, status) ~ LINC00665, 
                data = test)

m3 = coxph( Surv(time, status) ~ GS1, 
                data = test)

#lnc1 vs lnc1+lnc2
m4 =  coxph( Surv(time, status) ~ ZNF503 + LINC00665 + GS1, 
                data = test)
anova(m1, m4)

#lnc2 vs lnc1+lnc2
anova(m2, m4)

anova(m3,m4)


##PLOT lnc1 against lnc2 and draw lines on the medians to divide the plot
#into 4 qudrants 

medianScored <- lncsTostudy
med1 <- median(medianScored$ZNF503)
med2 <- median(medianScored$LINC00665)

pdf("lnc1lnc2correlations_quadrants.pdf")

ggplot(medianScored, aes_string(x = "ZNF503", y = "LINC00665")) + 
      geom_point() +
      geom_vline(xintercept = med1) + 
      geom_hline(yintercept = med2)
dev.off()































