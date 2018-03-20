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

#------PCAWG DATA---------------------------------------------------
pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March2018.rds")
pcawg_data = subset(pcawg_data, canc == "Kidney Adenocarcinoma, clear cell type")
z = which(colnames(pcawg_data) %in% c(colnames(canc_data), "time", "status"))
pcawg_data = pcawg_data[,z]

lncs = as.list(colnames(canc_data)[1:544])

pdf("lnc_exp_correlation_bw_pcawg_tcga_mar20.pdf")

check_cor_lnc = function(lnc){
	p = pcawg_data[,which(colnames(pcawg_data) == lnc)]
	p = log1p(p)
	p = as.data.frame(p) 
	colnames(p)[1] = "pcawg"
	z = dim(p)[1]
	t = canc_data[,which(colnames(canc_data) == lnc)]
	t = as.data.frame(t)
	colnames(t)[1] = "tcga"
	r = c()
	for(k in 1:5){
	train_ind <- sample(seq_len(nrow(t)), size = z)
	t = t[train_ind,]
	t = as.data.frame(t)
	colnames(t)[1] = "tcga"
	lnc_dat = cbind(p, t)
	library("ggpubr")
	g = ggscatter(lnc_dat, x = "pcawg", y = "tcga", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pcawg log1p(fpkm)", ylab = "tcga log1p(fpkm)", main = paste(lnc, "correlation"))
	g = g + geom_hline(yintercept = median(lnc_dat$tcga), linetype = 2, colour="red") + geom_vline(xintercept=median(lnc_dat$pcawg), 
		linetype=2, colour="orange")

	print(g)
	r = c(r, cor(lnc_dat$pcawg, lnc_dat$tcga))
	}
	r = mean(r)
	return(r)
}

rcors = llply(lncs, check_cor_lnc)
dev.off()
rcors = unlist(rcors)
which(rcors >=0.2)


