library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

#rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
#pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#--------This script ------------------------------------------------

#Train and evaluate the performance of a logistic regression model
#using expression of 166 lncRNAs to predict tumour type
#features --> expression of 166 lncRNAs
#outcome --> tumour type
#binary expression? yess.. 

#--------------------------------------------------------------------

z = which(cancers %in% allCands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "type"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

add_tags = function(dtt){
  print(dtt$Cancer[1])
  rownames(dtt) = dtt$patient
  dtt$patient = NULL

  #log1p 
  z = which(str_detect(colnames(dtt), "ENSG"))
  if(length(z)>1){
  medians = apply(dtt[,z], 2, median)}
  if(length(z)==1){
    medians = median(dtt[,z])
  }
  #add high low tag
  for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(dtt[,k] > 0)
    l2 = which(dtt[,k] ==0)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(dtt[,k] >= med)
    l2 = which(dtt[,k] < med)
    dtt[l1,k] = 1
    dtt[l2, k] = 0
    }
  }  
  return(dtt)
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")

#-------Combine all data into one dataframe-------------------------------


all_canc_data = ldply(filtered_data_tagged)


#-------Set up Cross-Validation------------------------------------------

run_cv = function(dtt){
  
  set.seed(101) 

  all_runs_results = as.data.frame(matrix(ncol=2))
  colnames(all_runs_results) = c("C-index", "round")

  for(j in 1:100){
    
    print(j)
    smp_size <- floor(0.7 * nrow(dtt))
    train_ind <- sample(seq_len(nrow(dtt)), size = smp_size)
    train <- dtt[train_ind, ]
    test <- dtt[-train_ind, ]

    #keep in the training set only the genes and tumour type 
    z1 = which(str_detect(colnames(train), "ENSG"))
    z2 = which(colnames(train) %in% "type")
    train = train[,c(z1,z2)]

}
}

#dtt = all_canc_data















