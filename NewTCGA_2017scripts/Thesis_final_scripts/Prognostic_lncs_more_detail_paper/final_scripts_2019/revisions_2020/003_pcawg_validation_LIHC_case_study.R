set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
saveRDS(allCands, "/u/kisaev/final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

tcga_results1 = readRDS("TCGA_results_multivariate_results_Oct3.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")
tcga_results1$clinical_only_concordance = NULL
tcga_results1$num_events = NULL
tcga_results1$perc_wevents = NULL
tcga_results1 = as.data.table(tcga_results1)
saveRDS(tcga_results1, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#--------------------------------------------------------------------
#DONT RUN BELOW#
#--------------------------------------------------------------------

pcawg_clin = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/top5_cancers_pcawg_thesis_jult2017/raw/Jan26_PCAWG_clinical")
pcawg_data = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNA_clinical_data_PCAWG_May2.rds")
pcawg_data$combo = paste(pcawg_data$canc, pcawg_data$histo)

pcawg_data$canc[pcawg_data$combo == "Uterus Adenocarcinoma, endometrioid"] = "Uterine Corpus Endometrial Carcinoma"
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, clear cell type"] = "Kidney renal clear cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Breast Infiltrating duct carcinoma"] = "Breast invasive carcinoma"
pcawg_data$canc[pcawg_data$combo == "Ovary Serous cystadenocarcinoma"] = "Ovarian serous cystadenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "Pancreas Pancreatic ductal carcinoma"] = "Pancreatic adenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "Liver Hepatocellular carcinoma"] = "Liver hepatocellular carcinoma"
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, papillary type"] = "Kidney renal papillary cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Squamous cell carcinoma"] = "Lung squamous cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Adenocarcinoma, invasive"] = "Lung adenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "CNS Glioblastoma"] = "Glioblastoma multiforme "
pcawg_data$canc[pcawg_data$combo == "Esophagus Adenocarcinoma"] = "Esophageal carcinoma "
pcawg_data$canc[pcawg_data$canc == "Colon/Rectum"] = "Colon adenocarcinoma"
pcawg_data$canc[pcawg_data$canc == "Stomach"] = "Stomach adenocarcinoma"
pcawg_data$canc[pcawg_data$canc == "Head/Neck"] = "Head and Neck squamous cell carcinoma"
pcawg_data$canc[pcawg_data$canc == "Cervix"] = "Cervical squamous cell carcinoma and endocervical adenocarcinoma"

#remove those already used in TCGA cohort - aka those with a TCGA id
ids = unique(pcawg_clin[,c("icgc_donor_id", "submitted_donor_id.y")])
colnames(ids) = c("patient", "TCGA_id")
pcawg_data = merge(pcawg_data, ids, by= "patient")
z = which(str_detect(pcawg_data$TCGA_id, "TCGA"))
pcawg_data = pcawg_data[-z,]

#add more pcawg people
cancers_tests = as.list(unique(tcga_results1$cancer[which(tcga_results1$cancer %in% pcawg_data$canc)]))

get_matched_data = function(cancer){
    dtt = subset(pcawg_data, canc == cancer)
    z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$canc[1]]), "canc",
    "histo", "time", "status", "sex", "donor_age_at_diagnosis", "patient"))
    dtt = dtt[,z]
    if(nrow(dtt) >= 20){
    return(dtt)}
}

#save LIHC data
lihc_dat = get_matched_data("Liver hepatocellular carcinoma")
saveRDS(lihc_dat, file="/u/kisaev/LIHC_PCAWG_data.rds")

#save TCGA LIHC data
z = which(rna$type == "LIHC")
lihc_tcga = rna[z,]
z = which(colnames(lihc_tcga) %in% c("ENSG00000231918",
  "ENSG00000178977", "ENSG00000264026", "ENSG00000258731",
  "ENSG00000227685", "ENSG00000230490", "ENSG00000236760",
  "ENSG00000236437", "ENSG00000237027",
 "OS", "OS.time", "patient"))
lihc_tcga = lihc_tcga[,..z]
saveRDS(lihc_tcga, file="/u/kisaev/LIHC_TCGA_data_2021.rds")
