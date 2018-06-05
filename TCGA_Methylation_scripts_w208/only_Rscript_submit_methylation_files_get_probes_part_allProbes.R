###methylation files 
###identify probes overlapping lncRNAs 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#############in command line###################################################################################################################################
#1
#head -n 1 OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > OV_patients.txt
#ov methylation data --> OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt 

#1in R
ov_pats = read.table("OV_patients.txt")
ov_pats = ov_pats[,3:ncol(ov_pats)]
ov_pats = as.character(unlist(ov_pats))
ov_meth = fread("OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(ov_meth) = c("probe", ov_pats)
saveRDS(ov_meth, file="OV_methylation_data_ALLprobes.rds")

#2
#head -n 1 PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > PAAD_patients.txt

#2in R
paad_pats = read.table("PAAD_patients.txt")
paad_pats = paad_pats[,3:ncol(paad_pats)]
paad_pats = as.character(unlist(paad_pats))
paad_meth = fread("PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(paad_meth) = c("probe", paad_pats)
saveRDS(paad_meth, file="PAAD_methylation_data_ALLprobes.rds")

#3
#head -n 1 KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRC_patients.txt

#in R
kirc_pats = read.table("KIRC_patients.txt")
kirc_pats = kirc_pats[,3:ncol(kirc_pats)]
kirc_pats = as.character(unlist(kirc_pats))
kirc_meth = fread("KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(kirc_meth) = c("probe", kirc_pats)
saveRDS(kirc_meth, file="KIRC_methylation_data_ALLprobes.rds")

#4
#head -n 1 LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LIHC_patients.txt

lihc_pats = read.table("LIHC_patients.txt")
lihc_pats = lihc_pats[,3:ncol(lihc_pats)]
lihc_pats = as.character(unlist(lihc_pats))
lihc_meth = fread("LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(lihc_meth) = c("probe", lihc_pats)
saveRDS(lihc_meth, file="LIHC_methylation_data_ALLprobes.rds")

#lung adeno
#LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt

#5
#head -n 1 LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUAD_patients.txt

luad_pats = read.table("LUAD_patients.txt")
luad_pats = luad_pats[,3:ncol(luad_pats)]
luad_pats = as.character(unlist(luad_pats))
luad_meth = fread("LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(luad_meth) = c("probe", luad_pats)
saveRDS(luad_meth, file="LUADmethylation_data_ALLprobes.rds")

#breast
#BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt

#6
#head -n 1 BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BRCA_patients.txt

brca_pats = read.table("BRCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="BRCAmethylation_data_ALLprobes.rds")

#7ACC

#head -n 1 ACC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ACC_patients.txt

brca_pats = read.table("ACC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("ACC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="ACCmethylation_data_ALLprobes.rds")

#8 ESCA

#head -n 1 ESCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ESCA_patients.txt

brca_pats = read.table("ESCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("ESCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="ESCAmethylation_data_ALLprobes.rds")

#9 GBM 

#head -n 1 GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > GBM_patients.txt

brca_pats = read.table("GBM_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="GBMmethylation_data_ALLprobes.rds")

#10 HNSC

#head -n 1 HNSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > HNSC_patients.txt

brca_pats = read.table("HNSC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("HNSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="HNSCmethylation_data_ALLprobes.rds")

#11 LGG 

#head -n 1 LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LGG_patients.txt
#grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LGG_lnc_probes.txt

brca_pats = read.table("LGG_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="LGGmethylation_data_ALLprobes.rds")

#12 MESO 

#head -n 1 MESO.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > MESO_patients.txt

brca_pats = read.table("MESO_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("MESO.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="MESOmethylation_data_ALLprobes.rds")

#13 READ

#head -n 1 READ.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > READ_patients.txt

brca_pats = read.table("READ_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("READ.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="READmethylation_data_ALLprobes.rds")

#14 THCA 

#head -n 1 THCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > THCA_patients.txt

brca_pats = read.table("THCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("THCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="THCAmethylation_data_ALLprobes.rds")

#15 UVM 

#head -n 1 UVM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UVM_patients.txt

brca_pats = read.table("UVM_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("UVM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="UVMmethylation_data_ALLprobes.rds")

#16 CESC

#head -n 1 CESC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > CESC_patients.txt

brca_pats = read.table("CESC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("CESC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="CESCmethylation_data_ALLprobes.rds")

#17 KIRP

#head -n 1 KIRP.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRP_patients.txt

brca_pats = read.table("KIRP_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("KIRP.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="KIRPmethylation_data_ALLprobes.rds")

#18 STAD

#head -n 1 STAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > STAD_patients.txt

brca_pats = read.table("STAD_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("STAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="STADmethylation_data_ALLprobes.rds")

#19 LUSC

#head -n 1 LUSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUSC_patients.txt

brca_pats = read.table("LUSC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("LUSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="LUSCmethylation_data_ALLprobes.rds")

#20 BLCA 

#head -n 1 BLCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BLCA_patients.txt

brca_pats = read.table("BLCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("BLCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="BLCAmethylation_data_ALLprobes.rds")

#21 COAD 

#head -n 1 COAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > COAD_patients.txt

brca_pats = read.table("COAD_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("COAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="COADmethylation_data_lncs_cands_ALLprobes.rds")

#22 UCEC

#head -n 1 UCEC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UCEC_patients.txt

brca_pats = read.table("UCEC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("UCEC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="UCECmethylation_data_ALLprobes.rds")

#23 SARC

#head -n 1 SARC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > SARC_patients.txt

brca_pats = read.table("SARC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("SARC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="SARCmethylation_data_ALLprobes.rds")








































































