###methylation files 
###identify probes overlapping lncRNAs 

library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#probes overlapping lncRNAs - all we need to get out of the methylation files for now
probes = fread("fantom_lncrnas_mapped_to450_probes.bed")
colnames(probes) = c("probe_chr", "cpg_start", "cpg_end", "cgid", "rand", "cpgstrand", "lncchr", "lncstart", "lncend", "rand2", "ensg", "lncstrand", "lncname")
#print just first column to intersect with the actual methylation files
probes = probes[,4]
probes = unique(probes$cgid)
write.table(probes, file="unique_cgprobes_mapping_to_lncRNA_cands.txt", quote=F, row.names=F, col.names=F)

#############in command line###################################################################################################################################
#1
head -n 1 OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > OV_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > OV_lnc_probes.txt

#1in R
ov_pats = read.table("OV_patients.txt")
ov_pats = ov_pats[,3:ncol(ov_pats)]
ov_pats = as.character(unlist(ov_pats))
ov_meth = fread("OV_lnc_probes.txt")
colnames(ov_meth) = c("probe", ov_pats)
saveRDS(ov_meth, file="OV_methylation_data_lncs_cands.rds")

#2
head -n 1 PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > PAAD_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > PAAD_lnc_probes.txt

#2in R
paad_pats = read.table("PAAD_patients.txt")
paad_pats = paad_pats[,3:ncol(paad_pats)]
paad_pats = as.character(unlist(paad_pats))
paad_meth = fread("PAAD_lnc_probes.txt")
colnames(paad_meth) = c("probe", paad_pats)
saveRDS(paad_meth, file="PAAD_methylation_data_lncs_cands.rds")

#3
head -n 1 KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRC_lnc_probes.txt

#in R
kirc_pats = read.table("KIRC_patients.txt")
kirc_pats = kirc_pats[,3:ncol(kirc_pats)]
kirc_pats = as.character(unlist(kirc_pats))
kirc_meth = fread("KIRC_lnc_probes.txt")
colnames(kirc_meth) = c("probe", kirc_pats)
saveRDS(kirc_meth, file="KIRC_methylation_data_lncs_cands.rds")

#4
head -n 1 LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LIHC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LIHC_lnc_probes.txt

lihc_pats = read.table("LIHC_patients.txt")
lihc_pats = lihc_pats[,3:ncol(lihc_pats)]
lihc_pats = as.character(unlist(lihc_pats))
lihc_meth = fread("LIHC_lnc_probes.txt")
colnames(lihc_meth) = c("probe", lihc_pats)
saveRDS(lihc_meth, file="LIHC_methylation_data_lncs_cands.rds")

#lung adeno
#LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt

#5
head -n 1 LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUAD_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LUAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUAD_lnc_probes.txt

luad_pats = read.table("LUAD_patients.txt")
luad_pats = luad_pats[,3:ncol(luad_pats)]
luad_pats = as.character(unlist(luad_pats))
luad_meth = fread("LUAD_lnc_probes.txt")
colnames(luad_meth) = c("probe", luad_pats)
saveRDS(luad_meth, file="LUADmethylation_data_lncs_cands.rds")

#breast
#BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt

#6
head -n 1 BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BRCA_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BRCA_lnc_probes.txt

brca_pats = read.table("BRCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("BRCA_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="BRCAmethylation_data_lncs_cands.rds")

#7ACC

head -n 1 ACC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ACC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt ACC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ACC_lnc_probes.txt

brca_pats = read.table("ACC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("ACC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="ACCmethylation_data_lncs_cands.rds")

#8 ESCA

head -n 1 ESCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ESCA_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt ESCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ESCA_lnc_probes.txt

brca_pats = read.table("ESCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("ESCA_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="ESCAmethylation_data_lncs_cands.rds")

#9 GBM 

head -n 1 GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > GBM_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > GBM_lnc_probes.txt

brca_pats = read.table("GBM_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("GBM_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="GBMmethylation_data_lncs_cands.rds")

#10 HNSC

head -n 1 HNSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > HNSC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt HNSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > HNSC_lnc_probes.txt

brca_pats = read.table("HNSC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("HNSC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="HNSCmethylation_data_lncs_cands.rds")

#11 LGG 

head -n 1 LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LGG_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LGG_lnc_probes.txt

brca_pats = read.table("LGG_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("LGG_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="LGGmethylation_data_lncs_cands.rds")

#12 MESO 

head -n 1 MESO.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > MESO_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt MESO.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > MESO_lnc_probes.txt

brca_pats = read.table("MESO_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("MESO_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="MESOmethylation_data_lncs_cands.rds")

#13 READ

head -n 1 READ.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > READ_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt READ.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > READ_lnc_probes.txt

brca_pats = read.table("READ_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("READ_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="READmethylation_data_lncs_cands.rds")

#14 THCA 

head -n 1 THCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > THCA_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt THCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > THCA_lnc_probes.txt

brca_pats = read.table("THCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("THCA_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="THCAmethylation_data_lncs_cands.rds")

#15 UVM 

head -n 1 UVM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UVM_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt UVM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UVM_lnc_probes.txt

brca_pats = read.table("UVM_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("UVM_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="UVMmethylation_data_lncs_cands.rds")

#16 CESC

head -n 1 CESC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > CESC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt CESC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > CESC_lnc_probes.txt

brca_pats = read.table("CESC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("CESC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="CESCmethylation_data_lncs_cands.rds")

#17 KIRP

head -n 1 KIRP.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRP_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt KIRP.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRP_lnc_probes.txt

brca_pats = read.table("KIRP_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("KIRP_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="KIRPmethylation_data_lncs_cands.rds")

#18 STAD

head -n 1 STAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > STAD_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt STAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > STAD_lnc_probes.txt

brca_pats = read.table("STAD_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("STAD_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="STADmethylation_data_lncs_cands.rds")

#19 LUSC

head -n 1 LUSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUSC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LUSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LUSC_lnc_probes.txt

brca_pats = read.table("LUSC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("LUSC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="LUSCmethylation_data_lncs_cands.rds")

#20 BLCA 

head -n 1 BLCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BLCA_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt BLCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > BLCA_lnc_probes.txt

brca_pats = read.table("BLCA_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("BLCA_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="BLCAmethylation_data_lncs_cands.rds")

#21 COAD 

head -n 1 COAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > COAD_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt COAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > COAD_lnc_probes.txt

brca_pats = read.table("COAD_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("COAD_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="COADmethylation_data_lncs_cands.rds")

#22 UCEC

head -n 1 UCEC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UCEC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt UCEC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > UCEC_lnc_probes.txt

brca_pats = read.table("UCEC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("UCEC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="UCECmethylation_data_lncs_cands.rds")

#23 SARC

head -n 1 SARC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > SARC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt SARC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > SARC_lnc_probes.txt

brca_pats = read.table("SARC_patients.txt")
brca_pats = brca_pats[,3:ncol(brca_pats)]
brca_pats = as.character(unlist(brca_pats))
brca_meth = fread("SARC_lnc_probes.txt")
colnames(brca_meth) = c("probe", brca_pats)
saveRDS(brca_meth, file="SARCmethylation_data_lncs_cands.rds")








































































