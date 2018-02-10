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
head -n 1 OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > OV_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > OV_lnc_probes.txt

#in R
ov_pats = read.table("OV_patients.txt")
ov_pats = ov_pats[,3:ncol(ov_pats)]
ov_pats = as.character(unlist(ov_pats))
ov_meth = fread("OV_lnc_probes.txt")
colnames(ov_meth) = c("probe", ov_pats)
saveRDS(ov_meth, file="OV_methylation_data_lncs_cands.rds")

head -n 1 PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > PAAD_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt PAAD.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > PAAD_lnc_probes.txt

#in R
paad_pats = read.table("PAAD_patients.txt")
paad_pats = paad_pats[,3:ncol(paad_pats)]
paad_pats = as.character(unlist(paad_pats))
paad_meth = fread("PAAD_lnc_probes.txt")
colnames(paad_meth) = c("probe", paad_pats)
saveRDS(paad_meth, file="PAAD_methylation_data_lncs_cands.rds")

head -n 1 KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt KIRC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > KIRC_lnc_probes.txt

#in R
kirc_pats = read.table("KIRC_patients.txt")
kirc_pats = kirc_pats[,3:ncol(kirc_pats)]
kirc_pats = as.character(unlist(kirc_pats))
kirc_meth = fread("KIRC_lnc_probes.txt")
colnames(kirc_meth) = c("probe", kirc_pats)
saveRDS(kirc_meth, file="KIRC_methylation_data_lncs_cands.rds")

head -n 1 LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LIHC_patients.txt
grep -f unique_cgprobes_mapping_to_lncRNA_cands.txt LIHC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > LIHC_lnc_probes.txt

lihc_pats = read.table("LIHC_patients.txt")
lihc_pats = lihc_pats[,3:ncol(lihc_pats)]
lihc_pats = as.character(unlist(lihc_pats))
lihc_meth = fread("LIHC_lnc_probes.txt")
colnames(lihc_meth) = c("probe", lihc_pats)
saveRDS(lihc_meth, file="LIHC_methylation_data_lncs_cands.rds")

