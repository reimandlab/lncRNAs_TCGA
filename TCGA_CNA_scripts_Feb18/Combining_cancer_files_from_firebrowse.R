###Combining_cancer_files_from_firebrowse.R
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#1. LIHC
lihc_cna = fread("LIHC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
lihc_cna$canc = "liver"
lihc_cna$Chromosome = paste("chr", lihc_cna$Chromosome, sep="")

#2. KIRC
kirc_cna = fread("KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
kirc_cna$canc = "kidney"
kirc_cna$Chromosome = paste("chr", kirc_cna$Chromosome, sep="")

#3. PAAD
paad_cna = fread("PAAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
paad_cna$canc = "pancreas"
paad_cna$Chromosome = paste("chr", paad_cna$Chromosome, sep="")

#4. OV 
ov_cna = fread("OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
ov_cna$canc = "ovary"
ov_cna$Chromosome = paste("chr", ov_cna$Chromosome, sep="")

#5. BRCA
brca_cna = fread("BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
brca_cna$canc = "brca"
brca_cna$Chromosome = paste("chr", brca_cna$Chromosome, sep="")

allcnas = rbind(lihc_cna, kirc_cna, paad_cna, ov_cna, brca_cna)

#6. LUAD
luad_cna = fread("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
luad_cna$canc = "luad"
luad_cna$Chromosome = paste("chr", luad_cna$Chromosome, sep="")

allcnas = rbind(lihc_cna, kirc_cna, paad_cna, ov_cna, brca_cna, luad_cna)

saveRDS(allcnas, file="6cancers_types_7candidates_May9th_COPYNUMBER_CNAs.rds")