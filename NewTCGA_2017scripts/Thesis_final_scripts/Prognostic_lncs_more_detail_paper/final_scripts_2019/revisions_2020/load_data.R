source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_packages.R")
source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_functions.R")

#1. lncRNA expression in different cancers
rna = readRDS("lncRNAs_2019_manuscript/5919_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds") #<- updated gbm cohort

#2. Protein coding genes expression in different cancers
pcg = readRDS("lncRNAs_2019_manuscript/19438_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds")

#3. Fantom data
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
colnames(fantom)[1] = "gene"

###---------------------------------------------------------------
###Get lncRNAs for each cancer
###---------------------------------------------------------------

#Add cancer type
canc_conversion = readRDS("lncRNAs_2019_manuscript/tcga_id_cancer_type_conversion.txt")
canc_conversion = as.data.frame(canc_conversion)
canc_conversion = canc_conversion[,c(2,4)]
colnames(canc_conversion)[2] = "patient"

rna = merge(rna, canc_conversion, by="patient")
pcg = merge(pcg, canc_conversion, by="patient")

cancers = as.list(unique(rna$Cancer)) #18 cancers wtih at least 90 patients in each cohort

#Remove any lncRNAs that are not expressed in any of the patients
sums = apply(rna[,2:(ncol(rna)-34)], 2, sum) #134 with 0 expression in ALL patients
zeroes = names(sums[which(sums ==0)]) #what are they?
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

rna = as.data.table(rna)

dim(rna)
dim(pcg)

#clean up data 

#1. remove discrepancy
z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]
z = which(is.na(as.numeric(rna$age_at_initial_pathologic_diagnosis)))
rna = rna[-z,]
z = which(is.na(as.numeric(rna$OS.time)))
rna = rna[-z,]
z = which(as.numeric(rna$OS.time) == 0)
rna = rna[-z,]

#2. list of cancers to apply function to
cancers = as.list(unique(rna$Cancer))

#3. function that splits data into cancers
get_canc = function(canc){
        canc_data = rna[which(rna$Cancer == canc),]
        return(canc_data)
}

canc_datas = llply(cancers, get_canc)

canc_conv = unique(rna[,c("type", "Cancer")])

z = which(pcg$patient %in% rna$patient)
pcg = pcg[z,]

cols = colnames(rna)[which(colnames(rna) %in% colnames(pcg))]
all = merge(rna, pcg, by = cols)

t = as.data.table(table(all$type))
t = t[order(N)]
print(t)

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

#remove cancer types with less than 50 patients
pats_num = as.data.table(table(all$type))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
all = all[-which(all$type %in% canc_rm),]

print(table(all$type))


#fix clinical variables so only one level per variable has NA
z = which(rna$race %in% c("[Not Available]", "[Not Evaluated]", "[Unknown]"))
rna$race[z] = "unknown"

z = which(rna$clinical_stage %in% c("[Not Applicable]", "[Not Available]"))
rna$clinical_stage[z] = "unknown"

z = which(rna$histological_grade %in% c("[Unknown]", "[Not Available]", "[Discrepancy]"))
rna$histological_grade[z] = "unknown"

print(table(rna$race))

#fix clinical variables so only one level per variable has NA
z = which(pcg$race %in% c("[Not Available]", "[Not Evaluated]", "[Unknown]"))
pcg$race[z] = "unknown"

z = which(pcg$clinical_stage %in% c("[Not Applicable]", "[Not Available]"))
pcg$clinical_stage[z] = "unknown"

z = which(pcg$histological_grade %in% c("[Unknown]", "[Not Available]", "[Discrepancy]"))
pcg$histological_grade[z] = "unknown"

print(table(pcg$race))

print("done loading everything YAY")
