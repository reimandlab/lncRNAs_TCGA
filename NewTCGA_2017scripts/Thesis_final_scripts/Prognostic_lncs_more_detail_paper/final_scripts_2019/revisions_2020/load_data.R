source("load_packages.R")
source("load_functions.R")

#1. lncRNA expression in different cancers
rna = readRDS("lncRNAs_2019_manuscript/5919_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds") #<- updated gbm cohort

#2. Fantom data
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

#3. Protein coding genes expression in different cancers
pcg = readRDS("lncRNAs_2019_manuscript/19438_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds")

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

#3. Remove any lncRNAs that are not expressed in any of the patients
sums = apply(rna[,2:(ncol(rna)-34)], 2, sum) #134 with 0 expression in ALL patients
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

#4. Now within each cancer get mean and variance for each gene
rna = as.data.table(rna)

dim(rna)
dim(pcg)

#function that tests each lncRNA's survival

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

#lncRNA candidates, n = 166, n=173 combos
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#cindices
r = readRDS("lncs_cross_validations_Results_nov20.rds")

z = which(duplicated(all$patient))
if(!(length(z)==0)){
  dups = all$patient[z]
  k = which(all$patient %in% dups)
  all = all[-k,]
}

print(table(all$type))
