#TCGA-Assembler_RNASeqDownload.R

source("./Module_A.R")
source("./Module_B.R")

#Ovary
sCancer <- "OV"
sPath1 <- "./QuickStartExample/Part1_DownloadedData"
sPath2 <- "./QuickStartExample/Part2_BasicDataProcessingResult"
sPath3 <- "./QuickStartExample/Part3_AdvancedDataProcessingResult"

path_geneExp <- DownloadRNASeqData(cancerType = sCancer, assayPlatform = "gene.normalized_RNAseq", saveFolderName = sPath1)


list_geneExp <- ProcessRNASeqData(inputFilePath = path_geneExp[1], outputFileName = paste(sCancer, "geneExp", sep = "__"), 
	dataType = "geneExp", outputFileFolder = sPath2)

#Pancreas




#Breast





#3 Kidney types 





#Liver 
