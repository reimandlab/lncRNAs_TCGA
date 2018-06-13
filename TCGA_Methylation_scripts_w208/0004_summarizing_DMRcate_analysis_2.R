options(stringsAsFactors=F)
#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(plyr)

# load packages required for analysis
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(GenomicRanges)


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
colnames(fantom)[1] = "lncRNA"

mypal = pal_npg("nrc", alpha = 0.7)(10)

#MAIN PACKAGE
library(DMRcate)

#all Methylation files - already processed (NAs and 0 probes removed)
#file converted into matrix 

all_results = list.files("files_generated_by_DMR_analysis_June8/", pattern="DMRs_June8.rds") #only cancer type that should be 
#missing is OV


#keep track of candidates
#Data---
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
cands = filter(cands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, 166 unique lncRNAs, 23 unique cancer types 

#which cancer types are the non-unique lncRNAs from?
cands$Combo = NULL
cands = cands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
cands = cands[!duplicated(cands), ]
cands_dups = unique(cands$gene[which(duplicated(cands$gene))])

#1. read each file and summarize number of DMRs per lncRNA 

read_dmr_results = function(filee){
	results = readRDS(paste("files_generated_by_DMR_analysis_June8/", filee, sep=""))
	numlncs = length(results)
	print(filee)

	read_each_lnc_dmr = function(lnc){
		num_dmr = dim(lnc)[1]
		if(!is.null(num_dmr)){
		print(lnc$lnc[1])
		#get num of DMRs 
		lnc_name = lnc$lnc[1]
		canc = lnc$canc[1]
		row = c(lnc_name, canc, num_dmr)
		names(row) = c("lncRNA", "cancer", "num_sig_DMRs")
		return(row)
		}
	}

	canc_results = llply(results, read_each_lnc_dmr)
	return(canc_results)
}

all_cancers_summary_dmrs = llply(all_results, read_dmr_results)

#remove NULLs (need to figure out later which lncRNAs were evaluated that didn't have any DMRs I guess would be all lncRNAs within cancer)

rm_nulls = function(res){
	print(length(res))
	res = Filter(Negate(is.null), res)
	if(!(length(res)==0)){
	res = do.call(rbind.data.frame, res)
	colnames(res) = c("lncRNA", "cancer", "num_sig_DMRs")
	print("success")
	return(res)
}
}

all_cancers_dmrs = llply(all_cancers_summary_dmrs, rm_nulls)
all_cancers_dmrs = Filter(Negate(is.null), all_cancers_dmrs)
all_cancers_dmrs = do.call(rbind.data.frame, all_cancers_dmrs)
all_cancers_dmrs$num_sig_DMRs = as.numeric(all_cancers_dmrs$num_sig_DMRs)

all_cancers_dmrs = as.data.table(all_cancers_dmrs)
all_cancers_dmrs = merge(all_cancers_dmrs, fantom, by="lncRNA")
all_cancers_dmrs = all_cancers_dmrs[order(num_sig_DMRs)]


#plot barplot to summarize 
pdf("DMR_analysis_summary_dmrs_per_cancertype.pdf", width=12, height=6)
g = ggbarplot(all_cancers_dmrs, x= "CAT_geneName", y="num_sig_DMRs", fill="cancer", 
	label = TRUE, label.pos = "in", lab.size=1.5) + theme_light() + ggtitle("Num DMRs per lncRNA") + xlab("lncRNA") + ylab("num FDR < 0.05 DMRs")
ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)
dev.off()


#-------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------PART2-----------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------

#2. get coordinates of DMRs and see how many lncRNAs fall into DMRs -- for those ones how many DMRs/lncRNA that it falls into

lnc_coords = fread("lncrna_coords_bed.bed")

read_dmr_results = function(filee){
	results = readRDS(paste("files_generated_by_DMR_analysis_June8/", filee, sep=""))
	numlncs = length(results)
	print(filee)

	read_each_lnc_dmr = function(lnc){
		num_dmr = dim(lnc)[1]
		if(!is.null(num_dmr)){
			#get coordinates of lncRNA
			lncrna = lnc$lnc[1]
			z = which(lnc_coords$V5 == lncrna)

			if(!(length(z) ==0)){

			lncrna_cord = lnc_coords[z,]
			lncrna_cord$V6 = "*"

			#look at only lncRNA promoter 
			lnc_prom_start = lncrna_cord$V2-2000
			lnc_prom_end = lncrna_cord$V2+2000

			lncrna_cord$V2 = lnc_prom_start
			lncrna_cord$V3 = lnc_prom_end

			#make granges object for lncRNA coord
			gr <- GRanges(
		     seqnames = Rle(lncrna_cord$V1, seq(1, length(lncrna_cord$V1))),
    		 ranges = IRanges(lncrna_cord$V2, end = lncrna_cord$V3, names = lncrna_cord$V5),
    		 strand = Rle(strand(lncrna_cord$V6), length(lncrna_cord$V1)),
    		 score = 1:length(lncrna_cord$V1))

			#make granges object for DMRs
			lnc$seqnames = as.character(lnc$seqnames)
			lnc$strand = as.character(lnc$strand)

			dmrs = makeGRangesFromDataFrame(lnc)

			#intersect with DMR coordinates
			overlap = as.data.frame(intersect(gr, dmrs))
			if(!(dim(overlap)[1]) == 0){
			#print(dim(as.data.frame(overlap))[1])
			overlap$lnc = lncrna
			overlap$canc = lnc$canc[1]
			overlap$num_DMRs_per_lncRNA = num_dmr
			return(overlap)
			}

			if((dim(overlap)[1]) == 0){
			overlap = c("no", "no", "no", "no", "no", lncrna, lnc$canc[1], num_dmr)
			names(overlap) = c("seqnames", "start", "end", "width", "strand", "lnc", "canc", "num_DMRs_per_lncRNA")
			return(overlap)
			}

		}
		}
	}

	canc_results = llply(results, read_each_lnc_dmr)
	canc_results = Filter(Negate(is.null), canc_results)
	if(!(length(canc_results))==0){
	canc_results = do.call(rbind.data.frame, canc_results)
	colnames(canc_results) = c("seqnames", "start", "end", "width", "strand", "lnc", "canc", "num_DMRs_per_lncRNA")
	return(canc_results)
	}
}


lncrna_dmr_overlap = llply(all_results, read_dmr_results)
lncrna_dmr_overlap = do.call(rbind.data.frame, lncrna_dmr_overlap)
colnames(lncrna_dmr_overlap) = c("seqnames", "start", "end", "width", "strand", "lnc", "canc", "num_DMRs_per_lncRNA")










