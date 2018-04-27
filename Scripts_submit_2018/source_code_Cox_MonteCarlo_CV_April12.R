###source_code_Cox_MonteCarlo_CV_Jan12.R

###Purpose--------------------------------------------------------------------

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. lncRNA expression in different cancers 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")

#2. 
lncrnas = colnames(rna)[1:(ncol(rna)-33)]

#4. Fantom data 
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
pcg = readRDS("19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")

#4. Matched normals expression 
norm = readRDS("all_genes_matched_normals_563_March13_wclinical_data.rds")

met = readRDS("all_genes_matched_metastatic_300_March13_wclinical_data.rds")


###---------------------------------------------------------------
###Get lncRNAs for each cancer 
###---------------------------------------------------------------

#Add cancer type
canc_conversion = readRDS("tcga_id_cancer_type_conversion.txt")
canc_conversion = as.data.frame(canc_conversion)
canc_conversion = canc_conversion[,c(2,4)]
colnames(canc_conversion)[2] = "patient"

norm_conversion = readRDS("tcga_id_NORMAL_samples_type_conversion.txt")
norm_conversion = as.data.frame(norm_conversion)
norm_conversion = norm_conversion[,c(2,4)]
colnames(norm_conversion)[2] = "patient"

met_conversion = readRDS("tcga_id_Metastatic_samples_type_conversion.txt")
met_conversion = as.data.frame(met_conversion)
met_conversion = met_conversion[,c(2,4)]
colnames(met_conversion)[2] = "patient"

rna = merge(rna, canc_conversion, by="patient")
pcg = merge(pcg, canc_conversion, by="patient")
norm = merge(norm, norm_conversion, by="patient")
met = merge(met, met_conversion, by="patient")

cancers = as.list(unique(rna$Cancer)) #18 cancers wtih at least 90 patients in each cohort

#3. Remove any lncRNAs that are not expressed in any of the patients 
sums = apply(rna[,2:(ncol(rna)-34)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

#4. Now within each cancer get mean and variance for each gene 
rna = as.data.table(rna)

###---------------------------------------------------------------
###Remove PCGs without expression 
###---------------------------------------------------------------

sums = apply(pcg[,2:(ncol(pcg)-34)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(pcg) %in% zeroes)
pcg = pcg[,-z]

#4. Now within each cancer get mean and variance for each gene 
pcg = as.data.table(pcg)

###---------------------------------------------------------------
###Remove all genes without expression in normal samples 
###---------------------------------------------------------------

sums = apply(norm[,2:(ncol(norm)-34)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(norm) %in% zeroes)
norm = norm[,-z]

#4. Now within each cancer get mean and variance for each gene 
norm = as.data.table(norm)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
