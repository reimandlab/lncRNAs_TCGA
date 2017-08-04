#top5_cancers_extraction_script3B.R

#Karina Isaev
#August 1st, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script3B - using the top 5 cancer types chosen
#PLOTS:

#1. Plot cancer specific lncRNA expression among the different cancer types
#with wilcoxon test to see if statistical significantly means between the different 
#cancers relative to the cancer that it's specific to 

#NOTE: CONTINUATION OF SCRIPT3 


#---------------------------------------------------------
#Analysis - how many/which lncRNAs are specific to each
#cancer type?
#---------------------------------------------------------

f <- fread("all_lncs_cancers.txt", data.table=F, sep="_")
#remove individual file headers
z <- which(f[,1] %in% "gene")
f <- f[-z,]

#subset to include only lncs covered by FANTOM
z <- which(f[,1] %in% fantom[,1])
f <- f[z,]

#remove 7SK gene from the list as its not cancer unique 
z <- which(f[,1] %in% fantom[which(fantom$CAT_geneName=="7SK"),1])
f <- f[-z,] #end up with 86 unique genes in the list

#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1]))
counts <- counts[order(N)]
#50/86 lncRNAs appear only once
#9/86 lncRNAs appear in all cancer types 

#unique to a single cancer type
specific <- counts$V1[counts$N ==1]
specific <- as.data.frame(specific)
#gene ids conversion ENSG-Hugo
specific$gene <- ""
for(i in 1:nrow(specific)){
	g <- specific[i,1]
	z <- which(fantom$CAT_geneID %in% g)
	specific$gene[i] <- fantom$CAT_geneName[z]
}

#common to all cancer types
all <- counts$V1[counts$N ==7]
all <- as.data.frame(all)
#gene ids conversion ENSG-Hugo
all$gene <- ""
for(i in 1:nrow(all)){
	g <- all[i,1]
	z <- which(fantom$CAT_geneID %in% g)
	all$gene[i] <- fantom$CAT_geneName[z]
}


#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1], f[,2]))
counts <- counts[order(N)]

#remove 0s 
counts <- counts[!(counts$N==0)]

#want to get list of genes that occur in only one cancer type what are they?
once <- as.data.table(table(f[,1]))
once <- once[order(N)]
once <- once[once$N==1] #50 unique genes that appear in one cancer 

counts <- counts[which(counts$V1 %in% once$V1), ]
counts <- as.data.frame(counts)

num_cancers <-  as.data.table(table(counts$V2))
colnames(num_cancers)[1] <- "new"

tum_types$new <- ""
for(i in 1:nrow(tum_types)){
	t1 <- tum_types$V1[i]
	t2 <- tum_types$V2[i]
	n <- paste(t1, t2, sep=" ")
	tum_types$new[i] <- n
}

num_cancers <- merge(tum_types, num_cancers, by="new")
num_cancers <- num_cancers[,-(2:3)]
num_cancers$rows <- ""
for(i in 1:nrow(num_cancers)){
	t1 <- num_cancers$V3[i]
	t2 <- num_cancers$N[i]
	n <- t1*t2
	num_cancers$rows[i] <- n
}

counts$num_patients <- ""
for(i in 1:nrow(counts)){
	t <- counts[i,2]
	z <- which(num_cancers$new %in% t)
	counts$num_patients[i] <- num_cancers$V3[z]
}

#---------------------------------------------------------
##Set up matrix with gene expression values only for patients 
##from the cancer for which the gene's expression is specific to 
#---------------------------------------------------------

rows <- sum(as.numeric(num_cancers$rows))
to_plot <- as.data.frame(matrix(ncol=4, nrow=rows))
colnames(to_plot) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(counts)){
	gene2 <- counts[i,1]
	tis <- counts$V2[i]
	if(i == 1){
		a <- as.numeric(counts$num_patients[i])
		to_plot[1:a,1] <- gene2
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]
		z <- which(dat$canc %in% tis)	
		dat <- dat[z,]
		to_plot[1:a,2] <- dat$canc
		to_plot[1:a,3] <- rownames(dat)
		to_plot[1:a,4] <- dat[,1]
	}
	if(!(i==1)){
		a <- as.numeric(counts$num_patients[i])-1
		z <- which(is.na(to_plot[,1]))[1]
		to_plot[z:(z+a), 1] <- gene2

		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		
		z3 <- which(dat$canc %in% tis)	
		dat <- dat[z3,]

		to_plot[z:(z+a),2] <- dat$canc
		to_plot[z:(z+a),3] <- rownames(dat)
		to_plot[z:(z+a),4] <- dat[,1]
	}
}

to_plot$Cancer <- as.factor(to_plot$Cancer)

##Plot violin
to_plot$gene2 <- ""
for(i in 1:nrow(to_plot)){
	g <- to_plot$Gene[i]	 
	z <- which(fantom$CAT_geneID %in% g)
	to_plot$gene2[i] <- fantom$CAT_geneName[z]
}
to_plot_logged <- to_plot  
to_plot_logged$GeneE <- log1p(to_plot_logged$GeneE)


#---------------------------------------------------------
##Set up matrix to compare the expression of the cancer specific lncRNAs
##across the different cancer types 
#---------------------------------------------------------

#num of rows = num_genes(50) * patients(497) = 24850
colnames(specific)[2] <- "gene2"

specific_genes <- as.data.frame(matrix(nrow=24850, ncol=4))
colnames(specific_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(specific)){
	gene2 <- specific[i,1]
	gene <- specific[i,2]
	if(i == 1){
		specific_genes[1:497,1] <- gene
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]	
		specific_genes[1:497,2] <- dat$canc
		specific_genes[1:497,3] <- rownames(dat)
		specific_genes[1:497,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(specific_genes[,1]))[1]
		specific_genes[z:(z+496), 1] <- gene
		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		specific_genes[z:(z+496),2] <- dat$canc
		specific_genes[z:(z+496),3] <- rownames(dat)
		specific_genes[z:(z+496),4] <- dat[,1]
	}
}

##Plot violin
specific_genes_logged <- specific_genes  
specific_genes_logged$GeneE <- log1p(specific_genes_logged$GeneE)

#GOAL: For each cancer specific lncRNA, plot the expression distribution in 
#each cancer type & using wilcoxon test compare the means of cancer-reference 
#so 6 tests in each gene plot 

gene_canc <- which(duplicated(to_plot_logged$Gene))
gene_canc <- to_plot_logged[-gene_canc, ]

specific_genes_logged$ref <- ""
for(i in 1:nrow(specific_genes_logged)){
	g <- specific_genes_logged$Gene[i]
	z <- which(gene_canc$gene2 %in% g)
	specific_genes_logged$ref[i] <- as.character(gene_canc$Cancer[z])
}

##***********************************
#Plot each cancer type seperatley ... 

for(i in 1:length(unique(specific_genes_logged$ref))){
	data <- subset(specific_genes_logged, specific_genes_logged$ref %in% unique(specific_genes_logged$ref)[i])	
	#make plot - reference group is the cancer type 
	ref <- unique(specific_genes_logged$ref)[i]
	my_comparisons <- list()
	k <- 1
	t <- ref

	for(j in 1:length(unique(specific_genes_logged$Cancer))){
		new <- unique(specific_genes_logged$ref)[j]	
		if(!(new==t)){
		add <- c(t, new)
		my_comparisons[[k]] <- add 
		k <- k+1
	}
	}

	data$Gene <- as.factor(data$Gene)
		file <- paste(ref, "cancer_specific_expressionamongOthers.pdf", sep="_")
		pdf(file, pointsize=8, width=14, height=13)
	
		for (y in seq(1, length(unique(data$Gene)), 4)) {
    
    	g <- ggviolin(data[data$Gene %in% levels(data$Gene)[y:(y+3)], ], 
                  x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene",
                  order = unique(data$Cancer)[c(1,3,7,2,4,5,6)])
   		g <- g + stat_compare_means(comparisons = my_comparisons)
    	g <- ggpar(g, font.legend = c(6, "plain", "black")) 
		g <- g + labs(title = ref, y="log1p(FPKM)") + 
     	 theme(plot.title = element_text(hjust = 0.5))
		g <- g + rremove("x.text")
		g <- g + rremove("x.ticks")  
		g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
		print(g)
		}
	
	    dev.off()

}

###Using the data for the 50 lncRNA expressed in a tissue specific manner 
###save it and use it for survival analysis and co-expression analysis ... 

write.table(specific_genes_logged, file="logged_geneExpression_50unique_specific_lncRNAs.txt", quote=F, row.names=F, sep=";")

#make the same type of dataframe for the 8 lncRNAs expressed in all cancers 

#---------------------------------------------------------------
##Set up matrix to compare the expression of the common lncRNAs
##across the different cancer types 
#---------------------------------------------------------------

#num of rows = num_genes(9) * patients(497) = 24850
colnames(all)[2] <- "gene2"

all_genes <- as.data.frame(matrix(nrow=4473, ncol=4))
colnames(all_genes) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(all)){
	gene2 <- all[i,1]
	gene <- all[i,2]
	if(i == 1){
		all_genes[1:497,1] <- gene
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]	
		all_genes[1:497,2] <- dat$canc
		all_genes[1:497,3] <- rownames(dat)
		all_genes[1:497,4] <- dat[,1]
	}
	if(!(i==1)){
		z <- which(is.na(all_genes[,1]))[1]
		all_genes[z:(z+496), 1] <- gene
		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		all_genes[z:(z+496),2] <- dat$canc
		all_genes[z:(z+496),3] <- rownames(dat)
		all_genes[z:(z+496),4] <- dat[,1]
	}
}

##Plot violin
all_genes_logged <- all_genes  
all_genes_logged$GeneE <- log1p(all_genes_logged$GeneE)

write.table(all_genes_logged, file="logged_geneExpression_9unique_common_lncRNAs.txt", quote=F, row.names=F, sep=";")















