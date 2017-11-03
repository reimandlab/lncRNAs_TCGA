#2. miRNA- Gene expression 
mi = fread("x3t2m1.mature.UQ.mirna.matrix.txt")

#need to clean up patients 
#row 150 for miR-142-3p

#3. linker file 
link = read.csv("PCAWGMay2016DataReleasev1_v1.1_v1.2_v1.3_v1.4release_may2016v1.4.csv")

for(i in 1:ncol(link)){
	z <- length(which(link[,i] %in% colnames(exp)))
	print(c(i, z))
}

columns_wIDs = c(10, 24, 65, 79)
z1 = which(colnames(mi) %in% link[,10])
z2 = which(colnames(mi) %in% link[,24])
z3 = which(colnames(mi) %in% link[,65])
z4 = which(colnames(mi) %in% link[,79])

z = unique(c(z1,z2,z3,z4))

mi = as.data.frame(mi)
mi = mi[,c(1,z)]

for(i in :ncol(mi)){
	z1 = which(link[,10] %in% colnames(mi)[i])
	z2 = which(link[,24] %in% colnames(mi)[i])
	z3 = which(link[,65] %in% colnames(mi)[i])
	z4 = which(link[,79] %in% colnames(mi)[i])

	z = unique(c(z1,z2,z3,z4))
	if(length(z) > 1){z = z[[1]]}
	id = link$icgc_donor_id[z]
	colnames(mi)[i] = id
}

mi = mi[150,]

#4. cancer subtype 
canc = fread("pcawg_specimen_histology_August2016_v6.tsv")

lymph = canc$icgc_donor_id[canc$histology_tier2=="Lymphoid"]
z <- which(colnames(mi) %in% lymph)
mi = mi[,c(1,z)]

mts = c("DO27833" , "DO52682" , "DO52685" , "DO27851" , "DO52672" , "DO52692", 
 "DO27813" , "DO221124", "DO52717" , "DO52675" , "DO52647" , "DO27855", 
"DO27847" , "DO222308" ,"DO222305" ,"DO27835" , "DO52679")

