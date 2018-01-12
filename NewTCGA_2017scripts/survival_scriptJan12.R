surv_test = function(gene){
  print(gene)
  results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")
  #1. Subset lnc_rna to those patients in cancer
  df <- subset(rna, rna$canc %in% cancer)
  z <- which(colnames(df) %in% gene)
  if(!(length(z)==0)){
  df = as.data.frame(df)
  df <- df[,c(z,5786:5790)]  
  df[,1] <- log1p(df[,1])

  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  median2 <- quantile(as.numeric(df[,1]), 0.5)
  if(median2 ==0){
  	median2 = mean(as.numeric(df[,1]))
  }

  #median2 <- median(df[,1])
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
  gene <- colnames(df)[1]
  df$status[df$status=="Alive"] <- 0
  df$status[df$status=="Dead"] <- 1
  df$status <- as.numeric(df$status)
  df$time <- as.numeric(df$time)
      
  #cox regression 
  res.cox <- coxph(Surv(time, status) ~ median, data = df)
  row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)], df$canc[1])
  names(row) <- names(results_cox)
}
 return(row)
}

