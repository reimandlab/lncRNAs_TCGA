#4. perform 1000CV survival LASSO on each cancer 

source("universal_LASSO_survival_script.R")

set.seed(911)

#lasso_cv = function(dat){
for(p in 6:10){
  dat = canc_datas[[p]]
  genes_results = list()
  cinds_clin = c()
  cinds_justlncs = c()
  cinds_combined = c()
  set.seed(101) 
  
  stage = as.data.table(table(dat$clinical_stage))
  rm = stage$V1[which(stage$N <10)]
  z = which(dat$clinical_stage %in% rm)
  if(!(length(z)==0)){
    dat = dat[-z,]}
  
  grade = as.data.table(table(dat$histological_grade))
  rm = grade$V1[which(grade$N <10)]
  z = which(dat$histological_grade %in% rm)
  if(!(length(z)==0)){
    dat = dat[-z,]}
  
  for(j in 1:100){
    
    smp_size <- floor(0.7 * nrow(dat))
    train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
    train <- dat[train_ind, ]
    test <- dat[-train_ind, ]
    
    train$OS.time = as.numeric(train$OS.time)
    z = which(is.na(train$OS.time))
    if(!(length(z)==0)){
      train = train[-z,]} 
    
    test$OS.time = as.numeric(test$OS.time)
    z = which(is.na(test$OS.time))
    if(!(length(z)==0)){
      test = test[-z,]} 
    
    ###---------------------------------------------------------------
    ###Extract significant predictors in training data 	
    ###---------------------------------------------------------------
    
    #use training medians to split the test set 
    medians = apply(train[,2:(ncol(train)-34)], 2, median)
    for(k in 1:length(medians)){
      med = medians[k]
      k = k +1
      if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(train[,k] > 0)
        l2 = which(train[,k] ==0)
        train[l1,k] = 1
        train[l2, k] = 0
      }
      
      if(!(med ==0)){
        l1 = which(train[,k] >= med)
        l2 = which(train[,k] < med)
        train[l1,k] = 1
        train[l2, k] = 0
      }
    }  
    
    #survival script
    #source("survival_script_march14_already_labelled_highlow.R")
    sums = apply(train[,2:(ncol(train)-34)], 2, sum) #134 with 0 expression in ALL patients 
    zeroes = names(sums[which(sums <10)]) #if less than 50, then unbalanced 
    z <- which(colnames(train) %in% zeroes)
    if(!(length(z)==0)){
      train = train[,-z]
    }
    
    print("p2")
    print("pass2")
    
    #-----------------------------------
    #TRAINING using just lncRNAs 
    #-----------------------------------
    
    z = which(str_detect(colnames(train), "ENSG"))
    if(!(length(z)) == 0){
    gene_data = train[,z]
    #gene_data = as.data.frame(train[,1:(ncol(train)-35)])
    genes = as.list(colnames(gene_data))
    x <- model.matrix( ~., gene_data)
    train$OS = as.numeric(train$OS)
    y <- Surv(train$OS.time, train$OS)

    if(table(train$OS)[2] >=10){

    cvfit = cv.glmnet(x, y, family = "cox", alpha =0.5) #uses cross validation to select
    #the best lamda and then use lambda to see which features remain in model 
    cvfit$lambda.min #left vertical line
    cvfit$lambda.1se #right vertical line 
    #active covariates and their coefficients 
    coef.min = coef(cvfit, s = "lambda.min") 
    active.min = which(coef.min != 0)
    index.min = coef.min[active.min]
    genes_keep = rownames(coef.min)[active.min]
    print(genes_keep)
    genes_results[j] = list(as.list(genes_keep))
    trainlncs = train[,c(which(colnames(train) %in% c(genes_keep,"OS", "OS.time")))]
    justlncs = try(coxph(Surv(OS.time, OS)  ~ ., data = trainlncs), TRUE)
    if(!(inherits(justlncs,"try-error"))){
    justlncs_cox_model = justlncs
    
    print("pass3")
    
    #-----------------------------------
    #TESTING using just lncRNAs 
    #-----------------------------------
    z = which(names(medians) %in% colnames(trainlncs))
    if(!(length(z)==0)){
      medians = medians[z]
      testlncs = test[,c(1,which(colnames(test) %in% colnames(trainlncs)))]
      rownames(testlncs) = testlncs$patient
      testlncs$patient = NULL
      if(colnames(testlncs)[1] == "patient.1"){
        testlncs = testlncs[,-1]
      }
      #Add high/low tags to each gene 
      for(k in 1:(ncol(testlncs)-2)){
        med = medians[which(names(medians) == colnames(testlncs)[k])]
        if(med ==0){
          #if median = 0 then anyone greater than zero is 1 
          for(m in 1:nrow(testlncs)){
            genexp <- testlncs[m,k]
            if(genexp > 0){
              testlncs[m,k] <- 1
            }
            if(genexp == 0 ){
              testlncs[m,k] <- 0
            }
          } 
        }
        if(!(med ==0)){
          for(m in 1:nrow(testlncs)){
            genexp <- testlncs[m,k]
            if(genexp >= med){
              testlncs[m,k] <- 1
            }
            if(genexp < med){
              testlncs[m,k] <- 0
            }
          } 
        }
      } #end loop assigning high/low to test set 
      
      print("pass4")
      
      colss = colnames(testlncs)
      colsskeep = which(str_detect(colss, "ENSG"))
      testlncs$OS = as.numeric(testlncs$OS)    
      pred_validation = predict(justlncs_cox_model, newdata = testlncs)
      # Determine concordance
      cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
                                            surv.event=testlncs$OS, method = "noether")
      
      cindex_validation = cindex_validation$c.index
      cinds_justlncs = c(cinds_justlncs, cindex_validation)
      
      print("pass5")

      #-------------end lncRNAs as predictors-------------------------------------------
      #-----------------------------------
      #TRAINING using just clinical variables  
      #-----------------------------------
      #remove all columns that are genes 
      z = which(str_detect(colnames(train), "ENSG"))
      clin_train = as.data.frame(train[,-z])
      rownames(clin_train) = clin_train$patient
      colss = c("age_at_initial_pathologic_diagnosis", "gender", "race", "histological_grade", "ajcc_pathologic_tumor_stage", "patient")
      clin_train = clin_train[,which(colnames(clin_train) %in% colss)]
      clin_train$age_at_initial_pathologic_diagnosis = as.numeric(clin_train$age_at_initial_pathologic_diagnosis)
      test_clin = as.data.frame(test[,which(colnames(test) %in% colnames(clin_train))])
      rownames(test_clin) = test_clin$patient
      test_clin = test_clin[,which(colnames(test_clin) %in% colnames(clin_train))]
      
      clin_train$patient = NULL
      test_clin$patient = NULL
      
            print("pass6")

      #remove columns with less than 2 contrasts 
      check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
      }
      keep_train = unlist(apply(clin_train, 2, check_contrasts))
      keep_test = unlist(apply(test_clin, 2, check_contrasts))
      
      keep = keep_test[which(keep_train %in% keep_test)]
      z = which(colnames(clin_train) %in% names(keep))
      clin_train = clin_train[,z]
      
            print("pass7")

      #make sure the levels are the same 
      for(l in 1:ncol(clin_train)){
        if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis")){
          clin_train[,l] = as.factor(clin_train[,l])
        }
      }
      
      test_clin = test_clin[,which(colnames(test_clin) %in% colnames(clin_train))]
      for(l in 1:ncol(test_clin)){
        if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis")){
          test_clin[,l] <- factor(test_clin[,l], levels=levels(clin_train[,which(colnames(clin_train) %in% colnames(clin_train)[l])]))
        }
      }
      
      #make sure there isn't only one class of some sort for example only one patient wtih "not available"
      print("pass8")

      trainclin = train[,c(which(colnames(train) %in% c(colnames(clin_train),"OS", "OS.time")))]
      trainclin$age_at_initial_pathologic_diagnosis = as.numeric(trainclin$age_at_initial_pathologic_diagnosis)
      
      #make sure there aren't any columns with just one value
      rm = which(sapply(trainclin, function(x) { length(unique(x)) }) == 1)
      if(!(length(rm))==0){
        trainclin = trainclin[,-rm]
      }

      justclin = coxph(Surv(OS.time, OS)  ~ ., data = trainclin)
      test_surv = test[,c(which(colnames(test) %in% c("patient", "OS", "OS.time")))]
      test_clin$patient = rownames(test_clin)
      testclin = merge(test_clin, test_surv, by="patient")
      #testclin = test[,which(colnames(test) %in% colnames(trainclin))]
      testclin$OS = as.numeric(testclin$OS)
      testclin$OS.time = as.numeric(testclin$OS.time)
      testclin$age_at_initial_pathologic_diagnosis = as.numeric(testclin$age_at_initial_pathologic_diagnosis)
      pred_validation = predict(justclin, newdata = testclin)
      # Determine concordance
      cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
                                            surv.event=testlncs$OS, method = "noether", na.rm=TRUE)
      cindex_validation = cindex_validation$c.index
      cinds_clin = c(cinds_clin, cindex_validation)
      print(colnames(trainclin))
            print("pass9")

      #-------------end clinical variables as predictors------------------------------------
      #-----------------------------------
      #TRAINING using lncRNAs and clinical variables  
      #-----------------------------------
      
      variables = unique(c(colnames(trainclin), colnames(trainlncs)))
      trainboth = train[,which(colnames(train) %in% variables)]
      trainboth$age_at_initial_pathologic_diagnosis = as.numeric(trainboth$age_at_initial_pathologic_diagnosis)

      #check if have any nas 
      trainboth$OS.time = as.numeric(trainboth$OS.time)
      trainboth$OS = as.numeric(trainboth$OS)

      remove = c()

      perc5 = 0.05*dim(trainboth)[1]

      for(b in 1:ncol(trainboth)){
        if(!(colnames(trainboth)[b] %in% c("age_at_initial_pathologic_diagnosis", "OS", "OS.time"))){
        dg = as.data.table(table(trainboth[,b]))
        print(dg)
        z = which(dg$N <perc5)
        if(!(length(z)==0)){
          var = dg$V1[z]
          namevar = colnames(trainboth)[b]
          remove = c(remove, namevar)
        }
        }
      }

      if(!(length(remove)==0)){
      trainboth = trainboth[,-(which(colnames(trainboth) %in% remove))]
      }
      print("pass10")


      both = try(coxph(Surv(OS.time, OS)  ~ ., data = trainboth), TRUE)
      if(!(inherits(both,"try-error"))){
      #testboth = test[,which(colnames(test) %in% colnames(trainboth))]
      testlncs$patient = rownames(testlncs)
      testboth = merge(testlncs, testclin, by=c("patient", "OS", "OS.time"))
      testboth$OS = as.numeric(testboth$OS)
      testboth$OS.time = as.numeric(testboth$OS.time)
      testboth$age_at_initial_pathologic_diagnosis = as.numeric(testboth$age_at_initial_pathologic_diagnosis)
      pred_validation = predict(both, newdata = testboth)
      # Determine concordance
      cindex_validation = concordance.index(pred_validation, surv.time = testboth$OS.time,
                                            surv.event=testboth$OS, method = "noether", na.rm=TRUE)
      cindex_validation = cindex_validation$c.index
      cinds_combined = c(cinds_combined, cindex_validation)
     }
   }
    }  
    }
    } #if(!(length(z)==0))
    
}#dim(pvalues)[1] >2

#end 1:1000 loop

if(!(length(cinds_combined)==0)){
  
  cinds_combined = as.data.frame(cinds_combined)
  cinds_combined$canc = dat$Cancer[1]
  colnames(cinds_combined)[1] = "cindex"
  cinds_combined$type = "cinds_combined"
  
  cinds_clin = as.data.frame(cinds_clin)
  cinds_clin$canc = dat$Cancer[1]
  colnames(cinds_clin)[1] = "cindex"
  cinds_clin$type = "cinds_clin"
  
  cinds_justlncs = as.data.frame(cinds_justlncs)
  cinds_justlncs$canc = dat$Cancer[1]
  colnames(cinds_justlncs)[1] = "cindex"
  cinds_justlncs$type = "cinds_justlncs"
  
  all_cinds = rbind(cinds_combined, cinds_clin, cinds_justlncs)
  filename_cinds = paste(dat$Cancer[1], "cindices_100CV_May2_ElasticNet_no_prefilter.rds", sep="_")
  saveRDS(all_cinds, file=filename_cinds)
  
  filename = paste(dat$Cancer[1], "lncRNAs_selected_100CV_May2_ElasticNet_no_prefiler.rds", sep="_")
  saveRDS(genes_results, file=filename)
}
}

#canc_datas_survival_results = llply(canc_datas, lasso_cv)




