setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")
print("hi")

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2020_manuscript")

print("done source script")

#survival script
#source("survival_script_march14_already_labelled_highlow.R")

library(doParallel)
library(curl)
library(RCurl)
library(selectiveInference)
args = commandArgs(trailingOnly = TRUE)
print(args[1])

#####################################################################
#READ ME
#####################################################################

#####################################################################
#THIS SCRIPT USES REAL DATA FOR EACH CANCER TYPE
#TO RUN ELASTIC NET USING 1000 CROSS-VALIDATIONS
#AND RETURN LNCRNAS THAT WERE SELECTED IN EACH ROUND
#ONCE THESE ARE ACQUIRED A FULL MODEL IS FIT
#SHOULD USE PROPER INFERENCE METRICS FOR MODELS THAT ALREADY UNDERWENT
#FEATURE SELECTION SINCE SAMPLE SPACE IS REDUCED

#RESULTS ARE SAVED IN "real_elastic_net_runs"

#####################################################################

#can we turn this into a function?

#registerDoParallel(cores=28)

#Get CANCER index as arguement
args = commandArgs(trailingOnly = TRUE)

set.seed(911)

#GET LIST OF CANCERS
cancers = as.data.table(table(rna$Cancer))
cancers = as.data.table(filter(cancers, N >=50))
cancers= cancers[order(N)]
cancers = unique(cancers$V1)
print(cancers)

print(args[1])
canc = cancers[as.integer(args[1])]
print(canc)
date = Sys.Date()

cols_keep = c("race", "clinical_stage", "histological_grade")

print(table(rna$histological_grade))
print(table(rna$clinical_stage))

#[1] ----- get cancer info ----------------------------
#---------get full cancer dataset ---------------------

get_canc = function(canc){
  canc_data = rna[which(rna$Cancer == canc),]
  print(dim(canc_data))
  return(canc_data)
}

#[2] ----- check clinical covariates ----------------------------

prepare_dat = function(dat){

  #check covariates have at least 10 patients in categories
  #stage = as.data.table(table(dat$clinical_stage))
  #rm = stage$V1[which(stage$N <10)]
  #z = which(dat$clinical_stage %in% rm)
  #if(!(length(z)==0)){
  #  dat = dat[-z,]}

  #grade = as.data.table(table(dat$histological_grade))
  #rm = grade$V1[which(grade$N <10)]
  #z = which(dat$histological_grade %in% rm)
  #if(!(length(z)==0)){
  #  dat = dat[-z,]}

  print("checked done")
  return(dat)

}

#[3] ----- elsatic net & splits ------------------------------------

main_elastic_net = function(dat){

  #set.seed(101)

  #given random dataset, split into training and test set
  smp_size <- floor(0.7 * nrow(dat))
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  train <- dat[train_ind, ]
  test <- dat[-train_ind, ]

  train$OS.time = as.numeric(train$OS.time)
  test$OS.time = as.numeric(test$OS.time)

  ###---------------------------------------------------------------
  ###Extract significant predictors in training data
  ###---------------------------------------------------------------

  z_train = which(str_detect(colnames(train), "ENSG"))
  train=as.data.frame(train)
  rownames(train) = train$patient
  train[,z_train] = apply(train[,z_train], 2, as.numeric)

  z = which(str_detect(colnames(test), "ENSG"))
  test=as.data.frame(test)
  rownames(test) = test$patient
  test[,z] = apply(test[,z], 2, as.numeric)

  #use training medians to split the test set
  medians = apply(train[,z_train], 2, median)

  #assign high/low labels to each patient for each lncRNA
  for(k in 1:length(medians)){
      med = medians[k]
      k = which(colnames(train) == names(med))
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
      } #end assigning labels to each patient for each gene

  lncs_dat = which(str_detect(colnames(dat), "ENSG"))
  sums = apply(dat[,..lncs_dat], 2, sum)

  #for the full cohort get lncRNAs that have expression in less than 15 patients
  #or less than 10 patients

  perc_10 = dim(dat)[1]*0.1

  #how many lncRNAs have only 10% of cohort or 15 patients (whichever is greater)
  #with > 0 expression in the cohort
  z1 = which(sums < perc_10)
  print(perc_10)
  z2 = which(sums <= 15)
  all_rm = names(sums)[unique(c(z1,z2))]

  z <- which(colnames(train) %in% all_rm)

  if(!(length(z)==0)){
      train = train[,-z]
    }

  z_train = which(str_detect(colnames(train), "ENSG"))
  genes = colnames(train)[z_train]
  print(length(genes))

  ###---------------------------------------------------------------
  ###Evaluate each gene using survival model
  ###---------------------------------------------------------------

  surv_test = function(gene){

  #print(gene)
  results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coefficient", "HR", "pval", "low95", "upper95")

  #1. Subset lnc_rna to those patients in cancer
  z =  which(colnames(train) %in% gene)

  if(!(length(z)==0)){

  df = as.data.frame(train)
  cols_keep = c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "OS", "OS.time")
  z2 = which(colnames(train) %in% cols_keep)

  df = df[,c(z,z2)]

  gene <- colnames(df)[1]
  df$OS <- as.numeric(df$OS)
  df$OS.time <- as.numeric(df$OS.time)
  df$age_at_initial_pathologic_diagnosis = as.numeric(df$age_at_initial_pathologic_diagnosis)

  colnames(df)[1] = "median"

  if((!(table(df$OS)[2] <5)) & (!(dim(table(df$OS)) == 1))){

  #cox regression
  res.cox <- coxph(Surv(OS.time, OS) ~ median + age_at_initial_pathologic_diagnosis, data = df)

  row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)])
  names(row) <- names(results_cox)
  return(row)

  }
  }
  }

  tests_survivaltwo = llply(genes, surv_test, .progress="text")
  print("p1")
  tests_survivalthree = as.data.table(ldply(tests_survivaltwo))
  print("p2")

  tests_survival3 = tests_survivalthree

  #as long as results are not empty, keep going...

  if(!(dim(tests_survival3)[1] ==0)){

      print("keep going")

      tests_survival3$pval = as.numeric(tests_survival3$pval)
      tests_survival3$fdr = p.adjust(tests_survival3$pval, method="fdr")

      #keep only lncRNAs that are univaraitly significant
      tests_survival3 = subset(tests_survival3, pval <=0.05)

      if(!(dim(tests_survival3)[1])==0){

        print(length(unique(tests_survival3$gene)))

        sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
        tests_survival3$HR = as.numeric(tests_survival3$HR)
        tests_survival3 = as.data.table(tests_survival3)
        sort.field = "HR"
        tests_survival3 = sortme(tests_survival3, sort.field)

        print("pass")

        tests_survival3 = as.data.table(tests_survival3)
        tests_survival3$coefficient = as.numeric(unlist(tests_survival3$coefficient))
        tests_survival3 = tests_survival3[order(-abs(coefficient))]

        z = table(train$OS)[2]
        tests_survival3 = as.data.frame(tests_survival3)

        #if number of deaths less than sig gene candidates
        #keep only the top genes, this was recommended by Nature Biotech paper
        #inspiration behind this project

        if(z < dim(tests_survival3)[1]){
          tests_survival3 = tests_survival3[1:z, ]
        }

        print("pass2")

        z1 = which(colnames(train) %in% tests_survival3$gene)
        z2 = which(!(str_detect(colnames(train), "ENSG")))

        #remake training dataset using only significant prognostic lncRNAs
        #and all remaining clinical variables

        train = train[,c(z1,z2)]

        #now make sure dataset is not empty
        #and start building actual models

        if(!(dim(tests_survival3)[1] <5)){

          #-----------------------------------
          #TRAINING using just lncRNAs
          #-----------------------------------
          z = which(str_detect(colnames(train), "ENSG"))
          gene_data = as.data.frame(train[,z])
          genes = as.list(colnames(gene_data))
          x <- model.matrix( ~., gene_data)
          train$OS = as.numeric(train$OS)
          y <- Surv(train$OS.time, train$OS)
          cvfit = cv.glmnet(x, y, family = "cox", alpha =0.5) #uses cross validation to select
          #the best lamda and then use lambda to see which features remain in model
          cvfit$lambda.min #left vertical line
          cvfit$lambda.1se #right vertical line
          #active covariates and their coefficients
          coef.min = coef(cvfit, s = "lambda.min")
          active.min = which(coef.min[,1] != 0)
          index.min = coef.min[active.min]
          genes_keep = rownames(coef.min)[active.min]

          genes_results = list(as.list(genes_keep))

          trainlncs = train[,c(which(colnames(train) %in% c(genes_keep,"OS", "OS.time")))]

          #fit model using lncRNAs selected by elastic net
          justlncs = try(coxph(Surv(OS.time, OS)  ~ ., data = trainlncs), TRUE)

          if(!(inherits(justlncs,"try-error"))){

            justlncs_cox_model = justlncs
            print("lnc based model done")

            #-----------------------------------
            #TESTING using just lncRNAs
            #-----------------------------------

            z = which(names(medians) %in% colnames(trainlncs))
            #label test set using gene labels from the training set for each gene

            if(!(length(z)==0)){

              medians = medians[z]
              testlncs = test[,c(which(colnames(test) %in% colnames(trainlncs)))]
              #rownames(testlncs) = testlncs$patient
              #testlncs$patient = NULL
              #Add high/low tags to each gene
              z = which(str_detect(colnames(testlncs), "ENSG"))
              for(k in 1:length(z)){
                k = z[k]
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
              cinds_justlncs = cindex_validation
              print(cinds_justlncs)

          #-------------end lncRNAs as predictors-------------------------------------------

          #---------------------------------------------------------------------------------
          #TRAINING using just clinical variables
          #---------------------------------------------------------------------------------

          z = which(!(str_detect(colnames(train), "ENSG")))
          clin_train = as.data.frame(train[,z])
          rownames(clin_train) = clin_train$patient
          cols_keep = c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "OS", "OS.time")
          clin_train = clin_train[,cols_keep]
          clin_train$age_at_initial_pathologic_diagnosis = as.numeric(clin_train$age_at_initial_pathologic_diagnosis)

          print(dim(clin_train))

          test_clin = test
          rownames(test_clin) = test_clin$patient
          test_clin = test_clin[,which(colnames(test_clin) %in% colnames(clin_train))]

          print(dim(test_clin))

              #remove columns with less than 2 contrasts
              check_contrasts = function(col){
                check = dim(table(col))
                if(check >1){
                  return("keep")
                }
              }
              keep_train = unlist(apply(clin_train, 2, check_contrasts))
              keep_test = unlist(apply(test_clin, 2, check_contrasts))

              keep = keep_train[which(keep_train %in% keep_test)]
              keep = c(keep)
              print(keep)

              if(length(which(names(keep) %in% "OS.time")) == 0){
                keep = c(keep, "keep")
                names(keep)[length(keep)] = "OS.time"
              }

              clin_train = as.data.frame(clin_train)
              print(head(clin_train))
              z = which(colnames(clin_train) %in% names(keep))
              print(length(z))
              clin_train = clin_train[,z]

              print(colnames(clin_train))

              clin_train$OS = as.numeric(as.character(clin_train$OS))
              clin_train$OS.time = as.numeric(as.character(clin_train$OS.time))

              print(head(clin_train))
              print(summary(clin_train$OS.time))
              print(summary(clin_train$OS))

              #make sure the levels are the same
              for(l in 1:ncol(clin_train)){
                if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis")){
                  clin_train[,l] = as.factor(clin_train[,l])
                }
              }

              print("done random for loop")

              test_clin = test_clin[,which(colnames(test_clin) %in% colnames(clin_train))]

              print(head(test_clin))

              z = which(!(str_detect(colnames(test_clin), "OS")))
              print(z)

              for(l in 1:length(z)){
                #print(colnames(clin_train)[z[l]])
                if(!(colnames(clin_train)[z[l]] == "age_at_initial_pathologic_diagnosis")){
                  test_clin[,z[l]] <- factor(test_clin[,z[l]], levels=levels(clin_train[,which(colnames(clin_train) %in% colnames(clin_train)[z[l]])]))
                }
              }

              trainclin = train[,c(which(colnames(train) %in% c(colnames(clin_train),"OS", "OS.time")))]
              trainclin$age_at_initial_pathologic_diagnosis = as.numeric(trainclin$age_at_initial_pathologic_diagnosis)
              #make sure there aren't any columns with just one value
              rm = which(sapply(trainclin, function(x) { length(unique(x)) }) == 1)
              if(!(length(rm))==0){
                trainclin = trainclin[,-rm]
              }

              justclin = try(coxph(Surv(OS.time, OS)  ~ ., data = trainclin), TRUE)
              #print(justclin)

              if(!(inherits(justclin,"try-error"))){

                testclin = test_clin
                testclin$OS = as.numeric(testclin$OS)
                testclin$OS.time = as.numeric(testclin$OS.time)
                testclin$age_at_initial_pathologic_diagnosis = as.numeric(testclin$age_at_initial_pathologic_diagnosis)
                pred_validation = predict(justclin, newdata = testclin)
                # Determine concordance
                cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
                                                      surv.event=testlncs$OS, method = "noether", na.rm=TRUE)
                cindex_validation = cindex_validation$c.index
                cinds_clin = cindex_validation

          #-------------end clinical variables as predictors------------------------------------

          #-------------------------------------------------------------------------------------
          #TRAINING using lncRNAs and clinical variables
          #-------------------------------------------------------------------------------------

          variables = unique(c(colnames(trainclin), colnames(trainlncs)))
          trainboth = train[,which(colnames(train) %in% variables)]
          trainboth$age_at_initial_pathologic_diagnosis = as.numeric(trainboth$age_at_initial_pathologic_diagnosis)

          #check if have any nas
          trainboth$OS.time = as.numeric(trainboth$OS.time)
          trainboth$OS = as.numeric(trainboth$OS)

          #remove = c()

          #perc5 = 0.05*dim(trainboth)[1]
          #make sure at least 5% of patients appear in each level of each variable

          #      for(b in 1:ncol(trainboth)){
          #        if(!(colnames(trainboth)[b] %in% c("age_at_initial_pathologic_diagnosis", "OS", "OS.time"))){
          #          dg = as.data.table(table(trainboth[,b]))
                    #print(dg)
          #          z = which(dg$N <perc5)
          #          if(!(length(z)==0)){
          #            var = dg$V1[z]
          #            namevar = colnames(trainboth)[b]
          #            remove = c(remove, namevar)
          #          }
          #        }
          #      }

          #      if(!(length(remove)==0)){
          #        trainboth = trainboth[,-(which(colnames(trainboth) %in% remove))]
          #      }
          #      print("pass10")


          both = try(coxph(Surv(OS.time, OS)  ~ ., data = trainboth), TRUE)
                if(!(inherits(both,"try-error"))){
                  testlncs$patient = rownames(testlncs)
                  testclin$patient = rownames(testclin)
                  testboth = merge(testlncs, testclin, by=c("patient", "OS", "OS.time"))
                  testboth$OS = as.numeric(testboth$OS)
                  testboth$OS.time = as.numeric(testboth$OS.time)
                  testboth$age_at_initial_pathologic_diagnosis = as.numeric(testboth$age_at_initial_pathologic_diagnosis)
                  pred_validation = predict(both, newdata = testboth)
                  # Determine concordance
                  cindex_validation = concordance.index(pred_validation, surv.time = testboth$OS.time,
                                                        surv.event=testboth$OS, method = "noether", na.rm=TRUE)
                  cindex_validation = cindex_validation$c.index
                  cinds_combined = cindex_validation

                  all_cinds = c(cinds_combined[1], cinds_justlncs[1], cinds_clin[1])
                  genes = genes_keep

                  #make a list, first component is genes kept, second component are the 3 cindices
                  res = list()
                  res[[1]] = genes
                  res[[2]] = all_cinds
                  print(all_cinds)
                  names(res) = c(paste(dat$Cancer[1], "genes", sep="_"), paste(dat$Cancer[1], "cinds", sep="_"))
                  #print(res)
                  return(res)

                } #both clin and lncs survival model no error
              } #just clin no error
            }
          }#dim(pvalues)[1] >2
        }
      }
    }
}#end main elastic net

print("main elastic net script loaded YAY")

#[4]  ----- prep full permutations of elastic net  ------------------
#main function that combines everything
#--------------------------------------------------------------------

source("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/survival_lncs_permutation_analysis_0214.R")

random_permutations = function(canc){ #main permutation cross-validation function, which basically runs elastic net function on 100 random datasets

  print(canc)

  #generate random dataset
  canc_datas = get_canc(canc)

  #check covariates
  dato = prepare_dat(canc_datas)
  print(dim(dato))

  print("start permutations")

  set.seed(101)
  run_res = replicate(1000, main_elastic_net(dato)) #DOUBLE CHECK number of replciations

  print("done permutations")

  #collect the lncRNAs that appear in more than 50% of the time
  lncrnas = unlist(run_res)
  z = which(str_detect(lncrnas, "ENSG"))
  lncrnas = lncrnas[z]
  #print(lncrnas)

  genes_list = as.data.table(table(unlist(lncrnas)))
  genes_list = genes_list[order(N)]
  print(as.data.frame(genes_list))
  print(paste(genes_list$V1, genes_list$N))
  genes_list_full = genes_list

  #genes_list = as.data.table(filter(genes_list, N >= 50)) #CHANGE after

  vals = as.numeric(unlist(run_res))
  cindices = vals[which(!(is.na(vals)))]
  cindices = matrix(cindices, ncol=3, byrow=T)
  colnames(cindices) = c("combined", "lncRNAs", "clinical")

  if(!(dim(genes_list)[1]==0)){

  genes_list$canc = dato$Cancer[1]
  genes_list$combo = paste(genes_list$V1, genes_list$canc, sep="_")

  combos = unique(genes_list$combo)
  print(combos)

  #survival function
  get_km_plot = function(comb){

  print(comb)
  gene = unlist(strsplit(comb, "_"))[1]
  canc = unlist(strsplit(comb, "_"))[2]
  cancer = rna$type[rna$Cancer == canc][1]

  newdat = dato

  z = which(colnames(newdat) %in% gene)

  if(!(length(z)==0)){

    cols_keep = c(gene, "patient", "age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "OS", "OS.time")
    newdat = newdat[,..cols_keep]
    newdat$age_at_initial_pathologic_diagnosis = as.numeric(newdat$age_at_initial_pathologic_diagnosis)

              #remove columns with less than 2 contrasts
              check_contrasts = function(col){
                check = dim(table(col))
                if(check >1){
                  return("keep")
                }
              }

    keep = unlist(apply(newdat, 2, check_contrasts))
    z = which(colnames(newdat) %in% names(keep))
    newdat = newdat[,..z]

    newdat$OS = as.numeric(as.character(newdat$OS))
    newdat$OS.time = as.numeric(as.character(newdat$OS.time))

    z = which(colnames(newdat) %in% gene)
    colnames(newdat)[z] = "gene"
    newdat$gene = as.numeric(newdat$gene)

    rownames(newdat) = newdat$patient
    newdat$patient = NULL
    newdat=as.data.frame(newdat)

    #split patients
    med = median(newdat$gene)

    #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1
    l1 = which(newdat[,z] > 0)
    l2 = which(newdat[,z] ==0)
    newdat[l1,z] = 1
    newdat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(newdat[,z] >= med)
    l2 = which(newdat[,z] < med)
    newdat[l1,z] = 1
    newdat[l2, z] = 0
    }

    cox_mod = coxph(Surv(OS.time, OS) ~ ., data = newdat)

    #fit <- survfit(Surv(newdat$OS.time, newdat$OS) ~ gene,
    #           data = newdat)
    # Visualize with survminer
    #ggsurvplot(fit, data = newdat, risk.table = TRUE, pval = TRUE)

    print(glance(cox_mod)$concordance)

    lnc_only = glance(coxph(Surv(OS.time, OS) ~ gene, data = newdat))$concordance

    conc = round(glance(cox_mod)$concordance, digits=2)
    wald_p = summary(cox_mod)$coefficients[1,5]
    hr = round(summary(cox_mod)$coefficients[1,2], digits=4)

    results = c(cancer, gene, canc, lnc_only, wald_p, hr)
    names(results) = c("type", "gene", "cancer", "concordance", "wald_p", "hr")

    return(results)
    print(paste("done", comb))
    }
    }

  surv_results = llply(combos, get_km_plot, .progress="text")

  perms_surv_results1 = as.data.table(ldply(surv_results))
  perms_surv_results1$wald_p = as.numeric(perms_surv_results1$wald_p)
  #perms_surv_results1 = perms_surv_results1[order(wald_p)]
  perms_surv_results1$fdr = p.adjust(perms_surv_results1$wald_p, method = "holm")
  print(paste("done", canc, "round of 1000 cross-validations"))
  perms_surv_results1$round = paste(sample(1:5000000, 1), sample(1:5000000, 1), sample(1:5000000, 1), sep="_")
  perms_surv_results1$inference_pvals = 1

  colnames(genes_list)[1:2] = c("gene", "num_rounds_selected")
  perms_surv_results1 = merge(perms_surv_results1, genes_list, by="gene")

  perms_surv_results1 = perms_surv_results1[order(num_rounds_selected)]
  cancer =canc_conv$type[canc_conv$Cancer==canc]
  output="/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2021_manuscript/real_elastic_net_runs_1000/"
  file = paste(output, "round", cancer, date, perms_surv_results1$round[1], "genes", ".rds", sep="_")
  saveRDS(perms_surv_results1, file)

  #save c-indices results
  cindices = as.data.frame(cindices)
  cindices$canc = canc
  file = paste(output, "round", cancer, date, "cindices", ".rds", sep="_")
  saveRDS(cindices, file)

  print(head(perms_surv_results1))
  print(head(cindices))

  #END
}
}#<---- make function


#RUN 1000 ELASTIC-NET CROSS VALIDATIONS
#FOR EACH CANCER TYPE

random_permutations(canc)
print(paste("done", canc))
