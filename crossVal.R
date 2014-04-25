# source('~/Documents/Education/ULB/Phd/Code/Start.R')
# path<-getPath()
source(paste(path,"Code/metrics/detection.R",sep=""))
source(paste(path,"Code/caret/preprocessCaret.R",sep=""))
source(paste(path,"Code/caret/confMatrix.R",sep=""))
source(paste(path,"Code/caret/trainCaret.R",sep=""))
source(paste(path,"Code/utility/Andrea.R",sep=""))
library(unbalanced)
library(caret)


#examples
# config <- list(nFolds = 10, 
#                 balance = list(type = "ubUnder", percOver = 200, percUnder = 200, k = 5, perc = 50, method = "percPos", w = NULL), 
#                 nCores = 5,
#                verbose = T)
# library(mlbench)
# data(Sonar)
# library(parallel)
# library(doParallel)
# cl <- makeCluster(config$nCores)
# registerDoParallel(cl)
# crossVal(algo="rf", data = Sonar, tgt = ncol(Sonar), positive = "R", config)
# stopCluster(cl)


#cross validation with positive class proportion equal in all the folds
#algo: ml algo used in caret pkg
#data: dataset
#tgt: index of the target/response variable in the dataset
#positive: positive (minority) class
#balanced: logical, if TRUE rebalance within each fold of the cross validation
crossVal <- function(algo = "rf", data, tgt = ncol(data), positive = 1, config, balanced = FALSE, ...){
  
  k <- config$nFolds
  y <- data[ ,tgt]
  
  fold.id <- createFolds(y, k, list = FALSE)
  
  #       results<- NULL
  #       for(i in 1:k){    
  #         res <- trainPredFold(algo, data, tgt, positive, fold.id, i, config, balanced, ...)
  #         results <- cbind(results, res)
  #       }
  
  #test all folds in parallel, need a registered cluster
  require(foreach)
  algoPkg <- getModelInfo(algo, regex = FALSE)[[1]]$library
  results <- foreach(i=1:k, .packages=c(algoPkg, "ROCR", "caret", "PerfMeas", "unbalanced"), .combine = 'cbind',
                     .export=c('trainPredFold','learn','singleGrid','getGrid','test',
                               'getClasMetrics','getProbMetrics','confusionMatrix','confusionMatrix.table',
                               'tabValues','mcc','print.confusionMatrix','classMetrics','confMatrix',
                               'getConfMatrixMetrics','getConfMatrixMetric','MCC','Chisquare','probMetrics','APk',
                               'normCumPR','AUCPR','brierScore','stratBrierScore','cumPR','minAuprc','probBias','AUCr','id.row.has.na'), 
                     .inorder = FALSE) %dopar% { trainPredFold(algo, data, tgt, positive, fold.id, i, config, balanced, ...) }
  
  return(results)
}



trainPredFold <- function(algo, data, tgt, positive, fold.id, i, config, balanced, ...){
  
  trainIndex <- which(fold.id != i)
  testIndex <- which(fold.id == i)
  training <- data[trainIndex, ]
  testing <- data[testIndex, ]
  
  #remove NAs
  training <- training[setdiff(1:nrow(training), id.row.has.na(training)), ]
  testing <- testing[setdiff(1:nrow(testing), id.row.has.na(testing)), ]
  
  if(balanced){
    #re-balance the dataset
    data <- ubBalance(training[ ,-tgt], training[ ,tgt], type = config$balance$type, positive, 
                      config$balance$percOver, config$balance$percUnder, config$balance$k, 
                      config$balance$perc, config$balance$method, config$balance$w,
                      config$verbose)
    training <- data.frame(data$X, Y = data$Y)
    testing <- data.frame(testing[ ,-tgt], Y = testing[ ,tgt])
    tgt <- ncol(training)
  }
  
  #make sure they have the same colnames if not balanced
  colnames(training)[tgt] <- "Y"
  colnames(testing) <- colnames(training)
  
  #transform output as factor
  training[ ,tgt] <- factor(training[ ,tgt] == positive, labels=c("NEG","POS"))
  testing[ ,tgt] <- factor(testing[ ,tgt] == positive, labels=c("NEG","POS"))
  
  #pre-process
  # processed <- preProc(Y ~ ., training, method = c("center", "scale"), dummy = TRUE)
  # training <- processed$data    
  # testing <- preProc(Y ~ ., testing,  processed$preProcValues, processed$preProcVars)$data
  
  #train a model 
  fit <- learn(Y ~ ., training, algo, verbose = config$verbose, ...)
  
  #predict
  resList <- test(Y ~ ., testing , fit, positive = "POS")
  results <- c(resList$clasMetrics$metrics, resList$probMetrics, trainTime = resList$trainTime, predTime = resList$predTime)
  
  return(results)
}