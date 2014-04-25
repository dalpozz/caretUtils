#*********************************************************************
# This is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under GPL conditions.
# Copyright (C) Andrea Dal Pozzolo, Feb 2014
# License: GPL (>= 2)                                           
#*********************************************************************
#train a model with caret
library(caret)

# # example
# library(mlbench)
# data(Sonar)
# #create a stratified random sample of the data into training and test sets
# inTraining <- createDataPartition(Sonar$Class, p = 0.75, list = FALSE)
# training <- Sonar[inTraining, ]
# testing <- Sonar[-inTraining, ]
# 
# #tuning parameter grid
# # gbmGrid <- createGrid("gbm") #does not exist anymore
# gbmGrid <- data.frame(interaction.depth = 4, n.trees = 100, shrinkage = .1)
# # 10-fold CV repeated 5 times
# gbmControl <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary)
# 
# fit <- learn(Class ~ ., training, "gbm")
# #fit2 <- learn(Class ~ ., training, "gbm", trControl = gbmControl, tuneGrid = gbmGrid)
# res <- test(Class ~ ., testing, fit, positive = "M")
# #res2 <- test(Class ~ ., testing, list(fit, fit2), positive = "M")



#return the tuning grid for a given algorithm and dataset
getGrid <- function(formula, data, algo, len = 3){ 
  tgt <- which(names(data) == as.character(formula[[2]]))
  check <- modelLookup(algo)
  modules <- getModelInfo(algo, regex = FALSE)[[1]]
  grid <- modules$grid(data[ ,-tgt], data[ ,tgt], len)  
  return(grid)  
}

#return a grid with a single value for each tuning parameters
#take the value that is the median of all possible values
singleGrid <- function(formula, data, algo, len = 3){
  grid <- getGrid(formula, data, algo, len)
  oneParam <- subset(grid, (1:nrow(grid)) == median(1:nrow(grid)))
  return(oneParam)
}

#**********************************
#train a model using caret function
# trControl: see caret::trainControl()
# ... is used to pass non-tuning parameters for the selected algorithm (i.e. ntree in random forest)
#**********************************
learn <- function(formula, training, algo = "gbm", 
                  tuneGrid = NULL,
                  #trControl = trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary), 
                  trControl = NULL,  
                  verbose = FALSE,
                  ...){
  
  
  tgt <- which(names(training) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable not defined")
  task <- ifelse(is.factor(training[ ,tgt]), "class", "reg")
  
  if(is.null(trControl)){
    #single prediction
    trControl <- trainControl(method = "none", 
                              classProbs = TRUE, 
                              summaryFunction = twoClassSummary,             
                              returnResamp = "none")
    if(is.null(tuneGrid))
      tuneGrid <- singleGrid(formula, training, algo)
    #stop("tuneGrid has to be provided when only one model is used for prediction")
    if(nrow(tuneGrid) > 1)
      stop("tuneGrid has to have only one value for each parameters when parameter tunining is not performed")
  }
  
  if(is.null(tuneGrid))
    tuneGrid <- getGrid(formula, training, algo, len = 3)
  
  #if trControl$classProbsm, use AUROC to find the best parameters
  metric <- ifelse(trControl$classProbs, "ROC", ifelse(task == "class", "Accuracy", "RMSE"))
  maximize <- ifelse(metric == "ROC", TRUE, ifelse(metric == "RMSE", FALSE, TRUE))
  
  #remove row of the training with at least one NA
  training <- training[setdiff(1:nrow(training), which(is.na(training)) %% nrow(training)), ]
  
  #Using the formula interface with a large number of predictors may slow the computations.
  #fit <- train(formula, data = training,
  fit <- train(x = training[ ,-tgt] , y = training[ ,tgt],
               method = algo,
               trControl = trControl,
               tuneGrid = tuneGrid,
               metric = metric, 
               maximize = maximize,
               #verbose = verbose, argument not matched with rpart
               ...)
  
  #model <- fit$finalModel
  
  return(fit)
}



#**********************************
#test a trained model on a testing set
# type: either "raw" or "prob" (Class probabilities are not available for all classification models)
# fit can be a single model or a list of models from caret, in the second case if a cluster is registered then it is done in parallel
#**********************************
test <- function(formula, testing, fit, type = "prob", positive = NULL, cutoff = 0.5, verbose = FALSE, ...){
  
  if(class(fit) != "train" & class(fit) != "list")
    stop("model fit of class ",class(fit)," not supported")
  
  tgt <- which(names(testing) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable not included")
  
  #remove row of the testing with at least one NA
  testing <- testing[setdiff(1:nrow(testing), which(is.na(testing)) %% nrow(testing)), ]
  
  predTime <- system.time({
  pred.prob <- NULL
  pred.class <- predict(fit, newdata = testing[ ,-tgt], type = "raw", verbose = verbose, ...)
  if(type == "prob")
    pred.prob <- predict(fit, newdata = testing[ ,-tgt], type = "prob", verbose = verbose, ...)
  })
  predTime <- as.numeric(predTime[3])
  
  Nmodel <- ifelse(is.list(pred.class), length(pred.class), 1)
  
  if(Nmodel == 1){
    if(is.list(pred.class)){
      pred.class <- unlist(pred.class)
      pred.prob <- data.frame(pred.prob[[1]])
    }
    clasMetrics <- getClasMetrics(pred.class, testing[ ,tgt], positive, verbose)
    probMetrics <- getProbMetrics(pred.prob, testing[ ,tgt], positive, verbose)
    trainTime <- fit$times$everything[[3]]
  }
  else{
    clasMetrics <- lapply(pred.class, getClasMetrics, testing[ ,tgt], positive, verbose)
    probMetrics <- lapply(pred.prob, getProbMetrics, testing[ ,tgt], positive, verbose) 
    trainTime <- sapply(fit, function(fit) fit$times$everything[[3]])
  }
  
  return(list(clasMetrics = clasMetrics, probMetrics = probMetrics, pred.prob = pred.prob, pred.class = pred.class, real.class = testing[ ,tgt], trainTime = trainTime, predTime = predTime))
}



getClasMetrics <- function(pred.class, real.class, positive, verbose = FALSE){

  stopifnot(is.factor(real.class), is.factor(pred.class))
  lev <- levels(real.class)
  if(all(lev %in% c(1, 0)))
    positive <- 1
  if(all(lev %in% c(1, -1)))
    positive <- 1
  if(all(lev %in% c(TRUE, FALSE)))
    positive <- TRUE
    
  if(is.null(positive))
    stop("need to specify the positive class")
    
  confMat <- suppressWarnings(confusionMatrix(data = pred.class, real.class, positive = positive))
  metrics <- c(confMat$tabValues, confMat$overall, confMat$byClass) #too many metrics
  #TO DO: select a subset of metrics from the previous results
  
#   pred.bin <- as.numeric(pred.class == positive)
#   real.bin <- as.numeric(real.class == positive)
#   metrics <- unlist(classMetrics(pred.bin, real.bin, verbose))
  
  if(verbose){
    print(confMat)
    print(round(metrics, digits = 3))
  }
  
  return(list(confMatrix = confMat, metrics = metrics))
  
}


getProbMetrics <- function(pred.prob, real.class, positive, verbose){
  
  stopifnot(is.factor(real.class))
  
  lev <- levels(real.class)
  if(all(lev %in% c(1, 0)))
    positive <- 1
  if(all(lev %in% c(1, -1)))
    positive <- 1
  if(all(lev %in% c(TRUE, FALSE)))
    positive <- TRUE
  
  if(is.null(pred.prob))
    retun(NULL)
  
  phat.1 <- pred.prob[ ,which(names(pred.prob) == positive)]
  real.bin <- as.numeric(real.class == positive)
  metrics <- probMetrics(phat.1, real.bin)
  
  if(verbose){
    probClassData <- data.frame(Y = real.class, phat.1)
    liftData <- lift(Y ~ phat.1, data = probClassData)
    pdf("lift.pdf")
    plot(liftData, values = 60, auto.key = list(columns = 3, lines = TRUE, points = FALSE))
    dev.off()
    calData <- calibration(Y ~ phat.1, data = probClassData, cuts = 13)
    pdf("calibration.pdf")
    plot(calData, type = "l", auto.key = list(columns = 3, lines = TRUE, points = FALSE))
    dev.off()
  }
  
  return(unlist(metrics))
}


