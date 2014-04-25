#*********************************************************************
# This is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under GPL conditions.
# Copyright (C) Andrea Dal Pozzolo, Feb 2014
# License: GPL (>= 2)                                           
#*********************************************************************
#functions to pre-process the data
library(caret)



#**********************************
# pre-process the data
# if preProcValues is not null, then apply preProcValues prepocessing steps to the data
#**********************************
## example 1:
# data(mdrr)
# data <- cbind(mdrrDescr,mdrrClass)
# newData <- preProc(mdrrClass~. , data)$data
## example 2: apply the same transformation to training and test set
# inTrain <- sample(seq(along = mdrrClass), length(mdrrClass)/2)
# training <- cbind(mdrrDescr[inTrain, ],  mdrrClass = mdrrClass[inTrain])
# test <- cbind(mdrrDescr[-inTrain, ], mdrrClass = mdrrClass[-inTrain])
# processed <- preProc(mdrrClass~. , training, method = c("center", "scale"))
# trainTransformed <- processed$data
# preProcValues <- processed$preProcValues
# preProcVars <- processed$preProcVars
# testTransformed <- preProc(mdrrClass~. , test,  preProcValues, preProcVars)$data
#**********************************
preProc <- function(formula, data, preProcValues = NULL, preProcVars = NULL, method = c("center", "scale"),
                    dummy = FALSE, nearZeroRemove = TRUE, highCorrRemove = FALSE, verbose = TRUE, ...){
  
  if(class(formula)!="formula")
    stop("formula not recognized")
  
  if(dummy)
    data <- createDummy(formula, data)
  
  if(is.null(preProcVars)){  
    if(nearZeroRemove)
      data <- nearZeroInputRemove(formula, data, verbose)
    
    if(highCorrRemove)
      data <- highCorrInputRemove(formula, data, cutoff = 0.95, verbose)  
  } else {
    subVars <- which(colnames(data) %in% preProcVars)
    data <- data[ ,subVars]
  }
  
  tgt <- which(colnames(data) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable ",as.character(formula[[2]])," not found")
  tgtName <- colnames(data)[tgt]
  input <- data[ ,-tgt]
  
  if(is.null(preProcValues)){
    if(verbose)
      cat("preProcess method: ", method, "\n")
    
    preProcValues <- caret::preProcess(input, method = method, verbose = verbose, ...)
  }
  
  #apply preProcValues to the data.
  newinput <- predict.preProcess(preProcValues, input)
  
  newdata <- data.frame(newinput, data[ ,tgt])
  colnames(newdata) <- c(colnames(newinput), tgtName)
  
  return(list(data = newdata, preProcValues = preProcValues, preProcVars = colnames(newdata)))
}


#**********************************
#create dummy variables for factors 
#one variable for each factor level
#**********************************
# example:
# library(earth)
# data(etitanic)
# createDummy(survived ~ ., etitanic)
#**********************************
createDummy <- function(formula, data){
  
  if(class(formula) != "formula")
    stop("formula not recognized")
  
  #the column where the target variable is
  tgt <- which(colnames(data) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable ",as.character(formula[[2]])," not found")
  tgtName<-colnames(data)[tgt]
  
  dummies <- caret::dummyVars(formula, data)
  newdata <- predict(dummies, newdata = data)
  
  numData <- data.frame(newdata, data[ ,tgt])
  colnames(numData) <- c(colnames(newdata),tgtName)
  
  return(numData)
}


#**********************************
# the function removes predictors with:
#   1) Near Zero-Variance
# DISCARGED: 2) empty conditional distributions in at least one class of y.
# Zero- and Near Zero-Variance Predictors:
#   -predictors that have one unique value (i.e. are zero variance predictors) or 
#   -predictors that are have both:
#       1)have very few unique values relative to the number of samples 
#       2)the ratio of the frequency of the most common value to the frequency of the second most common value is large.
#predictors may become zero-variance predictors when the data are split into cross-validation/bootstrap sub-samples
#or that a few samples may have an undue influence on the model. 
#These "near-zero-variance" predictors may need to be identified and eliminated prior to modeling
#**********************************
# example:
# data <- data.frame(x1 = rep(c(0, 1), 45),
#                 x2 = c(rep(0, 10), rep(1, 80)), 
#                 x3 = c(rep(0, 60), rep(1, 30)),
#                 x4 = rep(0, 90),
#                 y = factor(rep(letters[1:3], each = 30)))
# newdata <- nearZeroInputRemove(y ~. ,data)
#**********************************
nearZeroInputRemove <- function(formula,data,verbose=TRUE){
  
  if(class(formula)!="formula")
    stop("formula not recognized")
  
  tgt <- which(colnames(data) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable ",as.character(formula[[2]])," not found")
  tgtName <- colnames(data)[tgt]
  classes <- data[ ,tgt]
  input <- data[ ,-tgt]
  n<-ncol(input)
  
  nzv <- nearZeroVar(input)
  filteredInput <- input[ ,setdiff(1:n, nzv)]
  if(verbose & length(nzv) > 0){
    cat("predictors with near zero variance: ", length(nzv), "of", ncol(input), "\n")
    if(length(nzv) < 50) 
      cat(colnames(input)[nzv], "\n")
    #for(i in 1:length(nzv))
    #   print(table(input[ ,nzv[i]], y = classes))
  }
  
  #   #identifies factors that are sparse within groups of y
  #   posCond <- checkConditionalX(filteredInput, classes)
  #   if(ncol(filteredInput) > 0){
  #     nullCond <- setdiff(1:ncol(filteredInput), posCond)
  #     if(verbose & length(nullCond) > 0){
  #       cat("predictors with null conditional distribution removed: ", length(nullCond), "of", ncol(filteredInput), "\n")
  #       if(length(nullCond) < 50)
  #         cat(colnames(filteredInput)[nullCond], "\n")
  #       #for(i in 1:length(nullCond))
  #       #   print(table(filteredInput[ ,nullCond[i]], y = classes))
  #     }
  #   }
  #   
  #   filteredInput <- filteredInput[ ,posCond])
  
  newData <- data.frame(filteredInput, classes)
  colnames(newData) <- c(colnames(filteredInput), tgtName)
  
  return(newData)
}



#**********************************
#Identifying Correlated Predictors and remove them
#**********************************
# example:
# data(mdrr)
# data<-cbind(mdrrDescr,mdrrClass)
# newData<-highCorrInputRemove(mdrrClass~. , data)
#**********************************
highCorrInputRemove <- function(formula, data, cutoff = 0.95, verbose=TRUE){
  
  if(class(formula) != "formula")
    stop("formula not recognized")
  
  tgt <- which(names(data) == as.character(formula[[2]]))
  if(length(tgt) == 0)
    stop("target variable ",as.character(formula[[2]])," not found")
  tgtName <- names(data)[tgt]
  input <- data[ ,-tgt]
  
  inputCor <- cor(input)
  highlyCorDescr <- findCorrelation(inputCor, cutoff)
  
  if(verbose & length(highlyCorDescr)>0){
    cat("predictors highly corrrelated removed: ", length(highlyCorDescr), "of", ncol(input), "\n")
    if(length(highlyCorDescr) < 50) 
      cat(colnames(input)[highlyCorDescr], "\n")
  }
  newinput <- input[ ,setdiff(1:ncol(input), highlyCorDescr)]
  
  newData <- data.frame(newinput, data[ ,tgt])
  colnames(newData) <- c(colnames(newinput), tgtName)
  
  return(newData) 
}