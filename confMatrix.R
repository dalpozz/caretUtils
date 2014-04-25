#code originally from caret package - March 2014

#need caret loaded to use functions: 
# sensitivity.table(), 
# specificity.table(), 
# posPredValue.table(), 
# negPredValue.table()
# as.matrix.confusionMatrix()
# as.table.confusionMatrix()

#confusionMatrix is from confusionMatrix.default	
#code modified: confusionMatrix.table, print.confusionMatrix

confusionMatrix <- function(data, reference,
                            positive = NULL,
                            dnn = c("Prediction", "Reference"),
                            prevalence = NULL,
                            ...)
{
  library(e1071)
  if(!is.factor(data)) data <- factor(data)
  if(!is.factor(reference)) reference <- factor(reference)
  if(!is.character(positive) & !is.null(positive)) stop("positive argument must be character")
  
  if(length(levels(data)) != length(levels(reference)))
    stop("the data and reference factors must have the same number of levels")
  
  if(any(levels(data) != levels(reference)))
    stop("the data and reference values must have exactly the same levels")
  
  classLevels <- levels(data)
  numLevels <- length(classLevels)
  if(numLevels < 2) 
    stop("there must be at least 2 factors levels in the data")
  
  if(numLevels == 2 & is.null(positive))  positive <- levels(reference)[1]
  
  classTable <- table(data, reference, dnn = dnn, ...)
  
  confusionMatrix.table(classTable, positive, prevalence = prevalence)
}


confusionMatrix.table <- function(data, positive = NULL, prevalence = NULL, ...){
  
  library(caret)
  library(e1071)
  
  if(length(dim(data)) != 2) stop("the table must have two dimensions")
  if(!all.equal(nrow(data), ncol(data))) stop("the table must nrow = ncol")
  if(!all.equal(rownames(data), colnames(data))) stop("the table must the same classes in the same order")
  if(!is.character(positive) & !is.null(positive)) stop("positive argument must be character")
  
  classLevels <- rownames(data)
  numLevels <- length(classLevels)
  if(numLevels < 2) 
    stop("there must be at least 2 factors levels in the data")
  
  if(numLevels == 2 & is.null(positive))  positive <- rownames(data)[1]
  
  
  if(numLevels == 2 & !is.null(prevalence) && length(prevalence) != 1)
    stop("with two levels, one prevalence probability must be specified")
  
  if(numLevels > 2 & !is.null(prevalence) && length(prevalence) != numLevels)
    stop("the number of prevalence probability must be the same as the number of levels")
  
  if(numLevels > 2 & !is.null(prevalence) && is.null(names(prevalence)))
    stop("with >2 classes, the prevalence vector must have names")
  
  propCI <- function(tab)
  {
    n <- sum(tab)
    x <- sum(diag(tab))
    if(x > n){
      print(tab)
      stop("n:", n, "must be bigger than x:", x)
    }
    binom.test(x, n)$conf.int
  }
  
  propTest <- function(tab)
  {
    n <- sum(tab)
    x <- sum(diag(tab))
    if(x > n){
      print(tab)
      stop("n:", n, "must be bigger than x:", x)
    }
    out <- binom.test(x, n, p = max(apply(tab, 2, sum)/n), alternative = "greater")
    unlist(out[c("null.value", "p.value")])
  }
  
  overall <- c(
    unlist(classAgreement(data))[c("diag", "kappa")],
    #propCI(data),
    #propTest(data),
    mcnemar.test(data)$p.value)
  
  #names(overall) <- c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue")  
  names(overall) <- c("Accuracy", "Kappa", "McnemarPValue")  
  
  if(numLevels == 2){
    if(is.null(prevalence)) prevalence <- sum(data[, positive])/sum(data)
    negative <- classLevels[!(classLevels %in% positive)]
    # tableStats <- c(sensitivity.table(data, positive),
    # specificity.table(data, negative),
    # posPredValue.table(data, positive, prevalence = prevalence),
    # negPredValue.table(data, negative, prevalence = prevalence),
    # prevalence,
    # sum(data[positive, positive])/sum(data),
    # sum(data[positive, ])/sum(data))
    # names(tableStats) <- c("Sensitivity", "Specificity",
    # "Pos Pred Value", "Neg Pred Value",
    # "Prevalence", "Detection Rate",
    # "Detection Prevalence")   
    # tableStats["Balanced Accuracy"] <- (tableStats["Sensitivity"]+tableStats["Specificity"])/2
    
    Recall <- sensitivity.table(data, positive)
    Precision <-  posPredValue.table(data, positive, prevalence = prevalence)
    Fmeasure <- ifelse((Precision + Recall) != 0, 2*(Precision * Recall)/(Precision + Recall), NA)
    TrueNegativeRate <- specificity.table(data, negative)
    Mcc <- mcc(data, positive)
    Gmean <- sqrt(Recall * TrueNegativeRate)
    #Precision normalized by its random Precision which is equal to the prevalence
    normPrecision <- ifelse(prevalence != 0, Precision/prevalence, 0)
    BalancedAccuracy <- (Recall + TrueNegativeRate)/2
    BER <- ( 1- TrueNegativeRate + 1- Recall)/2  
    tableStats <- c(Fmeasure, 
                    Precision, 
                    Recall, 
                    TrueNegativeRate, 
                    Mcc, 
                    Gmean, 
                    normPrecision,
                    negPredValue.table(data, negative, prevalence = prevalence),
                    prevalence,
                    sum(data[positive, positive])/sum(data),
                    sum(data[positive, ])/sum(data),
                    BER)
    names(tableStats) <- c("Fmeasure", #"F-measure / F1-score"
                           "Precision", # Precision / PosPredValue"
                           "Recall", #"Recall / TruePositiveRate / Sensitivity"
                           "TrueNegativeRate", #"TrueNegativeRate / Specificity"
                           "Mcc", 
                           "Gmean", 
                           "normPrecision", 
                           "NegPredValue", 
                           "Prevalence", #"Prevalence / PosProportion"
                           "Detection Rate", 
                           "Detection Prevalence", 
                           "BER")   
    
    values <- tabValues(data, positive)
    
  } else {
    stop("code for numLevels > 2 implemented in Caret, but removed in this function for the moment")
  }
  
  structure(list(
    positive = positive,
    table = data, 
    tabValues = values,
    overall = overall, 
    byClass = tableStats,
    dots = list(...)), 
    class = "confusionMatrix")
}



mcc <- function(tab, pos = colnames(tab)[1])
{
  if(nrow(tab) != 2 | ncol(tab) != 2) stop("A 2x2 table is needed")
  neg <- colnames(tab)[colnames(tab) != pos]
  tp <- tab[pos, pos]
  tn <- tab[neg, neg]
  fp <- tab[pos,neg]
  fn <- tab[neg, pos]
  d1 <- tp + fp
  d2 <- tp + fn
  d3 <- tn + fp
  d4 <- tn + fn
  if(d1 == 0 | d2 == 0 | d3 == 0 | d4 == 0) return(0)
  ((tp * tn) - (fp * fn))/sqrt(d1*d2*d3*d4)  
}


tabValues <- function(tab, pos = colnames(tab)[1])
{
  if(nrow(tab) != 2 | ncol(tab) != 2) 
    stop("A 2x2 table is needed")
  neg <- colnames(tab)[colnames(tab) != pos]
  tp <- tab[pos, pos]
  tn <- tab[neg, neg]
  fp <- tab[pos, neg]
  fn <- tab[neg, pos]
  return(c(TP = tp, FP = fp, TN = tn, FN = fn))
}


print.confusionMatrix <- function(x, digits = 3, printStats = TRUE, ...){
  cat("Confusion Matrix and Statistics\n\n") 
  print(x$table, ...)
  
  #   cat("\n(columns are reference results, rows are predictions)\n")
  
  if(printStats){
    
    tmp <- round(x$overall, digits = digits)
    pIndex <- grep("PValue", names(x$overall))
    tmp[pIndex] <- format.pval(x$overall[pIndex], digits = digits)
    overall <- tmp    
    
    overallText <- c(
      paste(overall["Accuracy"]),
      paste(overall["Kappa"]),
      paste(overall["McnemarPValue"]))
    
    overallNames <- c("Accuracy", 
                      "Kappa",
                      "Mcnemar's Test P-Value")
    
    if(dim(x$table)[1] > 2)
    {
      cat("\nOverall Statistics\n")
      overallNames <- ifelse(
        overallNames == "",
        "",
        paste(overallNames, ":"))
      out <- cbind(
        format(overallNames, justify = "right"),
        overallText)
      colnames(out) <- rep("", ncol(out))
      rownames(out) <- rep("", nrow(out))
      
      print(out, quote = FALSE)
      
      cat("\nStatistics by Class:\n\n")
      print(t(x$byClass), digits = digits)
      
    } else {
      
      overallText <- c(format(x$byClass, digits = digits), overallText)
      overallNames <- c(names(x$byClass), overallNames)
      overallNames <- ifelse(overallNames == "", "", paste(overallNames, ":"))
      
      overallNames <- c("'Positive' Class :", overallNames)
      overallText <- c(x$positive, overallText)
      
      out <- cbind(format(overallNames, justify = "right"), overallText)
      colnames(out) <- rep("", ncol(out))
      rownames(out) <- rep("", nrow(out))
      
      out <- rbind(out, rep("", 2))
      
      print(out, quote = FALSE)
      
    }
    
    
  }
  invisible(x)   
}
