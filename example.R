# example
source('~/Documents/Education/ULB/Phd/Code/Start.R')
path<-getPath()
source(paste(path,"Code/detection.R",sep=""))
source(paste(path,"Code/caret/trainCaret.R",sep=""))
library(caret)
library(mlbench)
data(Sonar)
#create a stratified random sample of the data into training and test sets
inTraining <- createDataPartition(Sonar$Class, p = 0.75, list = FALSE)
training <- Sonar[inTraining, ]
testing <- Sonar[-inTraining, ]
tgt <- which(names(Sonar) == "Class")

#tuning parameter grid
# gbmGrid <- createGrid("gbm") #does not exist anymore
# gbmSingleGrid <- data.frame(interaction.depth = 4, n.trees = 100, shrinkage = .1)
# modules <- getModelInfo("gbm", regex = FALSE)[[1]]
# gbmGrid <- modules$grid(training[ ,-tgt], training[ ,tgt], len = 3)


#single prediction without parameter tuning
fit <- learn(Class ~ ., training, "gbm")
res <- test(Class ~ ., testing, fit)

# 10-fold CV with parameter tuning
gbmGrid <- getGrid(Class ~ ., training, "gbm", len = 3)
gbmControl <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary)
fit2 <- learn(Class ~ ., training, "gbm", trControl = gbmControl, tuneGrid = gbmGrid)
res <- test(Class ~ ., testing, fit2)



#speed comparison using RF, caret wrapper is slow

timecaretwrap<-system.time({
  fit <- learn(Class ~ ., training, "rf", tuneGrid = data.frame(mtry = 31), ntree = 200)	
})

timecaret<-system.time({
  fit2 <- train(Class ~ ., training, method = "rf", ntree = 200,
                trControl = trainControl(method = "none", returnData = FALSE, returnResamp = "none"),
                tuneGrid = data.frame(mtry = 31))
})

timerf<-system.time({
  rfmod <- randomForest(Class ~ ., training, mtry = 31, ntree = 200)
})

timecaretwrap2<-system.time({
  pred <- test(Class ~ ., testing, fit,  type = "raw")
})

timecaret2<-system.time({
  pred <-  predict.train(fit2, testing[ ,-tgt], type = "raw")
})

timerf2<-system.time({
  rfres <- predict(rfmod, testing, type ="response")
})

cat("Time to train with caret wrapper:",round(as.numeric(timecaretwrap[3]), digits=2),"\n")
cat("Time to train with caret:",round(as.numeric(timecaret[3]), digits=2),"\n")
cat("Time to train with rf:",round(as.numeric(timerf[3]), digits=2),"\n")
cat("Time to test with caret wrapper:",round(as.numeric(timecaretwrap2[3]), digits=2),"\n")
cat("Time to test with caret:",round(as.numeric(timecaret2[3]), digits=2),"\n")
cat("Time test rf:",round(as.numeric(timerf2[3]), digits=2),"\n")






