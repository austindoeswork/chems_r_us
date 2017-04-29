library("MASS") #Contains the Fisher LDA command 
library("readr")  #The library for read_csv
library("e1071") #SVM
# install.packages("e1071", repos="http://R-Forge.R-project.org")

# modify this for wherever your project is
path <- "/Users/austin/RPI/idm/proj/chems_r_us"
printf <- function(...) cat(sprintf(...))

# read the data
chemstrain <- read_csv(paste(path, "/data/chemsdata.csv", sep=""))
chemstest <- read_csv(paste(path, "/data/chemstest.csv", sep=""))

# QUESTION 2
# get the useful stuff
xtrain <- chemstrain[,-1]
xtest <- chemstest[,-1]
postrain <- xtrain[xtrain[,ncol(xtrain)] == 1,]
negtrain <- xtrain[xtrain[,ncol(xtrain)] == -1,]
trainclass <- chemstrain[, ncol(chemstrain)]

# QUESTION 2(b)
# boxplot(xtrain)

# QUESTION 2(c)
fishermethod <-function(D_train, D_test, y_train, y_test) {
  z <- lda(class ~ .,D_train,prior=c(1,1)/2)
  
  trainans <- predict(z,D_train)$class
  cctrain <- y_train - as.matrix(as.numeric(as.matrix(trainans)))
  trainfalseneg <- sum(cctrain >0) #Number of positive samples which the model classified as negative
  trainfalsepos <- sum(cctrain <0) #Number of negative samples which the model classified as positive
  trainerr <- (trainfalseneg+trainfalsepos)/nrow(D_train)
  printf("FISHER TRAIN ERR: %f\n", trainerr)
  
  testans <- predict(z,D_test)$class
  cctest <- y_test - as.matrix(as.numeric(as.matrix(testans)))
  testfalseneg <- sum(cctest >0) #Number of positive samples which the model classified as negative
  testfalsepos <- sum(cctest <0) #Number of negative samples which the model classified as positive
  testerr <- (testfalseneg+testfalsepos)/nrow(D_test)
  printf("FISHER TEST ERR: %f\n", testerr)
  
  thresh <- ((z$means[1,] + z$means[2,])/2)%*%z$scaling
  printf("FISHER THRESHOLD: %f\n", thresh)
  
  printf("FISHER NORM:\n")
  print(z$scaling)
}

# QUESTION 2(d)
supportvector <-function(D_train, D_test, y_train, y_test) {
  model <-svm(D_train, y_train)
  
  trainans <- predict(model,D_train)
  trainpred <- (trainans > 0) * 2 - 1
  cctrain <- y_train - as.matrix(trainpred)
  
  trainfalseneg <- sum(cctrain <0) #Number of positive samples which the model classified as negative
  trainfalsepos <- sum(cctrain >0) #Number of negative samples which the model classified as positive
  trainerr <- (trainfalseneg+trainfalsepos)/nrow(D_train)
  printf("SVM TRAIN ERR: %f\n", trainerr)
  
  testans <- predict(model,D_test)
  testpred <- (testans > 0) * 2 - 1
  cctest <- y_test - as.matrix(testpred)
  
  testfalseneg <- sum(cctest >0) #Number of positive samples which the model classified as negative
  testfalsepos <- sum(cctest <0) #Number of negative samples which the model classified as positive
  testerr <- (testfalseneg+testfalsepos)/nrow(D_test)
  printf("SVM TEST ERR: %f\n", testerr)
  
}

PERCENTTRAIN <- 0.8
ss <- nrow(xtrain) * PERCENTTRAIN 
trainrows <- sample(1:nrow(xtrain),ss)
D_train <- xtrain[trainrows ,]
D_test  <- xtrain[-trainrows,]
y_train <- trainclass[trainrows,]
y_test  <- trainclass[-trainrows,]
  
supportvector(D_train, D_test, y_train, y_test)
fishermethod(D_train, D_test, y_train, y_test)
