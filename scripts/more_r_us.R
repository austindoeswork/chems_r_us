library("MASS") #Contains the Fisher LDA command 
library("readr")  #The library for read_csv
library("e1071") #SVM
# install.packages("e1071", repos="http://R-Forge.R-project.org")

# modify this for wherever your project is
#path <- 
printf <- function(...) cat(sprintf(...))

# read the data
#chemstrain <- read_csv(paste(path, "/data/chemsdata.csv", sep=""))
#chemstest <- read_csv(paste(path, "/data/chemstest.csv", sep=""))

# QUESTION 2
# get the useful stuff
#xtrain <- chemstrain[,-1]
#xtest <- chemstest[,-1]
#postrain <- xtrain[xtrain[,ncol(xtrain)] == 1,]
#negtrain <- xtrain[xtrain[,ncol(xtrain)] == -1,]
#trainclass <- chemstrain[, ncol(chemstrain)]

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
  D_train <- chemstrain
  D_test <- chemstest
  y_train <- trainclass
  y_test <- testclass
  
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

# QUESTION 2(d)
supportvector_predict <-function(D_train, D_test, y_train) {
  D_train <- chemstrain
  D_test <- testpoints
  y_train <- trainclass
  
  model <-svm(D_train, y_train)
  
  trainans <- predict(model,D_train)
  trainpred <- (trainans > 0) * 2 - 1
  cctrain <- y_train - as.matrix(trainpred)
  
  trainfalseneg <- sum(cctrain <0) #Number of positive samples which the model classified as negative
  trainfalsepos <- sum(cctrain >0) #Number of negative samples which the model classified as positive
  trainerr <- (trainfalseneg+trainfalsepos)/nrow(D_train)
  printf("SVM TRAIN ERR: %f\n", trainerr)
  
  testans <- predict(model,as.data.frame(D_test))
  testpred <- (testans > 0) * 2 - 1
  return(testpred)
}

#PERCENTTRAIN <- 0.9
#ss <- nrow(xtrain) * PERCENTTRAIN 
#trainrows <- sample(1:nrow(xtrain),ss)
#D_train <- xtrain[trainrows ,]
#D_test  <- xtrain[-trainrows,]
#y_train <- trainclass[trainrows,]
#y_test  <- trainclass[-trainrows,]
  
#supportvector(D_train, D_test, y_train, y_test)
#fishermethod(D_train, D_test, y_train, y_test)

#########################################################
#A routine to plot a pair of histograms and a threshold
histopair <-function(pminus,pplus,thresh,yy=c(0,60),label2="Plus",label1="Minus") {
  require(ggplot2); require(gridExtra)
  hist1 <- ggplot(as.data.frame(pminus), aes(pminus)) + geom_histogram(col="blue",fill="blue")
  hist2 <- ggplot(as.data.frame(pplus),  aes(pplus))  + geom_histogram(col="red",fill="red")
  df <- data.frame(x1=c(thresh,thresh),y1=c(yy[1],yy[2]))
  pmin <- min(pminus,pplus)
  pmax<- max(pminus,pplus)
  me1 <- hist1 + expand_limits(x=c(pmin,pmax)) + geom_line(data=df,aes(x=x1,y=y1)) + xlab(label1)
  me2 <- hist2 + expand_limits(x=c(pmin,pmax)) + geom_line(data=df,aes(x=x1,y=y1)) + xlab(label2)
  pl <- vector("list",2)
  pl[[1]] <- me1;  pl[[2]]<-me2;
  grid.arrange(grobs=pl,ncols=1)
}

mean_method<-function(D_train, D_test, y_train, y_test){
  
  #D_train <- chemstrain
  #D_test <- chemstest
  #y_train <- trainclass
  #y_test <- testclass
  
  #SPLIT PLUS & MINUS
  plusonerows <- D_train[y_train==1,]
  minusonerows <- D_train[y_train==-1,]
  plusonerows_test <- D_test[y_test==1,]
  minusonerows_test <- D_test[y_test==-1,]
  
  pca <- prcomp(D_train, retx=TRUE, center=TRUE)

  plot(plusonerows%*%pca$rotation)
  plot(minusonerows%*%pca$rotation)
  plot(pca)
  
  ggplot(data=as.data.frame(pca$x), mapping=aes(x=`PC6`,y=`PC7`, color=trainclass)) + geom_point() + scale_color_gradientn(colors=rainbow(5))
  
  #mean(plusonerows) - mean(minusonerows)
  #
  
  ###GET CLASS MEANS
  mean_plusone <- apply(plusonerows, 2, mean)
  mean_minusone <- apply(minusonerows, 2, mean)
  midpnt <- (mean_plusone + mean_minusone)/2
  normal <- mean_plusone - mean_minusone
  
  
  ###GET TRAIN ERROR
  proj_por <- as.matrix(plusonerows)%*%as.matrix(normal)
  proj_mor <- as.matrix(minusonerows)%*%as.matrix(normal)
  proj_mdp <- t(as.matrix(midpnt))%*%as.matrix(normal)
  
  #classification error
  error <- 0
  for(i in 1:nrow(proj_por)){
    if(proj_por[i,]>=proj_mdp){
      error <- error + 1
    }
  }
  for(i in 1:nrow(proj_mor)){
    if(proj_mor[i]<proj_mdp){
      error <- error + 1
    }
  }
  train_error = error/nrow(D_train)

  
  
  ###GET TEST ERROR
  proj_por <- as.matrix(plusonerows_test)%*%as.matrix(normal)
  proj_mor <- as.matrix(minusonerows_test)%*%as.matrix(normal)
  
  #classification error
  error <- 0
  for(i in 1:nrow(proj_por)){
    if(proj_por[i,]-proj_mdp>=0){
      error <- error + 1
    }
  }
  for(i in 1:nrow(proj_mor)){
    if(proj_mor[i]-proj_mdp<0){
      error <- error + 1
    }
  }
  test_error = error/nrow(D_test)
  
  "MEAN METHOD TRAINING ERROR: "
  train_error
  "MEAN METHOD TEST ERROR: "
  test_error
  
}

strip_outliers <- function(dataset, ridiculousness){
  
  #dataset <- chemstrain
  #ridiculousness <- 3
  
  pca <- prcomp(dataset, scale=TRUE, retx = TRUE)
  
  ok <- c(1:nrow(dataset)) > 0
  
  for(pc in 1:nrow(pca$rotation)){
    for(i in 1:nrow(dataset)){
      if(pca$x[i,pc] > ridiculousness*pca$sdev[pc]){
        ok[i] <- FALSE
      }
    }
  }
  clean_dataset <- dataset[ok==TRUE,]
  return(clean_dataset)
}

nrow(chemstrain)
nrow(strip_outliers(chemstrain, 5))


######################################################################
library("MASS") #Contains the Fisher LDA command
library("readr")  #The library for read_csv


##################################################################
#Change the path to your folder or use rstudio to read in chemsdata
chemsdata <- read_csv("chemsdata.csv")

#################################################################
#Split the chemsdata into training and testing sets and classification columns and features
#ss will be the number of data in the training set
ss<- 950
#Caution: This next command will create a new random sample every time it is executed
train <- sample(1:nrow(chemsdata),ss)

#The first column is just the sample label, we can discard it for statics (This is the -1 in the second position)
chemstrain <- chemsdata[train , -1]  #The training data is just the training rows
chemstrain <- strip_outliers(chemstrain, 5)
chemstest <- chemsdata[-train, -1 ]  # Using -train gives us all rows except the training rows.

#The last column is the class label, which is 1 or -1, we can split it off as the class
trainclass <- chemstrain[, ncol(chemstrain)]
testclass <- chemstest[,ncol(chemstest)]

#We leave off the last column (the class label) to get the matrices of features
trainmatrix <- as.matrix(chemstrain[ ,c(1:ncol(chemstrain)-1)])
testmatrix <- as.matrix(chemstest[ ,c(1:ncol(chemstest)-1)])


#########################################################################
#The Fisher LDA command is as follows
#The classification column in chemstrain is called "class" which is used in the first argument. 
# The "prior" option specifies weighting between classes. This uses (1/2,1/2) saying they are weighted only by size.
z <- lda(class ~ .,chemstrain,prior=c(1,1)/2)

plusonerows <- chemstrain[trainclass==1,]
minusonerows <- chemstrain[trainclass==-1,]






#This returns z$scaling which is the vector normal to the separating hyperplane
#It also returns z$means which are the means of the two different classes. 
#Calculate the Fisher threshold from the means and the normal vector.
thresh <- ((z$means[1,] + z$means[2,])/2)%*%z$scaling


#########################################################################
#Plot the results of the training set
#z$scaling points in the normal of the hyperplane. 
#Project the training data onto this vector 
proj <- testmatrix%*%as.matrix(z$scaling)
pplus  <- proj[testclass[ ]>0] #All the class 1 projections
pminus <- proj[testclass[ ]<0] #All the class -1 projections
#Using the histopair command defined above make histogram plots
histopair(pminus,pplus,thresh,label1="Not Readily Biodegradable",label2="Readily Biodegradable") 

#########################################################################
#Find out how well the calculated hyperplane clasifies the training data
#The predict command understands the structure returned by lda and creates a list of classes predicted for the second argument
#In this case, we see how well we did with the training data
ans <- predict(z,chemstest)$class

#Find how many samples were misclassified in training data by subtracting ans from trainclass
cc <- testclass - as.matrix(as.numeric(as.matrix(ans)))
plusasminus <- sum(cc >0) #Number of positive samples which the model classified as negative
minusasplus <- sum(cc <0) #Number of negative samples which the model classified as positive

#Show some of the numbers and percentages
plusasminus 
minusasplus 
plusasminus/nrow(chemstest)  
minusasplus/nrow(chemstest)
(plusasminus+minusasplus)/nrow(chemstest)  #Percentage of total misclassified training points


################################################################## 2.D
#Change the path to your folder or use rstudio to read in chemsdata
chemsdata <- read_csv("chemsdata.csv")

#################################################################
#Split the chemsdata into training and testing sets and classification columns and features
#ss will be the number of data in the training set
ss <- 950
#Caution: This next command will create a new random sample every time it is executed
train <- sample(1:nrow(chemsdata),ss)

#The first column is just the sample label, we can discard it for statics (This is the -1 in the second position)
chemstrain <- chemsdata[train , -1]  #The training data is just the training rows
chemstrain <- strip_outliers(chemstrain, 10)
chemstest <- chemsdata[-train, -1 ]  # Using -train gives us all rows except the training rows.

#The last column is the class label, which is 1 or -1, we can split it off as the class
trainclass <- chemstrain[, ncol(chemstrain)]
testclass <- chemstest[,ncol(chemstest)]

testpoints <- read_csv("chemstest.csv")
testpoints <- testpoints[ ,c(1:ncol(testpoints)-1)]
testpoints
supportvector(chemstrain, chemstest, trainclass, testclass)
prediction <- supportvector_predict(chemstrain, chemstest, trainclass)
prediction


