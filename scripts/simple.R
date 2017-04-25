library("readr")  #The library for read_csv

# modify this for wherever your project is
path <- "/Users/austin/RPI/idm/proj/chems_r_us"

# read the data
chemstrain <- read_csv(paste(path, "/data/chemsdata.csv", sep=""))
chemstest <- read_csv(paste(path, "/data/chemstest.csv", sep=""))

# get the useful stuff
xtrain <- chemstrain[,-1]
xtest <- chemstest[,-1]
postrain <- xtrain[xtrain[,ncol(xtrain)] == 1,]
negtrain <- xtrain[xtrain[,ncol(xtrain)] == -1,]

# QUESTION 2(b)
# boxplot(xtrain)

# QUESTION 2(c)
meanmethod <-function(D, percenttrain = 0.9) {
  #Dtrain = 
}