library(e1071)

create_classifier_from_continuous <- function(column){
  avg <- mean(column)
  classified <- column
  for(i in 1:length(classified)){
    if(classified[i]>=avg){
      classified[i] = 1
    }else{
      classified[i] = -1
    }
  }
  return (classified)
}

scale_continuous <- function(column){
  avg <- mean(column)
  s <- sd(column)
  classified <- column
  for(i in 1:length(classified)){
    classified[i] = as.double(classified[i] - avg)/s
  }
  str(classified)
  return (classified)
}


run_svm_hyperplane <- function(predictors, target, cost){
  model <- svm(predictors, target, cost=cost)
  return (model)
}

run_svm_classifier <- function(predictors, target, cost){
    model <- svm(predictors, target, cost=cost)
    predictedTarget <- predict(model, predictors)
    return (predictedTarget)
