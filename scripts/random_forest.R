# RANDOM FOREST

library(randomForest)

#'@description Function which uses Random Forest to predict a multiclass response
#'
#'@param X_train Matrix of training data - predictors
#'@param y_train Vector of training data - responses
#'@param X_test Matrix of test data - predictors
#'@param mtry Specifies the number of predictors randomly selected at each split in the decision trees
#'@return A list containing the predicted class labels on the test set, and a matrix containing predicted class probabilities
rf_multiclass <- function(X_train, y_train, X_test, mtry) {
  
  y_train <- as.factor(y_train)
  train_data <- data.frame(X_train, y = y_train)
  
  # fit rf
  rf_model <- randomForest(y ~ ., data = train_data, mtry = mtry)
  
  # predict
  predictions <- predict(rf_model, newdata = X_test)
  probs <- predict(rf_model, newdata = X_test, type = "prob")
  
  # convert factor levels to original integer class labels (assuming 0-based)
  predictions_numeric <- as.numeric(predictions) - 1
  
  # store return object
  return_list <- list()
  return_list$pred_y <- predictions_numeric
  return_list$probs <- probs
  
  return(return_list)
}