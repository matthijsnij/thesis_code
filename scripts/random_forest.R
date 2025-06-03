# RANDOM FOREST

library(randomForest)
library(caret)

# ---------- RANDOM FOREST -------------
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
  return_list$post_probs <- probs
  
  return(return_list)
}

# ------- RANDOM FOREST CV -----------
#'@description Cross validation for random forest classifier
#'
#'@param X_train Matrix of training data - predictors
#'@param y_train Vector of training data - responses
#'@param X_test Matrix of test data - predictors
#'@param mtry_grid CV grid of mtry values
#'@param K Number of folds in the CV
#'@param seed Seed to control randomness of folds
rf_multiclass_cv <- function(X_train, y_train, X_test, mtry_grid, K = 5, seed = 123) {
  set.seed(seed)
  y_train <- as.factor(y_train)
  folds <- createFolds(y_train, k = K, list = TRUE, returnTrain = FALSE)
  
  accuracy_results <- data.frame(mtry = mtry_grid, accuracy = NA)
  
  for (i in seq_along(mtry_grid)) {
    mtry_val <- mtry_grid[i]
    acc_vec <- numeric(K)
    
    for (k in seq_along(folds)) {
      test_idx <- folds[[k]]
      train_idx <- setdiff(seq_along(y_train), test_idx)
      
      X_tr <- X_train[train_idx, , drop = FALSE]
      y_tr <- y_train[train_idx]
      X_val <- X_train[test_idx, , drop = FALSE]
      y_val <- y_train[test_idx]
      
      result <- rf_multiclass(X_tr, y_tr, X_val, mtry = mtry_val)
      pred_y <- result$pred_y
      y_val_num <- as.numeric(y_val) - 1
      
      acc_vec[k] <- mean(pred_y == y_val_num)
    }
    
    accuracy_results$accuracy[i] <- mean(acc_vec)
  }
  
  # get best mtry
  best_mtry <- accuracy_results$mtry[which.max(accuracy_results$accuracy)]
  
  # train final model and predict on test set
  final_result <- rf_multiclass(X_train, y_train, X_test, mtry = best_mtry)
  final_result$best_mtry <- best_mtry
  
  return(final_result)
}