# DATA APPLICATION FUNCTIONS

library(glue)
library(openxlsx)

# --------- READ DATA ---------

#'@description Function which reads the pre-processed data for a given data set and returns it as a list
#'
#'@param dataset The dataset to read
#'@return A list containing:
#'\item{X}{Matrix of normalized covariates}
#'\item{y}{Vector of observed class labels}
read_data <- function(dataset) {
  
  all_datasets <- c("glass", "vertebral", "iris")
  
  if (!(dataset %in% all_datasets)) {
    stop("Please input a correct dataset. Choose from 'glass', 'vertebral', 'iris'")
  }
  
  path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/{dataset}_preprocessed.csv")
  data <- read.csv(path, header = TRUE)
  y <- data[[ncol(data)]]
  X <- as.matrix(data[, 1:(ncol(data) - 1)])
  
  # return in a list
  list(
    X = X,
    y = y
  )
  
}

# --------- RUN METHOD ON TRAIN/TEST FOLDS OF DATA ----------

#'@description Function which creates train/test folds and runs a method on these folds. Writes the results to excel.
#'
#'@param method The method to run: smpbart, mpbart or rf
#'@param data The preprocessed data (covariates + class labels)
#'@param which_dataset The name of the dataset to run on
#'@param seed Seed for controlling randomness of the folds
run_method <- function(method, data, which_dataset, seed) {
  
  if (which_dataset == "glass") {
    num_classes <- 6
    num_folds <- 5
    mtry_grid <- c(2, 3, 4, 6, 7, 9)
  } else if (which_dataset == "vertebral") {
    num_classes <- 3
    num_folds <- 10
    mtry_grid <- c(2, 3, 4, 6)
  } else if (which_dataset == "iris") {
    num_classes <- 3
    num_folds <- 10
    mtry_grid <- c(1, 2, 3, 4)
  } else {
    stop("Run with a correct dataset. Choose from 'glass', 'vertebral', 'iris', ''")
  }
  
  # to store error rates and brier scores
  error_rates <- rep(0, num_folds)
  brier_scores <- rep(0, num_folds)
  
  # create train and test folds
  set.seed(seed = seed)
  y <- data$y
  X <- data$X
  folds <- createFolds(y, k = num_folds, list = TRUE, returnTrain = TRUE)
  
  # to store output
  wb_output <- createWorkbook()
  
  for (i in seq_along(folds)) {
    
    # get correct train and test data
    y_train <- y[folds[[i]]]
    X_train <- as.matrix(X[folds[[i]], ])
    y_test <- y[-folds[[i]]]
    X_test <- as.matrix(X[-folds[[i]], ])
    
    if (method == "smpbart") {
      # run mcmc
      mcmc_output <- soft_mpbart(y_train = y_train,
                                 X_train = X_train,
                                 X_test = X_test,
                                 num_classes = num_classes,
                                 num_burnin = 1500,
                                 num_sim = 1500
      )
      pred_output <- soft_mpbart_predict(predictions_z = mcmc_output$mu_test_draws)
    } else if (method == "mpbart") {
      stop()
    } else if (method == "rf") {
      pred_output <- rf_multiclass_cv(X_train = X_train, y_train = y_train, X_test = X_test, mtry_grid = mtry_grid)
    } else {
      stop("Please run with a correct method. Choose from ['smpbart', 'mpbart', 'rf']")
    }
    
    # compute test error rate and brier score
    error <- test_error_rate(y_actual = y_test, y_pred = pred_output$pred_y)
    brier_score <- brier_score_multiclass(y_actual = y_test, y_prob = pred_output$post_probs)
    error_rates[r] <- error
    brier_scores[r] <- brier_score
  }
  
  # save output
  addWorksheet(wb_output, "misclassification_rates")
  addWorksheet(wb_output, "brier_scores")
  writeData(wb_output, sheet = "misclassification_rates", x = error_rates)
  writeData(wb_output, sheet = "brier_scores", x = brier_scores)
  
  path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/{method}_{which_dataset}_output.xlsx")
  saveWorkbook(wb_output, path, overwrite = TRUE)
}