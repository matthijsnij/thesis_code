# SIMULATION FUNCTIONS

library(glue)
library(openxlsx)
library(MASS)
library(fields)
source("soft_mpbart.R")

# -------- DGP 1 DATA GENERATION -------

#'@description Function which generates data for a single replication of DGP 1
#'
#'@param n_train Number of training observations to be generated
#'@param n_test Number of test observations to be generated
#'@return A list containing the following objects:
#'\item{X_train}{Matrix of predictors for training data}
#'\item{y_train}{Vector of responses for training data}
#'\item{X_test}{Matrix of predictors for test data}
#'\item{y_test}{Vector of responses for test data}
generate_dgp1_data <- function(n_train, n_test) {
  
  # helper functions: one latent function per class
  compute_z0 <- function(x1, x2, eps) {
    ifelse(x1 <= 0.5,
           5 * sin(2 * pi * x1) + 3 * x2 + eps,
           5 * sin(2 * pi * x1) - 3 * x2 + eps)
  }
  
  compute_z1 <- function(x1, x2, eps) {
    ifelse(x2 <= 0.3,
           4 * cos(pi * x2) + x1 + eps,
           8 * cos(pi * x2) - x1 + eps)
  }
  
  compute_z2 <- function(x1, x2, eps) {
    ifelse(x1 + x2 <= 1,
           3 * sin(pi * (x1 + x2)) + eps,
           7 * sin(pi * (x1 + x2)) + 2) + eps
  }
  
  # train data
  X_train <- matrix(runif(n_train * 2), ncol = 2)
  eps_train <- matrix(rnorm(n_train * 3), ncol = 3) # 3 noises, one per class
  
  z_train <- cbind(
    compute_z0(X_train[, 1], X_train[, 2], eps_train[,1]),
    compute_z1(X_train[, 1], X_train[, 2], eps_train[,2]),
    compute_z2(X_train[, 1], X_train[, 2], eps_train[,3])
  )
  
  # test data
  X_test <- matrix(runif(n_test * 2), ncol = 2)
  eps_test <- matrix(rnorm(n_test * 3), ncol = 3)
  
  z_test <- cbind(
    compute_z0(X_test[, 1], X_test[, 2], eps_test[,1]),
    compute_z1(X_test[, 1], X_test[, 2], eps_test[,2]),
    compute_z2(X_test[, 1], X_test[, 2], eps_test[,3])
  )
  
  # assign class labels as the index of max latent value per row (subtract 1 to have classes 0,1,2)
  assign_class <- function(z_matrix) {
    apply(z_matrix, 1, function(z) which.max(z) - 1)
  }
  
  y_train <- assign_class(z_train)
  y_test <- assign_class(z_test)
  
  list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test
  )
}

# --------- DGP 2 DATA GENERATION -------

#'@description Function which generates data for a single replication of DGP 2 (possibly with extra noise predictors)
#'
#'@param n_train Number of training observations to be generated
#'@param n_test Number of test observations to be generated
#'@param p Number of predictors to include. Should be at least 10. Any additional predictors are noise predictors.
#'@return A list containing the following objects:
#'\item{X_train}{Matrix of predictors for training data}
#'\item{y_train}{Vector of responses for training data}
#'\item{X_test}{Matrix of predictors for test data}
#'\item{y_test}{Vector of responses for test data}
generate_dgp2_data <- function(n_train, n_test, p) {
  
  if (p < 10) {
    stop("p should be larger or equal to 10")
  }
  
  # generate predictors
  X_train <- matrix(runif(n_train * p), nrow = n_train, ncol = p)
  X_test  <- matrix(runif(n_test  * p), nrow = n_test,  ncol = p)
  
  # define latent functions
  compute_z0 <- function(X) {
    10 * sin(pi * X[,1] * X[,2]) + 
      20 * (X[,3] - 0.5)^2 + 
      10 * X[,4] + 
      5 * X[,5]
  }
  
  compute_z1 <- function(X) {
    8 * cos(pi * X[,6] * X[,7]) + 
      15 * (X[,8] - 0.4)^2 + 
      7 * X[,9] + 
      3 * X[,10]
  }
  
  compute_z2 <- function(X) {
    12 * sin(1.5 * pi * X[,3] * X[,8]) + 
      10 * (X[,5] - 0.6)^2 + 
      6 * X[,1] + 
      8 * X[,7]
  }
  
  # generate latent variables + noise for train and test
  z_train <- cbind(
    compute_z0(X_train) + rnorm(n_train),
    compute_z1(X_train) + rnorm(n_train),
    compute_z2(X_train) + rnorm(n_train)
  )
  
  z_test <- cbind(
    compute_z0(X_test) + rnorm(n_test),
    compute_z1(X_test) + rnorm(n_test),
    compute_z2(X_test) + rnorm(n_test)
  )

  assign_class <- function(z_matrix) {
    apply(z_matrix, 1, function(z) which.max(z) - 1)
  }
  
  y_train <- assign_class(z_train)
  y_test <- assign_class(z_test)
  
  list(
    X_train = X_train,
    y_train = y_train,
    X_test  = X_test,
    y_test  = y_test
  )
}

# --------- DGP 3 DATA GENERATION --------

#'@description Function which generates data for a single replication of DGP 3 
#'
#'@param n_train Number of training observations to be generated
#'@param n_test Number of test observations to be generated
#'@param n_predictors Number of predictors to include
#'@param n_classes Number of class labels
#'@return A list containing the following objects:
#'\item{X_train}{Matrix of predictors for training data}
#'\item{y_train}{Vector of responses for training data}
#'\item{X_test}{Matrix of predictors for test data}
#'\item{y_test}{Vector of responses for test data}
generate_dgp3_data <- function(n_train, n_test,
                                  n_predictors = 8,
                                  n_classes = 3) {
  
  n_total <- n_train + n_test
  
  # latent variable z governs all structure (smooth and 1D)
  z <- runif(n_total, min = 0, max = 2 * pi)
  
  # predictors: smooth nonlinear transformations of z
  X <- matrix(0, nrow = n_total, ncol = n_predictors)
  if (n_predictors >= 1) X[, 1] <- sin(z)
  if (n_predictors >= 2) X[, 2] <- cos(2 * z)
  if (n_predictors >= 3) X[, 3] <- z^2
  if (n_predictors >= 4) X[, 4] <- exp(-z)
  if (n_predictors >= 5) X[, 5] <- sin(3 * z)
  if (n_predictors >= 6) X[, 6] <- cos(5 * z)
  if (n_predictors >= 7) X[, 7] <- log(z + 1e-2)
  if (n_predictors >= 8) X[, 8] <- sqrt(z)
  
  # class-specific latent functions (smooth functions of z)
  F_all <- matrix(0, nrow = n_total, ncol = n_classes)
  F_all[, 1] <- sin(z)
  F_all[, 2] <- cos(z)
  F_all[, 3] <- exp(- (z - pi)^2)
  
  shifts <- c(-2, 0, 2)
  F_all <- sweep(F_all, 2, shifts, "+")
  
  # Softmax function
  softmax <- function(F) {
    expF <- exp(F - apply(F, 1, max))  # for numerical stability
    expF / rowSums(expF)
  }
  P_all <- softmax(F_all)
  
  # sample class labels
  sample_labels <- function(P) {
    apply(P, 1, function(p) sample(0:(n_classes - 1), size = 1, prob = p))
  }
  
  y_all <- sample_labels(P_all)
  y_train <- y_all[1:n_train]
  y_test  <- y_all[(n_train + 1):n_total]
  
  X_train <- X[1:n_train, , drop = FALSE]
  X_test  <- X[(n_train + 1):n_total, , drop = FALSE]
  
  # rank normalize predictors
  X_train <- rank_normalize(X_train)
  X_test <- rank_normalize(X_test)
  
  result <- list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test
  )
  
  return(result)
}

# --------- RUN METHOD ON SIMULATED DATA ----------

#'@description Function which runs a method on all simulated data (all replications). Writes the results to excel.
#'
#'@param method The method to run: smpbart, mpbart or rf
#'@param sim_data List of all replications of simulated data
#'@param which_dgp The DGP from which data is generated
run_method <- function(method, sim_data, which_dgp) {
  
  if (which_dgp == "dgp1") {
    num_classes <- 3
    mtry_grid <- c(1,2)
  } else if (which_dgp == "dgp2") {
    num_classes <- 3
    mtry_grid <- c(2, 3, 4, 6, 8, 10)
  } else if (which_dgp == "dgp2extranoise") {
    num_classes <- 3
    mtry_grid <- c(2, 4, 8, 15, 30, 60)
  } else if (which_dgp == "dgp3") {
    num_classes <- 3
    mtry_grid <- c(2, 3, 4, 6, 8)
  } else {
      stop("Run with a correct DGP. Choose from 'dgp1', 'dgp2', 'dgp2extranoise', 'dgp3'")
  }
  
  # to store error rates and brier scores
  reps <- length(sim_data)
  error_rates <- rep(0, reps)
  brier_scores <- rep(0, reps)
  
  # to store output
  wb_output <- createWorkbook()
  
  for (r in 1:reps) {
    data_list <- sim_data[[r]]
    
    X_train <- data_list$X_train
    y_train <- data_list$y_train
    X_test <- data_list$X_test
    y_test <- data_list$y_test
    
    # define column names
    p <- ncol(X_train)
    colnames(X_train) <- paste0("x", 1:p)
    colnames(X_test) <- paste0("x", 1:p)
    
    # run chosen method
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
        pred_output <- mpbart(X_train = X_train,
                              y_train = y_train,
                              X_test = X_test,
                              num_burnin = 1500,
                              num_sim = 1500,
                              num_classes = num_classes
      )
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
  
  path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/{method}_{which_dgp}_output.xlsx")
  saveWorkbook(wb_output, path, overwrite = TRUE)
}

# --------- WRITING GENERATED DATA TO EXCEL ---------

#'@description Function which writes generated data to excel
#'
#'@param n_rep Number of replications for which data is generated
#'@param all_data List of all generated data. Each replication of generated data is also a list.
#'@param which_dgp String of which DGP is generated
#'@param seed The seed used for generating the data
write_data <- function(n_rep, all_data, which_dgp, seed) {
  wb <- createWorkbook()
  
  for (i in 1:n_rep) {
    rep_data <- all_data[[i]]
    
    # combine X and y for train/test
    train_df <- as.data.frame(cbind(rep_data$X_train, y = rep_data$y_train))
    test_df <- as.data.frame(cbind(rep_data$X_test, y = rep_data$y_test))
    
    # define column names
    p <- ncol(rep_data$X_train)
    colnames(train_df) <- c(paste0("x_", 1:p), "y")
    colnames(test_df)  <- c(paste0("x_", 1:p), "y")
    
    addWorksheet(wb, paste0("rep_", i, "_train"))
    addWorksheet(wb, paste0("rep_", i, "_test"))
    
    writeData(wb, paste0("rep_", i, "_train"), train_df)
    writeData(wb, paste0("rep_", i, "_test"), test_df)
  }
  path <- glue("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/{which_dgp}_data_seed={seed}.xlsx")
  saveWorkbook(wb, path, overwrite = TRUE)
}