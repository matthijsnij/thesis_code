# SIMULATION FUNCTIONS

library(glue)
library(openxlsx)

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
  # helper function for generating y
  compute_z <- function(x1, x2, epsilon) {
    ifelse(x1 <= 0.5,
           5 * sin(2 * pi * x1) + 3 * x2 + epsilon,
           5 * sin(2 * pi * x1) - 3 * x2 + epsilon)
  }
  
  # train data
  X_train <- matrix(runif(n_train * 2), ncol = 2)
  eps_train <- rnorm(n_train)
  z_train <- compute_z(X_train[, 1], X_train[, 2], eps_train)
  
  # define quantile cut points based on z_train
  q1 <- quantile(z_train, probs = 1/3)
  q2 <- quantile(z_train, probs = 2/3)
  
  y_train <- cut(z_train, breaks = c(-Inf, q1, q2, Inf), labels = c(0, 1, 2))
  y_train <- as.integer(as.character(y_train))
  
  # test data
  X_test <- matrix(runif(n_test * 2), ncol = 2)
  eps_test <- rnorm(n_test)
  z_test <- compute_z(X_test[, 1], X_test[, 2], eps_test)
  
  y_test <- cut(z_test, breaks = c(-Inf, q1, q2, Inf), labels = c(0, 1, 2))
  y_test <- as.integer(as.character(y_test))
  
  list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test
  )
}

# --------- DGP 2 DATA GENERATION

#'@description Function which generates data for a single replication of DGP 2 (possibly with extra noise predictors)
#'
#'@param n_train Number of training observations to be generated
#'@param n_test Number of test observations to be generated
#'@param p Number of predictors to include. Should be at least 5. Any additional predictors are noise predictors.
#'@return A list containing the following objects:
#'\item{X_train}{Matrix of predictors for training data}
#'\item{y_train}{Vector of responses for training data}
#'\item{X_test}{Matrix of predictors for test data}
#'\item{y_test}{Vector of responses for test data}
generate_dgp2_data <- function(n_train, n_test, p) {
  
  if (p < 5) {
    stop("p should be larger or equal to 5")
  }
  
  # generate predictors
  X_train <- matrix(runif(n_train * p), nrow = n_train, ncol = p)
  X_test <- matrix(runif(n_test * p), nrow = n_test, ncol = p)
  
  # latent function for Friedman #1
  latent_fun <- function(X) {
    10 * sin(pi * X[,1] * X[,2]) + 
      20 * (X[,3] - 0.5)^2 + 
      10 * X[,4] + 
      5 * X[,5]
  }
  
  # generate latent variable with noise
  z_train <- latent_fun(X_train) + rnorm(n_train)
  z_test <- latent_fun(X_test) + rnorm(n_test)
  
  # compute quantile thresholds based on train latent variable
  q1 <- quantile(z_train, probs = 1/3)
  q2 <- quantile(z_train, probs = 2/3)
  
  # compute classes 0,1,2 based on thresholds
  y_train <- as.numeric(as.character(cut(z_train, breaks = c(-Inf, q1, q2, Inf), labels = c(0, 1, 2))))
  y_test  <- as.numeric(as.character(cut(z_test,  breaks = c(-Inf, q1, q2, Inf), labels = c(0, 1, 2))))
  
  list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test
  )
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
    mtry_grid <- c(2, 4, 8, 15, 30, 50)
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