# SOFT MPBART function

# packages
library(caret)
library(MASS)
library(MCMCpack)
library(TruncatedNormal)
library(SoftBart)

# ---------------- FUNCTION TO SAMPLE LATENT VARIABLES FROM TRUNCATED MULTIVARIATE NORMAL ---------------------

#'@description Function to sample latent variables from truncated multivariate normal distribution (TMVN)
#'
#'@param mu The Kx1 mean vector of the TMVN distribution
#'@param Sigma The KxK covariance matrix of the TMVN distribution
#'@param y_i The observed class label (scalar)
#'@param K Dimension of the latent vector to be sampled
#'@return A draw of the latent vector from a TMVN distribution
sample_latent_variables <- function(mu, # mean vector Kx1
                                      Sigma, # covariance matrix KxK
                                      y_i, # observed class label 
                                      K # dimension of latent vector to be sampled
) {
  
  # case: y_i is the reference category
  if (y_i == 0) {
    # then all z_ik < 0 for k = 1,...,K
    A <- diag(K)
    b <- rep(0, K)
  } else {
    # create constraints
    A <- matrix(NA_real_, nrow = K, ncol = K)
    b <- rep(0, K)
    
    # latent utility corresponding to observed y_i has to be larger then all others
    # constraints z_j - z_{y_i} < 0 for j != y_i
    row <- 1
    for (j in 1:K) {
      if (j != y_i) {
        A[row, y_i] <- -1
        A[row, j] <- 1
        row <- row + 1
      }
    }
    # constraint maximum utility should be larger than 0, -z_{y_i} < 0
    A[row, y_i] <- -1
    b <- rep(0, nrow(A))
  }
  
  # sample and return
  z_i <- rtmvnorm(n = 1, mu = mu, Sigma = Sigma, D = A, lower = rep(-Inf, K), upper = b)
  return(drop(z_i))
}

# ------------- SOFT MPBART FUNCTION -----------------

#'@description Implementation for soft MPBART
#'
#'@param y_train Vector of training data - class labels
#'@param X_train Matrix of training data - covariates
#'@param y_test Vector of test data - class labels
#'@param X_test Matrix of test data - covariates
#'@param num_burnin Number of burn-in iterations for the sampler
#'@param num_sim Number of simulation iterations for the sampler, excluding burn-in
#'@param num_trees Number of trees used in each sum-of-trees model
#'@param K Number of class labels - 1 (dimension of latent vector)
#'@param seed Seed used for controlling randomness
#'@return The following objects are returned:
#' \item{z_draws}{MCMC draws of latent variables from truncated normals for training observations}
#' \item{mu_train_draws}{MCMC draws of predictions of latent variables for training observations}
#' \item{mu_test_draws}{MCMC draws of predictions of latent variables for test observations}
#' \item{Sigma_draws}{MCMC draws of Sigma}
soft_mpbart <- function(y_train, # training data - outcomes
                        X_train, # training data - covariates
                        y_test, # test data - outcomes
                        X_test, # test data - covariates
                        K, # number of outcome categories - 1 (dim of latent vector)
                        seed, # seed used in set.seed to control randomness. If not passed, function will generate random seed
                        # parameters
                        num_burnin = 5000, # number of burn-in iterations
                        num_sim = 5000, # number of simulation iterations (excl. burn-in)
                        num_trees = 200, # number of trees in the sum-of-trees model
                        alpha = 1, # controls sparsity in Dirichlet prior
                        beta = 2, # branching process prior - penalize depth
                        gamma = 0.95, # branching process prior - penalize new nodes
                        e = 2, # leaf node param prior - controls prior variance
                        sigma_hat = NULL, # initial guess, will be estimated via simple linear regression
                        shape = 1, # shape parameter of gating probabilities????
                        width = 0.1, # bandwidth of gating probabilities
                        alpha_scale = NULL, # scale of hyperprior on alpha
                        alpha_shape_1 = 0.5, # shape 1 of hyperprior on alpha
                        alpha_shape_2 = 1, # shape 2 of hyperprior on alpha
                        tau_rate = 10, # rate of exponential prior on bandwidth parameters tau
                        quiet = FALSE # whether you want progress bar to show
) {
  
  # set some values
  num_obs_train = length(y_train)
  num_obs_test = length(y_test)
  p <- ncol(X_train)
  
  # set seed for reproduceability
  if (is.null(seed)) {
    seed <- as.integer(Sys.time()) %% .Machine$integer.max # pseudo-random seed if not provided
  }
  set.seed(seed = seed)
  
  # ------------ INITIALIZE PARAMS ----------------
  
  # create parameters and tree samplers for softBART, and initialize lists and matrices to store draws, predictions and errors
  z <- matrix(NA_real_, nrow = num_obs_train, ncol = K) # to store sampled latent variables (resampled each iteration)
  Sigma <- diag(K) # set Sigma to identity matrix
  
  nu_prior <- K + 1 # prior d.o.f inv-Wishart
  scalematr_prior <- nu_prior * diag(K) # prior scale matrix inv-Wishart
  
  hypers <- vector("list", K) # to store hyperparams for softBART models
  tree_samplers <- vector("list", K) # to store softBART models
  predictions_z_train <- matrix(NA_real_, nrow = num_obs_train, ncol = K) # to store training predictions of latent variables (re-generated each iteration)
  predictions_z_test <- matrix(NA_real_, nrow = num_obs_test, ncol = K) # to store test predictions of latent variables (re-generated each iteration)
  
  for (k in 1:K) {
    hypers[[k]] <- Hypers(X = X_train, # not used, does not matter 
                          y = y_train, # not used, does not matter
                          alpha = alpha, 
                          beta = beta, 
                          gamma = gamma, 
                          k = e, 
                          sigma_hat = sigma_hat,
                          shape = shape, 
                          width = width,
                          num_tree = num_trees,
                          alpha_scale = alpha_scale,
                          alpha_shape_1 = alpha_shape_1,
                          alpha_shape_2 = alpha_shape_2,
                          tau_rate = tau_rate
    )
    
    tree_samplers[[k]] <- MakeForest(hypers = hypers[[k]], opts = opts, warn = FALSE)
    
    # initialize predictions for latent variables
    predictions_z_train[,k] <- tree_samplers[[k]]$do_predict(X_train)
  }
  
  z_draws <- array(NA, dim = c(num_sim, num_obs_train, K)) # to store MCMC draws of latent variables
  mu_train_draws <- array(NA, dim = c(num_sim, num_obs_train, K)) # to store MCMC draws of training predictions
  mu_test_draws <- array(NA, dim = c(num_sim, num_obs_test, K)) # to store MCMC draws of test predictions
  Sigma_draws <- array(NA, dim = c(num_sim, K, K)) # to store MCMC draws of Sigma
  
  errors <- matrix(NA_real_, num_obs_train, K) # to store errors (recalculated each iteration)
  
  opts <- Opts(num_print = num_burnin + num_sim + 1) # no need to print here, rest are default settings
  
  # ------- MCMC --------
  
  # track progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = num_burnin + num_sim, style = 3)
  }
  
  # Gibbs sampler
  for (iter in 1:num_burnin+num_sim) {
    
    # sample latent variables from truncated multivariate normal
    z <- t(sapply(1:num_obs_train, function(i) {
      sample_latent_variables(mu = colMeans(predictions_z_train), Sigma = Sigma, y_i = y_train[i], K = K)
    }))
    
    # sample all tree model related parameters using softBART package, and generate predictions
    for (k in 1:K) {
      
      # compute input z for tree sampler
      temp_z <- z[, k] - (z[, -k] - predictions_z_train[, -k]) %*% solve(Sigma[-k, -k], Sigma[-k, k])
      
      # update sigma of the tree model
      tree_samplers[[k]]$set_sigma(sqrt(Sigma[k,k] - Sigma[k,-k]%*%(solve(Sigma[-k,-k]) %*% Sigma[-k,k]))) 
      
      # predict k'th component of z (training data)
      mu_train <- t(tree_samplers[[k]]$do_gibbs(X_train, temp_z, X_train, i = 1)) 
      
      # compute errors 
      errors[,k] <- z[,k] - mu_train
      
      # save predictions
      predictions_z_train[,k] <- mu_train 
      
      # generate test predictions for test observations and store them
      mu_test <- tree_samplers[[k]]$do_predict(X_test)
      predictions_z_test[, k] <- mu_test
    }
    
    # sample unconstrained Sigma from inverted-Wishart
    nu_posterior <- nu_prior + num_obs_train
    rss <- t(errors) %*% errors
    scalematr_posterior <- scalematr_prior + rss
    Sigma_star <- riwish(nu_posterior, scalematr_posterior)
    
    # scale to force trace restriction
    Sigma <- (K / sum(diag(Sigma_star))) * Sigma_star
    
    # save draws if after burn-in
    if (iter > num_burnin) {
      iter_min_burnin <- iter - num_burnin
      
      z_draws[iter_min_burnin,,] <- z
      mu_train_draws[iter_min_burnin,,] <- predictions_z_train
      mu_test_draws[iter_min_burnin,,] <- predictions_z_test
      Sigma_draws[iter_min_burnin,,] <- Sigma
    }
    
    # update progress bar
    if (!quiet) {
      setTxtProgressBar(pb, i)
    }
    
  } # end of sampler
  
  # close progress bar
  if (!quiet) {
    close(pb)
  }
  cat("\nGibbs sampling finished.\n")
  
  # return MCMC output
  return_list <- list()
  
  return_list$z_draws <- z_draws
  return_list$mu_train_draws <- mu_train_draws
  return_list$mu_test_draws <- mu_test_draws
  return_list$Sigma_draws <- Sigma_draws
  
  return(return_list)
}

# -------- FUNCTION TO GET POSTERIOR CLASS PROBABILITIES AND PREDICT REPONSES ----------

#'@description Function to compute posterior class probabilities and generate predictions for the response.
#' Assumes 0-based class indexing where 0 is the reference class.
#'
#'@param predictions_z 3D array of predicted latent variables from S-MPBART MCMC. Dimension (num_sim, num_obs, dim_z).
#'@return The following objects:
#'\item{post_probs}{Matrix of posterior class probabilities for each observation}
#'\item{pred_y}{Vector of predicted classes for the response}
soft_mpbart_predict <- function(predictions_z) {
  
  # extract dimensions
  num_sim <- dim(predictions_z)[1]
  num_obs <- dim(predictions_z)[2]
  num_classes <- dim(predictions_z)[3] + 1 # add one as latent vector  has dimension (num_classes - 1)
  
  class_counts <- matrix(0, nrow = num_obs, ncol = num_classes)
  
  for (iter in 1:num_sim) {
    pred_z_iter <- predictions_z[iter,,]
    
    # append reference category
    pred_z_iter_full <- cbind(0, pred_z_iter)
    
    # predict response for all observations
    pred_class <- apply(pred_z_iter_full, 1, which.max) - 1 # subtract one to go from 1-based to 0-based, as 0 is reference category.
  
    # update class counts
    for (i in 1:num_obs) {
      class_counts[i, pred_class[i] + 1] <- class_counts[i, pred_class[i] + 1] + 1
    }
  }
  
  # convert to posterior probabilities
  post_probs <- class_counts / num_sim
  colnames(post_probs) <- paste0("class_", 0:(n_classes - 1)) # add class labels
  
  # compute predicted class labels
  pred_y <- max.col(post_probs) - 1 # 0-based class indexing so subtract one
  
  # return posterior probabilities and predicted class labels
  return_list <- list()
  return_list$post_probs <- post_probs
  return_list$pred_y <- pred_y
  
  return(return_list)
}

# ------ TEST ERROR RATE FUNCTION -------

#'@description Function to compute test error rate
#'
#'@param y_actual Observed class labels of the test set
#'@param y_pred Predicted class labels of the test set
test_error_rate <- function(y_actual, y_pred) {
  return(mean(y_actual != y_pred))
}





