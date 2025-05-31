# SOFT MPBART function

# packages
library(caret)
library(coda)
library(MASS)
library(MCMCpack)
library(tmvmixnorm)
library(SoftBart)
library(parallel)

# ---------------- FUNCTION FOR COVARIATE NORMALIZATION ----------------------

#'@description applies rank normalization to a matrix of doubles 
#'
#'@param X matrix of doubles
#'@return Normalized matrix of doubles
rank_normalize <- function(X) {
  apply(X, 2, function(col) {
    ranks <- rank(col, ties.method = "average")
    (ranks - 1) / (length(ranks) - 1)  # maps to [0, 1]
  })
}

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
    lower <- rep(-20, K)
    upper <- rep(0, K)
  } else {
    # create constraints
    A <- matrix(0, nrow = K, ncol = K)
    lower <- rep(-20, K)
    upper <- rep(0, K)
    
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
    # constraint -z_{y_i} < 0
    A[row, y_i] <- -1
  }
  
  # sample and return
  z_i <- tmvmixnorm::rtmvn(n = 1, Mean = mu, Sigma = Sigma, lower = lower, upper = upper, D = A)
  return(drop(z_i))
}

# ------ FUNCTION TO SAFELY SAMPLE FROM INV-WISHART

#'@description Function to sample from inverted-Wishart, avoiding singular scale matrix
#'
#'@param nu Posterior d.o.f.
#'@param scale Posterior scale matrix 
#'@param jitter_init Initial jitter factor
#'@param jitter_mult Multiplier for jitter factor
#'@param max_attempts Maximum number of jitter increases before stopping
#'@return A draw from inverted-Wishart distribution
safe_riwish <- function(nu, scale, jitter_init = 1e-8, jitter_mult = 10, max_attempts = 5) {
  # check 
  if (any(!is.finite(scale))) stop("Scale has non-finite entries before jittering.")
  if (any(is.na(scale))) stop("Scale has NA entries before jittering.")
  
  # check positive definiteness 
  is_pos_def <- function(mat) {
    tryCatch({
      chol(mat)
      TRUE
    }, error = function(e) FALSE)
  }
  
  jitter <- jitter_init
  scale_jittered <- scale
  diag_n <- diag(nrow(scale))
  
  attempts <- 0
  while (!is_pos_def(scale_jittered) && attempts < max_attempts) {
    scale_jittered <- scale + jitter * diag_n
    jitter <- jitter * jitter_mult
    attempts <- attempts + 1
  }
  
  if (!is_pos_def(scale_jittered)) {
    stop("Scale matrix is not positive definite even after jittering.")
  }
  
  # sample from inverse Wishart
  Sigma_star <- riwish(nu, scale_jittered)
  return(Sigma_star)
}

# ------------- SOFT MPBART FUNCTION -----------------

#'@description Implementation for soft MPBART
#'
#'@param y_train Vector of training data - class labels
#'@param X_train Matrix of training data - covariates
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
                        X_test, # test data - covariates
                        num_classes, # number of outcome classes
                        seed = NULL, # seed used in set.seed to control randomness
                        # parameters
                        num_burnin = 1000, # number of burn-in iterations
                        num_sim = 1000, # number of simulation iterations (excl. burn-in)
                        num_trees = 200, # number of trees in the sum-of-trees model
                        alpha = 1, # controls sparsity in Dirichlet prior
                        beta = 2, # branching process prior - penalize depth
                        gamma = 0.95, # branching process prior - penalize new nodes
                        e = (1/3), # leaf node param prior - controls prior variance
                        sigma_hat = NULL, # initial guess, will be estimated via simple linear regression
                        shape = 1, # shape parameter of gating probabilities????
                        width = 0.1, # bandwidth of gating probabilities
                        alpha_scale = NULL, # scale of hyperprior on alpha
                        alpha_shape_1 = 0.5, # shape 1 of hyperprior on alpha
                        alpha_shape_2 = 1, # shape 2 of hyperprior on alpha
                        tau_rate = 10, # rate of exponential prior on bandwidth parameters tau
                        quiet = FALSE # whether you want progress bar to show
) {
  
  # set seed for reproduceability if passed
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  # set some values
  K <- num_classes - 1 # dimension latent vector
  num_obs_train = length(y_train)
  num_obs_test = nrow(X_test)
  p <- ncol(X_train)
  
  # ------------ INITIALIZE PARAMS ----------------
  
  # create parameters and tree samplers for softBART, and initialize lists and matrices to store draws, predictions and errors
  z <- matrix(NA_real_, nrow = num_obs_train, ncol = K) # to store sampled latent variables (resampled each iteration)
  Sigma <- diag(K) # set Sigma to identity matrix
  
  nu_prior <- K + 1 # prior d.o.f inv-Wishart
  scalematr_prior <- nu_prior * diag(K) # prior scale matrix inv-Wishart
  
  hypers <- vector("list", K) # to store hyperparams for softBART models
  tree_samplers <- vector("list", K) # to store softBART models
  predictions_z_train <- matrix(0, nrow = num_obs_train, ncol = K) # to store training predictions of latent variables (re-generated each iteration)
  predictions_z_test <- matrix(0, nrow = num_obs_test, ncol = K) # to store test predictions of latent variables (re-generated each iteration)
  
  opts <- Opts(num_print = num_burnin + num_sim + 1) # no need to print here, rest are default settings
  
  for (k in 1:K) {
    hypers[[k]] <- Hypers(X = X_train, # not used, does not matter 
                          Y = y_train, # not used, does not matter
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
  
  errors <- matrix(0, num_obs_train, K) # to store errors (recalculated each iteration)
  
  # ------- MCMC --------
  
  # track progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = num_burnin + num_sim, style = 3)
  }
  
  # Gibbs sampler
  for (iter in 1:(num_burnin+num_sim)) {
    
    if (any(!is.finite(Sigma)) || any(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values <= 1e-8)) {
      stop(paste("Numerical issue with Sigma at iteration", iter))
    }
    
    for (i in 1:num_obs_train) {
      mu_i <- predictions_z_train[i, ]
      
      if (any(!is.finite(mu_i))) {
        stop(paste("Non-finite mu at i =", i, "iter =", iter))
      }
      if (any(is.na(mu_i))) {
        stop(paste("NA mu at i =", i, "iter =", iter))
      }
      
      z_i <- sample_latent_variables(mu = mu_i, Sigma = Sigma, y_i = y_train[i], K = K)
      
      if (any(!is.finite(z_i))) {
        stop(paste("Non-finite z for observation", i, "at iteration", iter))
      }
      z[i, ] <- z_i
    }
    
    # sample all tree model related parameters using softBART package, and generate predictions
    for (k in 1:K) {
      
      # compute input z for tree sampler
      temp_z <- z[, k] - (z[, -k, drop = FALSE] - predictions_z_train[, -k, drop = FALSE]) %*% 
        solve(Sigma[-k, -k, drop = FALSE], Sigma[-k, k, drop = FALSE]
      )
      
      # update sigma of the tree model
      tree_samplers[[k]]$set_sigma(
        sqrt(Sigma[k, k] - Sigma[k, -k, drop = FALSE] %*% solve(Sigma[-k, -k, drop = FALSE]) %*% Sigma[-k, k, drop = FALSE])
      )
      
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
    
    if (any(!is.finite(errors))) stop("Non-finite errors at iteration ", iter)
    if (any(!is.finite(rss))) stop("Non-finite RSS at iteration ", iter)
    
    scalematr_posterior <- scalematr_prior + rss
    Sigma_star <- safe_riwish(nu = nu_posterior, scale = scalematr_posterior)
    
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
      setTxtProgressBar(pb, iter)
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
  colnames(post_probs) <- paste0("class_", 0:(num_classes - 1)) # add class labels
  
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

# ------ MULTICLASS BRIER SCORE FUNCTION ---------

#'@description Function which computes the multiclass Brier score
#'
#'@param y_actual nx1 vector of observed class labels (0-based)
#'@param y_prob nxK matrix of predicted posterior class probabilities
#'@return The multiclass Brier score
brier_score_multiclass <- function(y_actual, y_prob) {
  
  n <- length(y_actual)
  K <- ncol(y_prob)
  
  # one-hot encode y_actual
  y_onehot <- matrix(0, nrow = n, ncol = K)
  for (i in 1:n) {
    y_onehot[i, y_actual[i] + 1] <- 1 # take 0-based class labeling into account
  }
  
  # compute Brier score
  score <- mean(rowSums((y_prob - y_onehot)^2))
  return(score)
}

# ------ GEWEKE TEST FUNCTION ---------
geweke_test_z <- function(z_draws, n_check = 3) {
  
  n_iter <- dim(z_draws)[1]
  n_obs  <- dim(z_draws)[2]
  n_comp <- dim(z_draws)[3]
  
  obs_indices <- sample(1:n_obs, min(n_check, n_obs))
  
  for (obs in obs_indices) {
    comp <- sample(1:n_comp, 1)
    
    trace <- z_draws[, obs, comp]
    mcmc_trace <- mcmc(trace)
    geweke_result <- geweke.diag(mcmc_trace)
    z_score <- geweke_result$z
    
    cat(sprintf("Obs %d, Component %d --> Geweke z-score: %.4f\n", obs, comp, z_score))
    
    # plot trace and ACF
    par(mfrow = c(1, 2))
    plot(trace, type = "l", main = sprintf("Trace Plot\nObs %d, Comp %d", obs, comp), xlab = "Iteration", ylab = "Value")
    acf(trace, main = "ACF")
  }
  
  par(mfrow = c(1, 1))  # reset layout
}





