# SOFT MPBART


# packages
library(caret)
library(MASS)
library(MCMCpack)
library(TruncatedNormal)
library(SoftBart)

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

# ------------ READ DATA ---------------

glass_data <- read.csv('C:\Users\matth\OneDrive\Bureaublad\msc_thesis\Data\glass\glass.data', header = FALSE)
glass_y <- glass_data[[ncol(glass_data)]]
glass_X <- as.matrix(glass_data[, 2:(ncol(glass_data)-1)])
glass_X_norm <- rank_normalize(glass_X)

vertebral_data <- read.csv('C:\Users\matth\OneDrive\Bureaublad\msc_thesis\Data\vertebral\data.txt', header = FALSE)
vertebral_y <- glass_data[[ncol(glass_data)]]
vertebral_X <- as.matrix(glass_data[, 1:(ncol(glass_data)-1)])
vertebral_X_norm <- rank_normalize(vertebral_X)


# ------------ PREPROCESS DATA -------------

# glass
# clean the class labels such that they fall in the range [0,5]
# there is no class 4 in the data set
for (i in 1:length(glass_y)) {
  if (glass_y[i] < 4) {
    glass_y[i] <- glass_y[i] - 1
  } else {
    glass_y[i] <- glass_y[i] - 2
  }
}

# vertebral
# change class labels to 0 = Hernia, 1 = Spondylolisthesis, 2 = Normal
for (i in 1:length(vertebral_y)) {
  if (vertebral_y[i] == "Hernia") {
    vertebral_y[i] <- 0
  } else if (vertebral_y[i] == "Spondylolisthesis") {
    vertebral_y[i] <- 1
  } else {
    vertebral_y[i] <- 2
  }
}


# ------------ CREATE TRAIN AND TEST SETS

# generate train and test set folds
num_folds <- 5
folds <- createFolds(glass_y, k = num_folds, list = TRUE, returnTrain = TRUE)

# split the data based on these indices
glass_y_train <- glass_y[folds[[1]]]
glass_X_train <- glass_X_norm[folds[[1]], ]
glass_y_test <- glass_y[-folds[[1]]]
glass_X_test <- glass_X_norm[-folds[[1]], ]

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
#'@param K Number of class labels - 1 / Dimension of latent vector
#'@param seed Seed used for controlling randomness
#'
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
                        sigma_hat = NULL, #????
                        shape = 1, # shape parameter of gating probabilities????
                        width = 0.1, # bandwidth of gating probabilities
                        alpha_scale = NULL, # scale of hyperprior on alpha
                        alpha_shape_1 = 0.5, # shape 1 of hyperprior on alpha
                        alpha_shape_2 = 1, # shape 2 of hyperprior on alpha
                        tau_rate = 10, # rate of exponential prior on bandwidth parameters tau
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
  
  z <- matrix(NA_real_, nrow = num_obs_train, ncol = K) 
  mu_z <- rep(0, K) #?
  Sigma <- diag(K) # set Sigma to identity matrix
  
  nu_prior <- K + 1 # prior d.o.f inv-Wishart
  scalematr_prior <- nu_prior * diag(K)
  
  # create parameters and tree samplers for softBART, and # initialize lists and matrices to store draws, predictions and errors
  hypers <- vector("list", K)
  tree_samplers <- vector("list", K)
  predictions_z <- matrix(NA_real_, nrow = num_obs_train, ncol = K)
  
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
  
  z_draws <- vector("list", num_sim)
  test_z_preds <- matrix(NA_real_, num_obs_test, K) # maybe remove, or maybe not
  errors <- matrix(NA_real_, num_obs_train, K)
  
  opts <- Opts(update_sigma = FALSE, num_print = n.burnin + n.iter + 1) # check settings
  
  # ------- MCMC --------
  
  # Gibbs sampler
  for (iter in 1:num_burnin+num_sim) {
    
    # sample latent variables from truncated multivariate normal
    z <- t(sapply(1:num_obs_train, function(i) {
      sample_latent_variables(mu = mu_z, Sigma = Sigma, y_i = y_train[i], K = K)
    }))
    
    # check if correct dimension
    if (iter == 1){
      print(dim(z))
    }
    
    # sample all tree model related parameters using softBART package
    for (k in 1:K) {
      
      # compute input z for tree sampler (requires initial predictions?)
      temp_z <- z[,k] - t((t(Sigma[-k,k]) %*% solve(Sigma[-k,-k])) %*% t(z[,-k] - predictions_z_train[,-k]))
      
      #sampler.list[[k]]$set_sigma(sqrt(Sigma_mat[k,k] - Sigma_mat[k,-k]%*%(solve(Sigma_mat[-k,-k]) %*% Sigma_mat[-k,k] ))  ) #????? maybe remove
       
      # predict k'th component of z (training data)
      mu_train <- t(tree_samplers[[k]]$do_gibbs(X_train, temp_z, X_train, i = 1)) # returns predictions for all training obs of k'th component
      
      # compute errors 
      errors[,k] <- z[,k] - mu_train
      
      # update k'th component of mu_z????
      mu_z[k] <- mean(predictions_z_k)
      
      predictions_z_train[,k] <- mu_train # save predictions
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
      z_draws[[iter - num_burnin]] <- z
    }
    
  } # end of sampler
  
  # return MCMC output
  return(z_draws)
}



  





