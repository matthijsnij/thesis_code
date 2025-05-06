# SOFT MPBART


# packages
library(caret)
library(MASS)
library(MCMCpack)
library(TruncatedNormal)

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


# ------------ PREPROCESS DATA -------------

# clean the class labels such that they fall in the range [0,5]
# there is no class 4 in the data set
for (i in 1:length(glass_y)) {
  if (glass_y[i] < 4) {
    glass_y[i] <- glass_y[i] - 1
  } else {
    glass_y[i] <- glass_y[i] - 2
  }
}


# ------------ CREATE TRAIN AND TEST SETS

# generate train and test set folds
num_folds <- 5
folds <- createFolds(glass_y, k = num_folds, list = TRUE, returnTrain = TRUE)

# split the data based on these indices
glass_y_train <- glass_y[folds[[1]]]
glass_X_train <- glass_X[folds[[1]], ]
glass_y_test <- glass_y[-folds[[1]]]
glass_X_test <- glass_X[-folds[[1]], ]

# ------------- SOFT MPBART FUNCTION -----------------

#'@description MCMC algorithm for soft MPBART
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
                        num_burnin, # number of burn-in iterations
                        num_sim, # number of simulations (excl. burn-in)
                        num_trees, # number of trees in the sum-of-trees model
                        K, # number of outcome categories - 1 (dim of latent vector)
                        seed, # seed used in set.seed to control randomness. If not passed, function will generate random seed
                        ) {
  
  # set some values
  num_obs_train = length(y_train)
  num_obs_test = length(y_test)
  
  # set seed for reproduceability
  if (is.null(seed)) {
    seed <- as.integer(Sys.time()) %% .Machine$integer.max # pseudo-random seed if not provided
  }
  set.seed(seed = seed)
  
  # ------------ INITIALIZE PARAMS ----------------
  
  z <- matrix(NA_real_, nrow = num_obs_train, ncol = K) 
  mu_z <- rep(0, K)
  Sigma <- diag(K) # set Sigma to identity matrix
  
  z <- mvrnorm(num_obs_train, mu = mu_z, Sigma = Sigma) # initialize latent variables from standard mv normal # remove?
  
  nu_prior <- K + 1 # prior d.o.f inv-Wishart
  scalematr_prior <- nu_prior * diag(K)
  
  # initialize each sum-of-trees model with 'num_trees' single node trees + what to do with leaf node params
  
  
  

  # ------- MCMC --------
  
  # initialize lists to store draws
  z_draws <- vector("list", num_sim)
  
  
  for (iter in 1:num_burnin+num_sim) {
    
    # sample latent variables from truncated multivariate normal
    z <- t(sapply(1:num_obs_train, function(i) {
      sample_latent_variables(mu = mu_z, Sigma = Sigma, y_i = y_train[i], K = K)
    }))
    
    for (k in 1:K) {
      
       # sample all trees, leaf node params and tau via softBART package
       for (j in 1:num_trees) {
         
         
       }
      
      
       # sample splitting probabilities from Dirichlet, should also come out of softBART
      
      
    }
    
    # sample unconstrained Sigma from inverted-Wishart
    nu_posterior <- nu_prior + num_obs_train
    scalematr_posterior <- scalematr_prior + #RSS treemodel
    Sigma_star <- riwish(nu_posterior, scalematr_posterior)
    
    # scale to force trace restriction
    Sigma <- (K / sum(diag(Sigma_star))) * Sigma_star
    
    
    # save draws if after burn-in
    if (iter > num_burnin) {
      z_draws[[iter - num_burnin]] <- z
    }
    
  }
  
  
  
}


  





