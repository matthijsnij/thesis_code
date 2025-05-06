# SOFT MPBART


# packages
library(MASS)
library(MCMCpack)
library(TruncatedNormal)

# ---------------- FUNCTION TO SAMPLE FROM TRUNCATED MULTIVARIATE NORMAL ---------------------

# the function assumes category 0 is the reference level
sample_truncated_mvnormal <- function(mu, # mean vector Kx1
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


# ------------ PREPROCESS DATA -------------



# ------------- SOFT MPBART FUNCTION -----------------
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
    seed <- as.integer(Sys.time()) %% .Machine$integer.max
  }
  set.seed(seed = seed)
  
  # ------------ INITIALIZE PARAMS ----------------
  
  z <- matrix(NA_real_, nrow = num_obs_train, ncol = K) 
  mu_z <- rep(0, K)
  Sigma <- diag(K) # set Sigma to identity matrix
  
  z <- mvrnorm(num_obs_train, mu = mu_z, Sigma = Sigma) # initialize latent variables from standard mv normal
  
  nu_prior <- K + 1 # prior d.o.f inv-Wishart
  scalematr_prior <- nu_prior * diag(K)
  
  # initialize each sum-of-trees model with 'num_trees' single node trees + what to do with leaf node params
  
  
  

  # ------- MCMC --------
  for (iter in 1:num_burnin+num_sim) {
    
    # sample latent variables from truncated multivariate normal
    z <- t(sapply(1:num_obs_train, function(i) {
      sample_truncated_mvnorm(mu = mu_z, Sigma = Sigma, y_i = y_train[i], K = K)
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
      
    }
    
  }
  
  
  
}


  





