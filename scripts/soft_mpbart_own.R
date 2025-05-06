# SOFT MPBART


# packages
library(MASS)
library(MCMCpack)

# ------------ READ DATA ---------------


# ------------ PREPROCESS DATA -------------




soft_mpbart <- function(y_train, # training data - outcomes
                        X_train, # training data - covariates
                        y_test, # test data - outcomes
                        X_test, # test data - covariates
                        num_burnin, # number of burn-in iterations
                        num_sim, # number of simulations (excl. burn-in)
                        num_trees, # number of trees in the sum-of-trees model
                        K, # number of outcome categories - 1 (dim of latent vector)
                        ) {
  
  # set some values
  num_obs_train = length(y_train)
  num_obs_test = length(y_test)
  
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
    for (i in 1:num_obs_train) {
      
    }
    
    
    
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




