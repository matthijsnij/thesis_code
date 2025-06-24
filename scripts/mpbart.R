# MPBART with hard splits (Xu et al (2024))

library(GcompBART)

# ------ MPBART FUNCTION -------

mpbart <- function(X_train, # matrix of covariates - training data
                   y_train, # vector of class labels - training data
                   X_test, # matrix of covariates - test data
                   y_test, # vector of class labels - test data
                   num_burnin = 1000,
                   num_sim = 1000,
                   num_trees = 200,
                   num_classes) {
   
  
  p <- num_classes - 1 # dimension latent
  
  # prior settings
  Prior = list(nu=p+2, # Sigma prior; d.o.f.
               V= diag(p), # Sigma prior; scale matrix
               ntrees=200, # number of trees
               kfac=2.0, # leaf node param prior constant
               pswap=0,  # prob of swap move in tree prior; set to 0
               pbd=1.0, # prob of birth/death; set to one
               pb=0.5 , # prob of birth given birth/death
               beta = 2.0, # branching process prior
               alpha = 0.95, # branching process prior
               nc = 100, # number of equally spaced cutpoints between min and max for each covariate
               minobsnode = 10) # minimum number of observations in bottom nodes for birth in simulating trees
  
  # mcmc settings
  Mcmc <- list(burn = num_burnin, ndraws = num_sim)
  
  df_train <- data.frame(y = y_train, X_train)
  reference_class <- 0 # reference level
  
  model_fit <- model_bart(formula = y ~ ., 
                      data = df_train, 
                      type = "multinomial",
                      base = reference_class,
                      Prior = Prior,
                      Mcmc = Mcmc,
                      correction = FALSE, 
                      )
  
  predict_obj <- predict_bart(obj = model_fit, newdata = as.data.frame(X_test))
  
}

