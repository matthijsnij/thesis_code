# MPBART with hard splits (Xu et al (2024))

library(GcompBART)

# ------ MPBART FUNCTION -------

#'@description Function for running MPBART as implemented by Xu et al. (2024)
#'
#'@param X_train Matrix of covariates - train data
#'@param y_train Vector of class labels - train data
#'@param X_test Matrix of covariates - test data
#'@param num_burnin Number of burn-in iterations
#'@param num_sim Number of iterations, post burn-in
#'@param num_trees Number of trees in the BART model
#'@param num_classes Number of outcome classes
#'@return List containing
#'\item{pred_y}{Vector of predicted class labels on test set}
#'\item{post_probs}{Matrix of predicted class probabilities on test set}
mpbart <- function(X_train, # matrix of covariates - training data
                   y_train, # vector of class labels - training data
                   X_test, # matrix of covariates - test data
                   num_burnin = 1500,
                   num_sim = 1500,
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
  
  df_train <- data.frame(y = factor(y_train, levels = 0:(num_classes - 1)), X_train)
  reference_class <- "0" # reference level
  
  # fit MPBART
  model_fit <- model_bart(formula = y ~ ., 
                      data = df_train, 
                      type = "multinomial",
                      base = reference_class,
                      Prior = Prior,
                      Mcmc = Mcmc,
                      correction = FALSE, 
                      )
  
  # predict on test data
  predict_obj <- predict_bart(obj = model_fit, newdata = as.data.frame(X_test))
  
  # extract results
  samp_y <- t(predict_obj$samp_y)  # ndraws x n matrix of predicted class labels
  classes <- 0:(num_classes - 1)
  
  # compute posterior class probabilities: n_test x num_classes
  class_probs <- sapply(classes, function(k) {
    colMeans(samp_y == k)
  })
  class_probs <- as.matrix(class_probs)
  colnames(class_probs) <- paste0("class_", classes)
  
  # predicted class = class with highest posterior probability
  pred_y <- apply(class_probs, 1, function(row) classes[which.max(row)])
  
  return_list = list(
    pred_y = pred_y,
    post_probs = class_probs
  )
  
  return(return_list)
}

