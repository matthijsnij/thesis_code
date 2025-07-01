# MPBART with hard splits (Xu et al (2024))

library(GcompBART)

# ------ MPBART FUNCTION -------

mpbart <- function(X_train, # matrix of covariates - training data
                   y_train, # vector of class labels - training data
                   X_test, # matrix of covariates - test data
                   y_test, # vector of class labels - test data
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
  
  df_train <- data.frame(y = y_train, X_train)
  reference_class <- 0 # reference level
  
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
  samp_y <- predict_obj$samp_y  # ndraws x n matrix
  ndraws <- nrow(samp_y)
  n_test <- ncol(samp_y)
  classes <- sort(unique(y_train))  # class labels, e.g. 0, 1, 2
  
  # Compute posterior class probabilities: n_test x num_classes
  class_probs <- sapply(classes, function(k) {
    colMeans(samp_y == k)
  })
  
  # Ensure class_probs is a matrix with rows = test samples, columns = classes
  class_probs <- as.matrix(class_probs)
  colnames(class_probs) <- paste0("class_", classes)
  
  # Predicted class = class with highest posterior probability
  y_pred <- apply(class_probs, 1, function(row) classes[which.max(row)])
  
  # Compute metrics
  err <- test_error_rate(y_test, y_pred)
  brier <- brier_score_multiclass(y_test, class_probs)
  
  # Return everything
  return(list(
    model_fit = model_fit,
    predicted_class = y_pred,
    predicted_prob = class_probs,
    test_error = err,
    brier_score = brier
  ))
  
  
}

