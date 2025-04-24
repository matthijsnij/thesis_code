# soft mpbart
# original code: https://github.com/EoghanONeill/SURBART/blob/master/R/soft_surbart_original.R
# which I edited for multinomial probit framework, i.e. sampling latent variables in MCMC, and taking categorical input


make_01_norm <- function(x) {
  a <- min(x)
  b <- max(x)
  return(function(y0) (y0 - a) / (b - a))
}

# --------------- DESCRIPTION

#' @title Soft Multinomial Probit Bayesian Additive Regression Trees implemented using MCMC by changing the SURBART approach from Chakraborty (2016).
#'
#' @description Soft Multinomial Probit Bayesian Additive Regression Trees implemented using MCMC by changing the SURBART approach from Chakraborty (2016).
#' In particular, code is changed to accommodate categorical input + sample the latent variables in the MCMC.
#' @import dbarts
#' @import truncnorm
#' @import LaplacesDemon
#' @import MASS
#' @import collapse
#' @import SoftBart
#' what follows now are descriptions of soft surbart by eoghan, so do not reflect MNP framework
#' @param x.train The training covariate data for all training observations. Matrix or list. If one matrix (not a list), then the same set of covariates are used for each outcome. Each element of the list is a covariate matrix that  corresponds to a different outcome variable. Number of rows equal to the number of observations. Number of columns equal to the number of covariates.
#' @param x.test The test covariate data for all test observations. Matrix or list. If one matrix (not a list), then the same set of covariates are used for each outcome. Each element of the list is a covariate matrix that  corresponds to a different outcome variable. Number of rows equal to the number of observations. Number of columns equal to the number of covariates.
#' @param y The training data list of vectors of outcomes
#' @param num_outcomes The number of outcome variables.
#' @param num_obs The number of observations per outcome.
#' @param num_test_obs THe number of test observations per outcome.
#' @param n.iter Number of iterations excluding burnin.
#' @param n.burnin Number of burnin iterations.
#' @param n.trees (dbarts control option) A positive integer giving the number of trees used in the sum-of-trees formulation. Assuming each outcome variable is modelled using the same number of trees. Each outcome is modelled using a distinct sum of trees (as opposed to shared trees).
#' @param n.chains (dbarts control option) A positive integer detailing the number of independent chains for the dbarts sampler to use (more than one chain is unlikely to improve speed because only one sample for each call to dbarts).
#' @param n.threads  (dbarts control option) A positive integer controlling how many threads will be used for various internal calculations, as well as the number of chains. Internal calculations are highly optimized so that single-threaded performance tends to be superior unless the number of observations is very large (>10k), so that it is often not necessary to have the number of threads exceed the number of chains.
#' @param printEvery (dbarts control option)If verbose is TRUE, every printEvery potential samples (after thinning) will issue a verbal statement. Must be a positive integer.
#' @param printCutoffs (dbarts control option) A non-negative integer specifying how many of the decision rules for a variable are printed in verbose mode
#' @param rngKind (dbarts control option) Random number generator kind, as used in set.seed. For type "default", the built-in generator will be used if possible. Otherwise, will attempt to match the built-in generator’s type. Success depends on the number of threads.
#' @param rngNormalKind (dbarts control option) Random number generator normal kind, as used in set.seed. For type "default", the built-in generator will be used if possible. Otherwise, will attempt to match the built-in generator’s type. Success depends on the number of threads and the rngKind
#' @param rngSeed (dbarts control option) Random number generator seed, as used in set.seed. If the sampler is running single-threaded or has one chain, the behavior will be as any other sequential algorithm. If the sampler is multithreaded, the seed will be used to create an additional pRNG object, which in turn will be used sequentially seed the threadspecific pRNGs. If equal to NA, the clock will be used to seed pRNGs when applicable.
#' @param updateState (dbarts control option) Logical setting the default behavior for many sampler methods with regards to the immediate updating of the cached state of the object. A current, cached state is only useful when saving/loading the sampler.
#' @param tree.prior (dbarts option) An expression of the form dbarts:::cgm or dbarts:::cgm(power,base) setting the tree prior used in fitting.
#' @param node.prior (dbarts option) An expression of the form dbarts:::normal or dbarts:::normal(k) that sets the prior used on the averages within nodes.
#' @param resid.prior (dbarts option) An expression of the form dbarts:::chisq or dbarts:::chisq(df,quant) that sets the prior used on the residual/error variance
#' @param proposal.probs (dbarts option) Named numeric vector or NULL, optionally specifying the proposal rules and their probabilities. Elements should be "birth_death", "change", and "swap" to control tree change proposals, and "birth" to give the relative frequency of birth/death in the "birth_death" step.
#' @param sigmadbarts (dbarts option) A positive numeric estimate of the residual standard deviation. If NA, a linear model is used with all of the predictors to obtain one.
#' @param print.opt Print every print.opt number of Gibbs samples.
#' @param quiet Does not show progress bar if TRUE
#' @param outcome_draws If TRUE, output draws of the outcome.
#' @export
#' @return The following objects are returned:
#' \item{mutrain_draws}{Matrix of MCMC draws of expected values for training observations. Number of rows equal to the number of training observations multiplied by the number of outcomes. The rows are ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{mutest_draws}{Matrix of MCMC draws of expected values for test observations. Number of rows equal to the number of test observations multiplied by the number of outcomes. The rows are ordered by beginning with all Ntest (units') observations for the first outcome variable, then all Ntest for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{ytrain_draws}{Matrix of MCMC draws of outcomes for training observations. Number of rows equal to the number of observations multiplied by the number of outcomes. The rows should be ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{ytest_draws}{Matrix of MCMC draws of outcomes for training observations. Number of rows equal to the number of observations multiplied by the number of outcomes. The rows should be ordered by beginning with all N (units') observations for the first outcome variable, then all N for the second outcome variable and so on. Number of columns equals n.iter.}
#' \item{Sigma_draws}{3 dimensional array of MCMC draws of the covariance matrix for the outcome-specific error terms. The numbers of rows and columns equal are equal to the number of outcome variables. The number of slices is }
#' @export

# ------------- SOFT MPBART
soft_mpbart_original <- function(x.train, # predictor matrix training set
                                  x.test, # predictor matrix test set
                                  y.train, # categorical response training set
                                  y.test, # categorical response test set
                                  K, # number of response variables categories
                                  num_latent, # number of latent variables (was num_outcomes)
                                  num_obs, # number of train obs 
                                  num_test_obs, # number of test obs 
                                  n.iter=1000, # number of iterations excluding burnin
                                  n.burnin=100, # number of burnin iterations
                                  n.trees = 50L, # (dbarts control) number of trees in sum-of-trees. Each latent variable will have a distinct sum-of-trees, but the number of trees in them is the same.
                                  n.burn = 0L, # not used? (oc)
                                  n.samples = 1L, # not used? (oc)
                                  n.thin = 1L, # not used? (oc)
                                  n.chains = 1, # (dbarts control) number of indep chains for dbarts sampler to use (oc)
                                  n.threads = guessNumCores(), # (dbarts control) number of CPU cores for computation, not sure if that relevant
                                  printEvery = 100L, # (dbarts control) not used? (oc). Progress updates every 100 iterations
                                  printCutoffs = 0L, # (dbarts control) not used? (oc). Print cutoff values in splitting. Probably can go due to soft splits?
                                  rngKind = "default", # (dbarts control) type of RNG. not used? (oc)
                                  rngNormalKind = "default", # (dbarts control) how to generate std. normals. not used? (oc)
                                  rngSeed = NA_integer_, # (dbarts control) set.seed for dbarts algorithms. not used? (oc)
                                  updateState = FALSE, # (dbarts control) TRUE = full posterior inference but slower, FALSE = only for prediction, faster. TRUE is default. not used? (oc)
                                  tree.prior = dbarts:::cgm, # (dbarts option) the tree prior. Default cgm(0.95, 2). Pass as dbarts:::cgm(0.95, 2). not used? (oc)
                                  node.prior = dbarts:::normal, # (dbarts option) the terminal node param prior. Pass as dbarts:::normal(k). As data centered around 0. Here k controls regularization in the prior variance sigma^2 / k.
                                  resid.prior = dbarts:::chisq, # (dbarts option) residual variance prior. Pass as dbarts:::chisq(3, 0.95).
                                  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5), # (dbarts option) MH proposal probabilities. Might want to set swap to 0.
                                  sigmadbarts = NA_real_, # (dbarts option) residual stddev estimate. If NA, linear model with all predictors is used.
                                  print.opt = 100, # (dbarts control) Controls output formatting and verbosity.
                                  quiet = FALSE, # does not show progress bar if TRUE.
                                  outcome_draws = FALSE, # if TRUE, output draws of the outcome.
                                  SB_group = NULL,
                                  SB_alpha = 1,
                                  SB_beta = 2,
                                  SB_gamma = 0.95,
                                  SB_k = 2,
                                  SB_sigma_hat = NULL,
                                  SB_shape = 1,
                                  SB_width = 0.1,
                                  # SB_num_tree = 20,
                                  SB_alpha_scale = NULL,
                                  SB_alpha_shape_1 = 0.5,
                                  SB_alpha_shape_2 = 1,
                                  SB_tau_rate = 10,
                                  SB_num_tree_prob = NULL,
                                  SB_temperature = 1,
                                  SB_weights = NULL,
                                  SB_normalize_Y = TRUE){
  
  # -------------------- INPUT DATA CHECKS -------------------------
  
  # check if covariates are in a matrix
  if(!is.matrix(x.train)) stop("argument x.train must be a double matrix.")
  if(!is.matrix(x.test)) stop("argument x.test must be a double matrix.")
  
  # check if response is in a vector
  if(!is.vector(y.train)) stop("argument y.train must be a vector.")
  if(!is.vector(y.test)) stop("argument y.test must be a vector.")
  
  # check if correct number of observations
  if(nrow(x.train) != num_obs) stop("number of rows of x.train must equal num_obs.")
  if(nrow(x.test) != num_test_obs) stop("number of rows of x.test must equal num_test_obs.")
  if(length(y.train) != num_obs) stop("length of y.train must equal number of observations.")
  if(length(y.train) != num_test_obs) stop("length of y.train must equal number of observations.")
  
  # check if number of columns for x.test and x.train are the same
  if((ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train.")
  
  # check if all categories exist in y
  if(length(unique(y.train))!= K) stop("number of distinct values y.train does not equal number of categories K.")
  if(length(unique(y.test))!= K) stop("number of distinct values y.test does not equal number of categories K.")
  
  
  for(matind in 1:num_outcomes){
    
    x.traintemp <- Xlist[[matind]]
    x.testtemp <- Xtestlist[[matind]]
    
    ecdfs   <- list()
    for(i in 1:ncol(x.traintemp)) {
      ecdfs[[i]] <- ecdf(x.traintemp[,i])
      if(length(unique(x.traintemp[,i])) == 1) ecdfs[[i]] <- identity
      if(length(unique(x.traintemp[,i])) == 2) ecdfs[[i]] <- make_01_norm(x.traintemp[,i])
    }
    for(i in 1:ncol(x.traintemp)) {
      x.traintemp[,i] <- ecdfs[[i]](x.traintemp[,i])
      if(nrow(x.testtemp) > 0){
        x.testtemp[,i] <- ecdfs[[i]](x.testtemp[,i])
      }
    }
    
    rm(ecdfs)
    
    Xlist[[matind]] <- x.traintemp
    Xtestlist[[matind]] <- x.testtemp
    
  }
  
  
  n <- length(y)
  
  # draw = list(
  #   Z.mat = array(NA, dim = c(n, n.iter)),
  #   Z.matcens = array(NA, dim = c(n0, n.iter)),
  #   #Z.matuncens = array(NA, dim = c(n1, n.iter)),
  #   Z.matcensbelow = array(NA, dim = c(n_censbelow, n.iter)),
  #   Z.matcensabove = array(NA, dim = c(n_censabove, n.iter)),
  #   mu = array(NA, dim = c(n, n.iter)),#,
  #   mucens = array(NA, dim = c(n0, n.iter)),#,
  #   muuncens = array(NA, dim = c(n1, n.iter)),#,
  #   mucensbelow = array(NA, dim = c(n_censbelow, n.iter)),#,
  #   mucensabove = array(NA, dim = c(n_censabove, n.iter)),#,
  #   ystar = array(NA, dim = c(n, n.iter)),#,
  #   ystarcens = array(NA, dim = c(n0, n.iter)),#,
  #   ystaruncens = array(NA, dim = c(n1, n.iter)),#,
  #   ystarcensbelow = array(NA, dim = c(n_censbelow, n.iter)),#,
  #   ystarcensabove = array(NA, dim = c(n_censabove, n.iter)),#,
  #   test.mu =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.y_nocensoring =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.y_withcensoring =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.probcensbelow =  array(NA, dim = c(ntest, n.iter)),#,
  #   test.probcensabove =  array(NA, dim = c(ntest, n.iter)),
  #   sigma = rep(NA, n.iter)
  # )
  
  
  
  Ymu <- rep(NA, num_outcomes)
  Ysd <- rep(NA, num_outcomes)
  
  
  # standardize outcomes
  for(i in 1:num_outcomes){
    Ymu[i] <- mean(y[[i]])
    Ysd[i] <- sd(y[[i]])
    y[[i]] <- (y[[i]]-Ymu[i])/Ysd[i]
  }
  
  ymat <- matrix(NA, nrow = num_obs, ncol = num_outcomes)
  
  for(i in 1:num_outcomes){
    ymat[,i] <- y[[i]]
  }
  
  # set various prior parameter values
  
  rprior <- num_outcomes +1
  Rprior <- rprior*diag(num_outcomes)
  
  # initialize sigma draw
  
  rss_initial <- t(ymat) %*% ymat
  
  print("rss_initial = ")
  print(rss_initial)
  
  Sigma_mat <- rinvwishart(rprior + num_obs, Rprior + rss_initial)
  
  # ch <- chol(Sigma_mat)
  # dd <- diag(ch)
  # Lmat <- t(ch/dd)
  # Hvec <- dd
  
  #Set tree sampler parameters
  
  # control <- dbartsControl(updateState = updateState, verbose = FALSE,  keepTrainingFits = TRUE,
  #                          keepTrees = TRUE,
  #                          n.trees = n.trees,
  #                          n.burn = n.burn,
  #                          n.samples = n.samples,
  #                          n.thin = n.thin,
  #                          n.chains = n.chains,
  #                          n.threads = n.threads,
  #                          printEvery = printEvery,
  #                          printCutoffs = printCutoffs,
  #                          rngKind = rngKind,
  #                          rngNormalKind = rngNormalKind,
  #                          rngSeed = rngSeed)
  
  
  opts <- Opts(update_sigma = FALSE, num_print = n.burnin + n.iter + 1)
  
  
  hyperslist <- list()
  
  
  for(i in 1:num_outcomes){
    
    hyperslist[[i]] <- Hypers(Xlist[[i]], y[[i]],
                              num_tree = n.trees, #sigma_hat = 1,
                              group = SB_group,
                              alpha = SB_alpha,
                              beta = SB_beta,
                              gamma = SB_gamma,
                              k = SB_k,
                              sigma_hat = NULL, #sighat,
                              shape = SB_shape,
                              width = SB_width,
                              # num_tree = 20,
                              alpha_scale = SB_alpha_scale,
                              alpha_shape_1 = SB_alpha_shape_1,
                              alpha_shape_2 = SB_alpha_shape_2,
                              tau_rate = SB_tau_rate,
                              num_tree_prob = SB_num_tree_prob,
                              temperature = SB_temperature,
                              weights = SB_weights,
                              normalize_Y = SB_normalize_Y
    )
    
  }
  
  
  
  
  # print(colnames(Xmat.train))
  
  # print("begin dbarts")
  
  
  
  
  #take initial tree draws
  
  sampler.list <- list()
  
  # if(nrow(Xtestlist[[1]])==0){
  preds.train <- matrix(0, num_obs, num_outcomes)
  
  
  for (jj in 1:num_outcomes){
    
    sampler.list[[jj]] <- MakeForest(hyperslist[[jj]], opts, warn = FALSE)
    
    preds.train[,jj] <- sampler.list[[jj]]$do_predict( Xlist[[jj]])
    
    # Xmat.train <- data.frame(y = y[[jj]], x = Xlist[[jj]] )
    #
    # sampler <- dbarts(y ~ .,
    #                   data = Xmat.train,
    #                   #test = x.test,
    #                   control = control,
    #                   tree.prior = tree.prior,
    #                   node.prior = node.prior,
    #                   resid.prior = resid.prior,
    #                   proposal.probs = proposal.probs,
    #                   sigma = sigmadbarts)
    #
    # sampler.list[[jj]] <- sampler
    
  }
  
  
  # }else{
  #   for (jj in 1:num_outcomes){
  #
  #     Xmat.train <- data.frame(y = y[[jj]], x = Xlist[[jj]] )
  #     Xmat.test <- data.frame(x = Xtestlist[[jj]] )
  #
  #     sampler <- dbarts(y ~ .,
  #                       data = Xmat.train,
  #                       test = Xmat.test,
  #                       control = control,
  #                       tree.prior = tree.prior,
  #                       node.prior = node.prior,
  #                       resid.prior = resid.prior,
  #                       proposal.probs = proposal.probs,
  #                       sigma = sigmadbarts
  #     )
  #
  #     sampler.list[[jj]] <- sampler
  #
  #   }
  #
  # }
  
  sampler.run <- list()
  
  
  
  preds.test <- matrix(NA, num_test_obs, num_outcomes)
  
  # sigma.mat <- matrix(NA, num_outcomes, num_outcomes)
  
  # count.mat <- list() #stored as list to allow for varying numbers of splitting variables per outcome
  # require separate counts for each tree
  
  
  Y_store <- array(NA,dim=c(n.iter,num_obs,num_outcomes))
  mu_store <- array(NA,dim=c(n.iter,num_obs,num_outcomes))
  
  Ytest_store <- array(NA,dim=c(n.iter,num_test_obs,num_outcomes))
  mutest_store <- array(NA,dim=c(n.iter,num_test_obs,num_outcomes))
  
  Sigma_store <- array(NA,dim=c(num_outcomes,num_outcomes, n.iter))
  
  
  
  eta <- matrix(NA, nrow = num_obs, ncol =  num_outcomes)
  
  # start Gibbs sampler
  # show progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = n.iter+n.burnin, style = 3)
    start  <- Sys.time()
  }
  
  
  
  #loop through the Gibbs sampler iterations
  for(iter in 1:(n.iter+n.burnin)){
    
    
    for (mm in 1:num_outcomes){
      # if (mm > 1){
      # Z_mm <- eta[,1:(mm-1), drop=F]
      # A0_mm <- A0_draw[mm,1:(mm-1)]
      
      #LDL decomposition of covariance matrix
      #D in LDL deomposition, variances of each error in transformed equations
      
      # A0_mm <- Lmat[mm,1:(mm-1)]
      # print("Z_mm =")
      # print(Z_mm)
      #
      # print("A0_mm =")
      # print(A0_mm)
      
      # if(iter ==1){
      #   # sampler.list[[mm]]$setResponse(y = ymat[,mm] )
      #   tempy <- ymat[,mm] #- (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  )
      # }else{
      tempy <- ymat[,mm] - (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  )
      
      # sampler.list[[mm]]$setResponse(y = ymat[,mm] - (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  ) )
      sampler.list[[mm]]$set_sigma(sqrt(Sigma_mat[mm,mm] - Sigma_mat[mm,-mm]%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm] ))  )
      # }
      
      # }else{
      #   tempy <- ymat[,mm] - (ymat[,-mm]-preds.train[,-mm])%*%(solve(Sigma_mat[-mm,-mm]) %*% Sigma_mat[-mm,mm]  )
      # }
      
      
      
      #Draw all trees for outcome mm in iteration iter
      # rep_mm <- sampler.list[[mm]]$run(0L, 1L) # construct BART sample using dbarts (V. Dorie)
      mutrain <- t(sampler.list[[mm]]$do_gibbs(Xlist[[mm]], tempy, Xlist[[mm]], 1))
      # mutest <- sampler.list[[jj]]$do_predict(x.test)
      
      # sampler.run[[mm]] <- rep_mm
      
      #line commented out because sampling from inverse wishart, not the prior from BAVART paper
      # sigma.mat[mm,] <- rep_mm$sigma
      
      if (any(is.na(mutrain))){
        stop("Line 468. NA in predictions")
      }
      eta[,mm] <- ymat[,mm] - mutrain
      # A_draw[,mm] <- X.ginv%*%rep_mm$train
      preds.train[,mm] <-  mutrain
      
      
      if(nrow(Xtestlist[[1]])>0){
        mutest <- sampler.list[[mm]]$do_predict(Xtestlist[[mm]])
        preds.test[,mm] <- mutest
      }
      
      
      # print("rep_mm$test[,1] =")
      # print(rep_mm$test[,1])
      # preds.test[,mm] <-  rep_mm$test[,1]
      
      # count.mat[,mm] <- rep_mm$varcount
      # if (mm > 1){
      #   norm_mm <- as.numeric(exp(-.5*sv_latent[[mm]]) * 1/sigma.mat[mm,])
      #   u_mm <- eta[,1:(mm-1),drop=F]*norm_mm
      #   eta_mm <- eta[,mm]*norm_mm
      #   if (mm == 2) v0.inv <- 1/theta_A0[mm,1] else v0.inv <- diag(1/theta_A0[mm,1:(mm-1)])
      #   V.cov <- solve(crossprod(u_mm) + v0.inv)
      #   mu.cov <- V.cov %*% crossprod(u_mm, eta_mm)
      #   mu.cov.draw <- mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov))
      #   A0_draw[mm,1:(mm-1)] <- mu.cov.draw
      # }
      
      
    } # end loop over mm
    
    
    
    # draw next sigma matrix
    rss <- t(eta) %*% eta
    
    Sigma_mat <- rinvwishart(nu = rprior + num_obs, S = Rprior + rss)
    
    # ch <- chol(Sigma_mat)
    # dd <- diag(ch)
    # Lmat <- t(ch/dd)
    # Hvec <- dd
    
    
    #now save training and test draws if past burn-in
    
    
    if(iter>n.burnin){
      iter_min_burnin <- iter-n.burnin
      
      #save mu draws
      for(mm in 1:num_outcomes){
        mu_store[iter_min_burnin,1:num_obs,mm] <- preds.train[,mm]*Ysd[mm] + Ymu[mm]
        mutest_store[iter_min_burnin,1:num_test_obs,mm] <- preds.test[,mm]*Ysd[mm] + Ymu[mm]
      }
      #save y draws
      
      
      if(outcome_draws==TRUE){
        #draw each individual's vector of observations from a multivariate normal
        for(obs_ind in 1:num_obs){
          temp_sample <- mvrnorm(n = 1,
                                 mu = preds.train[obs_ind,],
                                 Sigma = Sigma_mat)
          
          #Note: element-wise vector multiplication
          Y_store[iter_min_burnin,obs_ind , 1:num_outcomes] <- temp_sample*Ysd + Ymu
          
        }
        
        for(obs_ind in 1:num_test_obs){
          temp_sample <- mvrnorm(n = 1,
                                 mu = preds.test[obs_ind,],
                                 Sigma = Sigma_mat)
          
          #Note: element-wise vector  multiplication
          Ytest_store[iter_min_burnin,obs_ind , 1:num_outcomes] <- temp_sample*Ysd + Ymu
          
        }
      }
      
      #save sigma matrix draws
      #if mm==1, then diagonal with sigma given by independent value from dbarts initial value? or save next draw?
      
      
      
      Sigma_store[,,iter_min_burnin] <- Sigma_mat
      
      
    } # end if n.iter > n.burnin
    
    if(!quiet) setTxtProgressBar(pb, iter)
    
    # if(iter %% print.opt == 0){
    #   print(paste("Gibbs Iteration", iter))
    #   # print(c(sigma2.alpha, sigma2.beta))
    # }
    
    
  }#end iterations of Giibs sampler
  
  
  ret_list <- list()
  
  ret_list$mu_draws <- mu_store
  ret_list$mutest_draws <- mutest_store
  
  if(outcome_draws==TRUE){
    ret_list$Ytrain_draws <- Y_store
    ret_list$Ytest_draws <- Ytest_store
  }
  
  ret_list$Sigma_draws <- Sigma_store
  
  return(ret_list)
  
  
  
}