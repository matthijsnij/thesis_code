# Smooth DGP simulation for S-MPBART - Friedman (1991) inspired
# Benchmark methods are Random Forest and MPBART with hard splits as implemented by Xu et al. (2024)

# ---------- FRIEDMAN SMOOTH NONLINEAR DGP ------------

# set seed for reproduceability
set.seed(123)

# set replications
reps <- 10
friedman_error_rates <- rep(0, reps)
friedman_brier_scores <- rep(0, reps)

# create workbooks
wb <- createWorkbook() # for data
wb_output <- createWorkbook() # for output

for (r in 1:reps) {
  # generate predictors from U(0,1)
  n <- 1000      # number of observations
  p <- 10        # number of predictors
  friedman_X <- matrix(runif(n * p), nrow = n, ncol = p) # no rank normalization necessary
  colnames(friedman_X) <- paste0("X", 1:p)
  
  # latent functions for each class (smooth and nonlinear)
  f1 <- 10 * sin(pi * friedman_X[,1] * friedman_X[,2]) + 20 * (friedman_X[,3] - 0.5)^2
  f2 <- 10 * friedman_X[,4] + 5 * friedman_X[,5] + 2 * friedman_X[,6]
  f3 <- 15 * friedman_X[,7] * friedman_X[,8] - 10 * (friedman_X[,9] - 0.5)^2
  
  # combine and add noise
  F_mat <- cbind(f1, f2, f3) + matrix(rnorm(n * 3, sd = 1), nrow = n)
  
  # assign class by max latent function
  friedman_y <- apply(F_mat, 1, which.max) - 1  # adjust for 0-based class indexing
  
  # get train and test sets deterministically, as data is already generated randomly (70/30 split)
  n_train <- floor(0.7 * n)
  friedman_X_train <- friedman_X[1:n_train, ]
  friedman_y_train <- friedman_y[1:n_train]
  friedman_X_test <- friedman_X[(n_train + 1):n, ]
  friedman_y_test <- friedman_y[(n_train + 1):n]
  
  # save data
  friedman_train_full <- cbind(friedman_X_train, friedman_y_train)
  friedman_test_full <- cbind(friedman_X_test, friedman_y_test)
  addWorksheet(wb, paste0("rep_", r, "_train"))
  addWorksheet(wb, paste0("rep_", r, "_test"))
  writeData(wb, sheet = paste0("rep_", r, "_train"), wave_train_full, colNames = TRUE)
  writeData(wb, sheet = paste0("rep_", r, "_test"), wave_test_full, colNames = TRUE)
  
  # ----- SOFT MPBART ------
  # run mcmc
  mcmc_output <- soft_mpbart(y_train = wave_y_train,
                             X_train = wave_X_train,
                             X_test = wave_X_test,
                             num_classes = 3,
                             num_burnin = 1500,
                             num_sim = 1500
  )
  pred_output <- soft_mpbart_predict(predictions_z = mcmc_output$mu_test_draws)
  
  # ----- RF -----
  # run rf
  #p <- ncol(wave_X_train)
  #mtry_grid <- unique(c(1, floor(sqrt(p)) - 1, floor(sqrt(p)), floor(sqrt(p)) + 1, 5, 10, p))
  #pred_output <- rf_multiclass_cv(X_train = wave_X_train, y_train = wave_y_train, X_test = wave_X_test, mtry_grid = mtry_grid)
  
  # ----- MPBART WITH HARD SPLITS -----
  
  # compute test error rate and brier score
  error <- test_error_rate(y_actual = wave_y_test, y_pred = pred_output$pred_y)
  brier_score <- brier_score_multiclass(y_actual = wave_y_test, y_prob = pred_output$post_probs)
  wave_error_rates[r] <- error
  wave_brier_scores[r] <- brier_score
}
  
# save data and/or output
saveWorkbook(wb, file = "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/friedman_data.xlsx", overwrite = TRUE)

addWorksheet(wb_output, "misclassification_rates")
addWorksheet(wb_output, "brier_scores")
writeData(wb_output, sheet = "misclassification_rates", x = wave_error_rates)
writeData(wb_output, sheet = "brier_scores", x = wave_brier_scores)

# (only write to one file, outcomment the other two)
saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/smpbart_friedman_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/rf_friedman_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/mpbart_friedman_output.xlsx", overwrite = TRUE)

