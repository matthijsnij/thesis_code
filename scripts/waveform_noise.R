# Simulation study for S-MPBART, with random forest and MBART with hard splits as benchmarks
# Waveform recognition with noise variables.

library(mlbench)
library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")

# ------- WAVEFORM RECOGNITION WITH NOISE VARIABLES -------

# perform 10 replications and save the error rates
reps <- 10
wave_error_rates <- rep(0, reps)
wave_brier_scores <- rep(0, reps)
num_noise <- 60 # number of noise variables to add

# set seed for reproduceability
set.seed(123)

# create workbooks
wb <- createWorkbook() # for data
wb_output <- createWorkbook() # for output

for (r in 1:reps) {
  # generate waveform data
  train_wave = mlbench.waveform(300)
  test_wave = mlbench.waveform(500)
  
  wave_X_train_real <- as.matrix(train_wave$x)
  wave_X_test_real <- as.matrix(test_wave$x)
  
  wave_X_train_noise <- matrix(rnorm(nrow(wave_X_train_real) * n_noise), nrow = nrow(wave_X_train_real), ncol = n_noise)
  wave_X_test_noise <- matrix(rnorm(nrow(wave_X_test_real) * n_noise), nrow = nrow(wave_X_test_real), ncol = n_noise)
  
  # combine and normalize
  wave_X_train <- rank_normalize(cbind(wave_X_train_real, wave_X_train_noise))
  wave_X_test  <- rank_normalize(cbind(wave_X_test_real, wave_X_test_noise))
  
  wave_y_train <- as.numeric(train_wave$classes) - 1 # convert to 0-based class labels
  wave_y_test <- as.numeric(test_wave$classes) - 1 # convert to 0-based class labels
  
  # define column names
  p <- ncol(wave_X_train)
  colnames(wave_X_train) <- paste0("x", 1:p)
  colnames(wave_X_test) <- paste0("x", 1:p)
  
  wave_train_full <- cbind(wave_X_train, y = wave_y_train)
  wave_test_full <- cbind(wave_X_test, y = wave_y_test)
  
  # save data
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
saveWorkbook(wb, file = "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/waveform_noise_data.xlsx", overwrite = TRUE)

addWorksheet(wb_output, "misclassification_rates")
addWorksheet(wb_output, "brier_scores")
writeData(wb_output, sheet = "misclassification_rates", x = wave_error_rates)
writeData(wb_output, sheet = "brier_scores", x = wave_brier_scores)

# (only write to one file, outcomment the other two)
saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/smpbart_waveform_noise_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/rf_waveform_noise_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/mpbart_waveform_noise_output.xlsx", overwrite = TRUE)
