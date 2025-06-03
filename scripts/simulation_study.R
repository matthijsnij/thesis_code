# Simulation study for S-MPBART, with random forest and MBART with hard splits as benchmarks.

library(mlbench)
library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")

# ------- WAVEFORM RECOGNITION -------

# perform 10 replications and save the error rates
reps <- 10
wave_error_rates <- rep(0, reps)
wave_brier_scores <- rep(0, reps)

# set seed for reproduceability
set.seed(123)

# create workbook to store output
wb_output <- createWorkbook()

for (r in 1:reps) {
  # generate waveform data
  train_wave = mlbench.waveform(300)
  test_wave = mlbench.waveform(500)
  
  # normalize covariates and define input
  wave_X_train <- rank_normalize(as.matrix(train_wave$x))
  wave_y_train <- as.numeric(train_wave$classes) - 1 # convert to 0-based class labels
  wave_X_test  <- rank_normalize(as.matrix(test_wave$x))
  wave_y_test <- as.numeric(test_wave$classes) - 1 # convert to 0-based class labels
  
  # define column names
  p <- ncol(wave_X_train)
  colnames(wave_X_train) <- paste0("x", 1:p)
  colnames(wave_X_test) <- paste0("x", 1:p)
  
  #wave_train_full <- cbind(wave_X_train, y = wave_y_train)
  #wave_test_full <- cbind(wave_X_test, y = wave_y_test)
  
  # save data
  #addWorksheet(wb, paste0("rep_", r, "_train"))
  #addWorksheet(wb, paste0("rep_", r, "_test"))
  #writeData(wb, sheet = paste0("rep_", r, "_train"), wave_train_full, colNames = TRUE)
  #writeData(wb, sheet = paste0("rep_", r, "_test"), wave_test_full, colNames = TRUE)
  
  # ----- SOFT MPBART ------
  # run mcmc
  #mcmc_output <- soft_mpbart(y_train = wave_y_train,
                             #X_train = wave_X_train,
                             #X_test = wave_X_test,
                             #num_classes = 3,
                             #num_burnin = 1500,
                             #num_sim = 1500
  #)
  #pred_output <- soft_mpbart_predict(predictions_z = mcmc_output$mu_test_draws)
  
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
#saveWorkbook(wb, file = "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/waveform_data.xlsx", overwrite = TRUE)

addWorksheet(wb_output, "misclassification_rates")
addWorksheet(wb_output, "brier_scores")
writeData(wb_output, sheet = "misclassification_rates", x = wave_error_rates)
writeData(wb_output, sheet = "brier_scores", x = wave_brier_scores)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/smpbart_waveform_output.xlsx", overwrite = TRUE)
saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/rf_waveform_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/mpbart_waveform_output.xlsx", overwrite = TRUE)


