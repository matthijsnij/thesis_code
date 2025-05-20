# Simulation study for S-MPBART

library(mlbench)
source("soft_mpbart.R")

# ------- WAVEFORM RECOGNITION -------

# perform 20 replications and save the error rates
reps <- 20
wave_error_rates <- rep(0, reps)

for (r in 1:reps) {
  # generate waveform data
  train_wave = mlbench.waveform(300)
  test_wave = mlbench.waveform(500)
  
  # normalize covariates and define input for S-MPBART
  wave_X_train <- rank_normalize(as.matrix(train_wave$x))
  wave_y_train <- as.numeric(train_wave$classes) - 1 # convert to 0-based class labels
  wave_X_test  <- rank_normalize(as.matrix(test_wave$x))
  wave_y_test <- as.numeric(test_wave$classes) - 1 # convert to 0-based class labels
  
  # run mcmc
  mcmc_output <- soft_mpbart(y_train = wave_y_train,
                             X_train = wave_X_train,
                             X_test = wave_X_test,
                             K = 2)
  
  # predict
  pred_output <- soft_mpbart_predict(predictions_z = mcmc_output$mu_test_draws)
  
  # compute test error rate
  error <- test_error_rate(y_actual = wave_y_test, y_pred = pred_output$pred_y)
  wave_error_rates[r] <- error
}


