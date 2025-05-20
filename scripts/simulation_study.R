# Simulation study for S-MPBART

library(mlbench)
source("soft_mpbart.R")

# ------- WAVEFORM RECOGNITION -------
train_wave = mlbench.waveform(300)
test_wave = mlbench.waveform(500)

wave_X_train <- rank_normalize(as.matrix(train_wave$x))
wave_y_train <- as.numeric(train_wave$classes) - 1 # convert to 0-based class labels
wave_X_test  <- rank_normalize(as.matrix(test_wave$x))
wave_y_test <- as.numeric(test_wave$classes) - 1 # convert to 0-based class labels

mcmc_output <- soft_mpbart(y_train = wave_y_train,
                           X_train = wave_X_train,
                           X_test = wave_X_test,
                           K = 2,
                           seed = 123)

