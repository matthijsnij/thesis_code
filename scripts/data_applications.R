# Data application for S-MPBART

source("soft_mpbart.R")

# --------- READ DATA ----------

# read preprocessed, ready-to-use data
glass <- read.csv("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/glass_preprocessed.csv", header = TRUE)
glass_y <- glass[[ncol(glass)]]
glass_X <- glass[, 1:ncol(glass) - 1]

vertebral <- read.csv("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/vertebral_preprocessed.csv", header = TRUE)
vertebral_y <- vertebral[[ncol(vertebral)]]
vertebral_X <- vertebral[, 1:ncol(vertebral) - 1]

# ------------ CREATE TRAIN AND TEST SETS

# generate train and test set folds
num_folds <- 5
folds <- createFolds(glass_y, k = num_folds, list = TRUE, returnTrain = TRUE)

# split the data based on these indices
glass_y_train <- glass_y[folds[[1]]]
glass_X_train <- glass_X_norm[folds[[1]], ]
glass_y_test <- glass_y[-folds[[1]]]
glass_X_test <- glass_X_norm[-folds[[1]], ]