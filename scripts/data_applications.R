# Data application for S-MPBART

# --------- READ DATA ----------

# read preprocessed, ready-to-use data


# ------------ CREATE TRAIN AND TEST SETS

# generate train and test set folds
num_folds <- 5
folds <- createFolds(glass_y, k = num_folds, list = TRUE, returnTrain = TRUE)

# split the data based on these indices
glass_y_train <- glass_y[folds[[1]]]
glass_X_train <- glass_X_norm[folds[[1]], ]
glass_y_test <- glass_y[-folds[[1]]]
glass_X_test <- glass_X_norm[-folds[[1]], ]