# Data application for S-MPBART - Forensic glass data set
# Benchmark methods are Random Forest and MPBART with hard splits as implemented by Xu et al. (2024)

library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")

# --------- READ DATA ----------

# read preprocessed, ready-to-use data
glass <- read.csv("C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/glass_preprocessed.csv", header = TRUE)
glass_y <- glass[[ncol(glass)]]
glass_X <- as.matrix(glass[, 1:(ncol(glass) - 1)])

# ------------ CREATE TRAIN AND TEST SETS ------------

# generate train and test set folds
set.seed(123)
num_folds <- 5
folds <- createFolds(glass_y, k = num_folds, list = TRUE, returnTrain = TRUE)

# to store error rates and brier scores
glass_error_rates <- rep(0, num_folds)
glass_brier_scores <- rep(0, num_folds)

# create workbook to store output
wb_output <- createWorkbook() 

# ------- ANALYSIS, OUTCOMMENT THE 2 METHODS YOU DO NOT WANT TO RUN -------
for (i in seq_along(folds)) {
  # split the data based on fold
  glass_y_train <- glass_y[folds[[i]]]
  glass_X_train <- as.matrix(glass_X[folds[[i]], ])
  glass_y_test <- glass_y[-folds[[i]]]
  glass_X_test <- as.matrix(glass_X[-folds[[i]], ])
  
  # ----- S-MPBART -----
  # run mcmc
  mcmc_output <- soft_mpbart(y_train = glass_y_train,
                             X_train = glass_X_train,
                             X_test = glass_X_test,
                             num_classes = 6,
                             num_burnin = 1500,
                             num_sim = 1500
  )
  pred_output <- soft_mpbart_predict(predictions_z = mcmc_output$mu_test_draws)
  
  # ----- RANDOM FOREST -----
  # run rf
  #p <- ncol(glass_X_train)
  #mtry_grid <- 1:p
  #pred_output <- rf_multiclass_cv(X_train = glass_X_train, y_train = glass_y_train, X_test = glass_X_test, mtry_grid = mtry_grid)
  
  # ----- MPBART -----
  
  # compute and store misclassification error and brier score
  error <- test_error_rate(y_actual = glass_y_test, y_pred = pred_output$pred_y)
  brier_score <- brier_score_multiclass(y_actual = glass_y_test, y_prob = pred_output$post_probs)
  glass_error_rates[i] <- error
  glass_brier_scores[i] <- brier_score
}

# save output 
addWorksheet(wb_output, "misclassification_rates")
addWorksheet(wb_output, "brier_scores")
writeData(wb_output, sheet = "misclassification_rates", x = glass_error_rates)
writeData(wb_output, sheet = "brier_scores", x = glass_brier_scores)

# (only write to one file, outcomment the other two, depending on which method is used)
saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/smpbart_glass_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/rf_glass_output.xlsx", overwrite = TRUE)
#saveWorkbook(wb_output, "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/output/mpbart_glass_output.xlsx", overwrite = TRUE)

