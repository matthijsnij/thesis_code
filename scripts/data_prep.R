# Loading and preprocessing data sets for S-MPBART

# ---------------- FUNCTION FOR COVARIATE NORMALIZATION ----------------------

#'@description applies rank normalization to a matrix of doubles 
#'
#'@param X matrix of doubles
#'@return Normalized matrix of doubles
rank_normalize <- function(X) {
  apply(X, 2, function(col) {
    ranks <- rank(col, ties.method = "average")
    (ranks - 1) / (length(ranks) - 1)  # maps to [0, 1]
  })
}

# ------------ READ DATA ---------------

glass_data <- read.csv('C:\Users\matth\OneDrive\Bureaublad\msc_thesis\Data\glass\glass.data', header = FALSE)
glass_y <- glass_data[[ncol(glass_data)]]
glass_X <- as.matrix(glass_data[, 2:(ncol(glass_data)-1)])

vertebral_data <- read.csv('C:\Users\matth\OneDrive\Bureaublad\msc_thesis\Data\vertebral\data.txt', header = FALSE)
vertebral_y <- glass_data[[ncol(glass_data)]]
vertebral_X <- as.matrix(glass_data[, 1:(ncol(glass_data)-1)])

# ------------ PREPROCESS DATA -------------

# ---- GLASS -----
# clean the class labels such that they fall in the range [0,5]
# there is no class 4 in the data set
for (i in 1:length(glass_y)) {
  if (glass_y[i] < 4) {
    glass_y[i] <- glass_y[i] - 1
  } else {
    glass_y[i] <- glass_y[i] - 2
  }
}

# normalize covariates
glass_X_norm <- rank_normalize(glass_X)

# ---- VERTEBRAL -----
# change class labels to 0 = Hernia, 1 = Spondylolisthesis, 2 = Normal
for (i in 1:length(vertebral_y)) {
  if (vertebral_y[i] == "Hernia") {
    vertebral_y[i] <- 0
  } else if (vertebral_y[i] == "Spondylolisthesis") {
    vertebral_y[i] <- 1
  } else {
    vertebral_y[i] <- 2
  }
}

# normalize covariates
vertebral_X_norm <- rank_normalize(vertebral_X)


# --------- SAVE PREPROCESSED DATA SETS -------------

