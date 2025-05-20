# Loading and preprocessing data sets for S-MPBART

source("soft_mpbart.R")

# ------------ READ DATA ---------------

glass_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/glass/glass.data', header = TRUE)
glass_y <- glass_data[[ncol(glass_data)]]
glass_X <- as.matrix(glass_data[, 2:(ncol(glass_data)-1)])

vertebral_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/vertebral/data.txt', header = TRUE)
vertebral_y <- vertebral_data[[ncol(vertebral_data)]]
vertebral_X <- as.matrix(vertebral_data[, 1:(ncol(vertebral_data)-1)])

# ------------ PREPROCESS DATA -------------

# -- GLASS --
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

# -- VERTEBRAL --
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
write.table(as.data.frame(cbind(glass_X_norm, glass_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/glass_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(vertebral_X_norm, vertebral_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/vertebral_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)







