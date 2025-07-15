# Loading and preprocessing data sets for S-MPBART

source("soft_mpbart.R")
source("data_application_functions.R")
library(AER)
library(mlogit)

# ------------ READ DATA ---------------

glass_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/glass/glass.data', header = TRUE)
glass_y <- glass_data[[ncol(glass_data)]]
glass_X <- as.matrix(glass_data[, 2:(ncol(glass_data)-1)])

vertebral_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/vertebral/data.txt', header = TRUE)
vertebral_y <- vertebral_data[[ncol(vertebral_data)]]
vertebral_X <- as.matrix(vertebral_data[, 1:(ncol(vertebral_data)-1)])

iris_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/iris/iris.data', header = TRUE)
iris_y <- iris_data[[ncol(iris_data)]]
iris_X <- as.matrix(iris_data[, 1:(ncol(iris_data)-1)])

wine_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/wine/wine.data', header = TRUE, sep = ",")
wine_y <- wine_data[[1]] 
wine_X <- as.matrix(wine_data[, 2:ncol(wine_data)])

vehicle_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/vehicle/vehicle_data.txt', header = TRUE)
vehicle_y <- vehicle_data[[ncol(vehicle_data)]]
vehicle_X <- as.matrix(vehicle_data[, 1:(ncol(vehicle_data)-1)])

vowel_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/vowel/vowel-context.data', header = TRUE)
vowel_y <- vowel_data[[ncol(vowel_data)]]
vowel_X <- as.matrix(vowel_data[, 1:(ncol(vowel_data)-1)])

waveform_data <- read.csv('C:/Users/matth/OneDrive/Bureaublad/msc_thesis/Data/waveform/waveform.data', header = TRUE, sep = ",")
waveform_y <- waveform_data[[ncol(waveform_data)]]
waveform_X <- as.matrix(waveform_data[, 1:(ncol(waveform_data)-1)])

data("TravelMode")
travel_data <- prep_travelmode(TravelMode)
travel_y <- travel_data$y
travel_X <- travel_data$X

data("Fishing", package = "mlogit")
fishing_y <- Fishing[[1]]
fishing_X <- as.matrix(Fishing[, 2:ncol(Fishing)])



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

# -- IRIS --
# change class labels to 0 = Iris-setosa, 1 = Iris-versicolour, 2 = Iris-virginica
for (i in 1:length(iris_y)) {
  if (iris_y[i] == "Iris-setosa") {
    iris_y[i] <- 0
  } else if (iris_y[i] == "Iris-versicolor") {
    iris_y[i] <- 1
  } else {
    iris_y[i] <- 2
  }
}

# normalize covariates
iris_X_norm <- rank_normalize(iris_X)

# -- WINE --
# 0 based
wine_y <- wine_y - 1
# normalize covariates
wine_X_norm <- rank_normalize(wine_X)

# -- VEHICLE --
# 0-based
for (i in 1:length(vehicle_y)) {
  if (vehicle_y[i] == "opel") {
    vehicle_y[i] <- 0
  } else if (vehicle_y[i] == "saab") {
    vehicle_y[i] <- 1
  } else if (vehicle_y[i] == "bus") {
    vehicle_y[i] <- 2
  } else {
    vehicle_y[i] <- 3
  }
}

# normalize covariates
vehicle_y <- as.numeric(vehicle_y)
vehicle_X_norm <- rank_normalize(vehicle_X)

# -- VOWEL --

# normalize covariates
vowel_X_norm <- rank_normalize(vowel_X)

# -- WAVEFORM --
waveform_X_norm <- rank_normalize(waveform_X)

# -- TRAVEL MODE --
travel_X_norm <- rank_normalize(travel_X)

# -- FISHING --
# change class labels to 0 = beach, 1 = pier, 2 = boat, 3 = charter
fishing_y <- as.integer(fishing_y) - 1

fishing_X_norm <- rank_normalize(fishing_X)


# --------- SAVE PREPROCESSED DATA SETS -------------
write.table(as.data.frame(cbind(glass_X_norm, glass_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/glass_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(vertebral_X_norm, vertebral_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/vertebral_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(iris_X_norm, iris_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/iris_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(wine_X_norm, wine_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/wine_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(vehicle_X_norm, vehicle_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/vehicle_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(vowel_X_norm, vowel_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/vowel_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(waveform_X_norm, waveform_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/waveform_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(travel_X_norm, travel_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/travel_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(cbind(fishing_X_norm, fishing_y)), "C:/Users/matth/OneDrive/Bureaublad/msc_thesis/thesis_code/data/fishing_preprocessed.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)






