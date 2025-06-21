# Run S-MPBART on real data set
# Benchmark methods are Random Forest and MPBART with hard splits as implemented by Xu et al. (2024)

library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")
source("data_application_functions.R")

# ---- READ DATA -----

# specify which data set to run
dataset <- "iris"   # choose from "glass", "vertebral", "iris"
data <- read_data(dataset = dataset)

# specify which method to run
method <- "smpbart"  # choose from "smpbart", "mpbart","rf"

# ---- RUN METHOD ON TRAIN/TEST SPLITS OF REAL DATA -----
seed <- 123
run_method(method = method, data = data, which_dataset = dataset, seed = seed)

