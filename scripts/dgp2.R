# DGP 2: Friedman 
# Benchmark methods are Random Forest and MPBART with hard splits as implemented by Xu et al. (2024)

library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")
source("simulation_functions.R")

# -------- SET SEED AND PARAMS --------

seed <- 123
set.seed(seed)  # for reproducibility

n_train <- 500
n_test <- 1000
n_rep <- 10
method <- "rf"   # choose from "smpbart", "mpbart","rf"

# ------ GENERATE DATA -----

# generate replications
n_extra_noise <- 0 # number of extra noise predictors to include
dgp2_data <- lapply(1:n_rep, function(i) generate_dgp2_data(n_train, n_test, p = (10 + n_extra_noise)))

# write all generated data to excel
write_data(n_rep = n_rep, all_data = dgp2_data, which_dgp = "dgp2", seed = seed)

# ------ RUN METHOD -------
run_method(method = method, sim_data = dgp2_data, which_dgp = "dgp2")