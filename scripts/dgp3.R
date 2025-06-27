# DGP 3: Gaussian process 
# Benchmark methods are Random Forest and MPBART with hard splits as implemented by Xu et al. (2024)

library(openxlsx)
source("soft_mpbart.R")
source("random_forest.R")
source("simulation_functions.R")

# -------- SET SEED AND PARAMS --------

seed <- 456
set.seed(seed)  # for reproducibility

n_train <- 500
n_test <- 1000
n_rep <- 10
method <- "rf"   # choose from "smpbart", "mpbart","rf"

# ------ GENERATE DATA -----

# generate replications
dgp3_data <- lapply(1:n_rep, function(i) generate_dgp3_data(n_train, n_test))

# write all generated data to excel
#write_data(n_rep = n_rep, all_data = dgp3_data, which_dgp = "dgp3", seed = seed)

# ------ RUN METHOD -------
run_method(method = method, sim_data = dgp3_data, which_dgp = "dgp3")
