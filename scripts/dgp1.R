# DGP 1: Piecewise sinusoid + discontinuity 

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
dgp1_data <- lapply(1:n_rep, function(i) generate_dgp1_data(n_train, n_test))

# write all generated data to excel
#write_data(n_rep = n_rep, all_data = dgp1_data, which_dgp = "dgp1", seed = seed)

# ------ RUN METHOD -------
run_method(method = method, sim_data = dgp1_data, which_dgp = "dgp1")




